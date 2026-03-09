import streamlit as st
import os
import subprocess
import re
import numpy as np
import py3Dmol
import streamlit.components.v1 as components

# --- ページ設定 ---
st.set_page_config(page_title="XTB-CREST Web Console v8.3", layout="wide")

st.title("🧪 XTB-CREST Web Console v8.3")
st.write("Multi-Format Support: XYZ, MOL, PDB, SDF")

# --- 状態リセット ---
def reset_all():
    for f in ["xtb.trj", "input.xyz", "crest_conformers.xyz", "crestopt.log", "xtbopt.xyz"]:
        if os.path.exists(f): os.remove(f)
    if "current_file" in st.session_state:
        del st.session_state.current_file
    st.rerun()

# --- サイドバー ---
with st.sidebar:
    st.header("⚙️ Settings")
    comp_name = st.text_input("Compound Name", value="", placeholder="Enter name")
    level = st.radio("Calculation Level", ["GFN2-xTB (High Accuracy)", "GFN-FF (Ultra Robust)"])
    solvent = st.selectbox("Solvent (ALPB)", ["methanol", "water", "chcl3", "benzene", "none"])
    cores = st.slider("CPU Cores", 1, 12, 4)
    
    if st.button("🗑️ Clear All & Reset"):
        reset_all()

# --- メインパネル ---
col1, col2 = st.columns([1, 1])

with col1:
    st.subheader("📥 Export Results")
    # v8.3: 対応形式を拡大
    uploaded_file = st.file_uploader("Upload Structure File", type=["xyz", "mol", "pdb", "sdf"])

    if uploaded_file:
        file_ext = uploaded_file.name.split('.')[-1].lower()
        
        if "current_file" not in st.session_state or st.session_state.current_file != uploaded_file.name:
            reset_all() # 新しいファイルが来たらリセット
            st.session_state.current_file = uploaded_file.name

        # 一旦そのまま保存
        temp_input = f"raw_input.{file_ext}"
        with open(temp_input, "wb") as f:
            f.write(uploaded_file.getbuffer())

        # --- v8.3: フォーマット変換ロジック ---
        if file_ext != "xyz":
            st.info(f"Converting {file_ext.upper()} to XYZ...")
            # xtbを使用してXYZに変換
            subprocess.run(["xtb", temp_input, "--steps", "0"], capture_output=True)
            if os.path.exists("xtbopt.xyz"):
                os.replace("xtbopt.xyz", "input.xyz")
            else:
                # 変換に失敗した場合の予備（単純リネームは危険だが一応）
                st.warning("Direct conversion failed. Attempting standard parse.")
        else:
            # XYZの場合は洗浄処理を適用
            os.replace(temp_input, "input.xyz")

        # 洗浄とチェック
        if os.path.exists("input.xyz"):
            with open("input.xyz", "r") as f:
                content = f.read()
            # v8.2の洗浄ロジックを適用
            lines = [l.strip() for l in content.splitlines() if l.strip()]
            try:
                num_atoms = int(lines[0])
                st.success(f"✅ Structure loaded: {num_atoms} atoms.")
            except:
                st.error("Format Error: Could not determine atom count.")

        if st.button("🚀 Run Analysis", type="primary"):
            with st.spinner("Calculating..."):
                method = "--gfn2" if level == "GFN2-xTB (High Accuracy)" else "--gfnff"
                cmd = ["crest", "input.xyz", method, "-T", str(cores), "--quick"]
                if solvent != "none": cmd.extend(["--alpb", solvent])
                
                subprocess.run(cmd)
                
                if os.path.exists("crest_conformers.xyz"):
                    os.replace("crest_conformers.xyz", "xtb.trj")
                    st.success("Analysis Complete!")
                    st.rerun()
                else:
                    st.error("Calculation failed. Try GFN-FF level.")

    # 解析データの読み込み
    frames = []
    if os.path.exists("xtb.trj"):
        try:
            with open("xtb.trj", 'r') as f:
                lines = f.readlines()
            if lines:
                at_num = int(lines[0].strip())
                chunk_size = at_num + 2
                for i in range(0, len(lines), chunk_size):
                    chunk = lines[i:i+chunk_size]
                    if len(chunk) < chunk_size: break
                    frames.append({"energy": 0.0, "coords": chunk}) # 簡略化
        except: pass

    is_ready = True if (frames and comp_name) else False
    st.download_button("Download SDF", data="", disabled=not is_ready)

with col2:
    st.subheader("🔍 3D Viewer")
    if frames:
        rank = st.slider("Select Conformer", 1, len(frames), 1) if len(frames) > 1 else 1
        xyz_data = "".join(frames[rank-1]["coords"])
        view = py3Dmol.view(width=450, height=450)
        view.addModel(xyz_data, 'xyz')
        view.setStyle({'stick': {'radius': 0.15}, 'sphere': {'scale': 0.25}})
        view.zoomTo()
        components.html(view._make_html(), height=460)
    else:
        st.info("Upload any supported format (MOL, PDB, XYZ) to start.")