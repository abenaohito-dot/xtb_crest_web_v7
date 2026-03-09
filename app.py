import streamlit as st
import pandas as pd
import numpy as np
import re
import os
import subprocess
import py3Dmol
import streamlit.components.v1 as components

# --- ページ設定 ---
st.set_page_config(page_title="XTB-CREST Web Console v7.7", layout="wide")

st.title("🧪 XTB-CREST Web Console v7.7")
st.write("Natural Products Conformer Ensemble Analyzer")

# --- 状態リセット用関数 ---
def reset_all():
    if os.path.exists("xtb.trj"):
        os.remove("xtb.trj")
    if "current_file" in st.session_state:
        del st.session_state.current_file
    st.rerun()

# --- サイドバー：設定 ---
with st.sidebar:
    st.header("⚙️ Settings")
    comp_name = st.text_input("Compound Name", value="", placeholder="Enter name to unlock downloads")
    
    calc_mode = st.radio("Calculation Method", ["CREST (Conformer Search)", "xTB (Optimization Only)"])
    solvent = st.selectbox("Solvent (ALPB)", ["methanol", "water", "chcl3", "benzene", "none"])
    cores = st.slider("CPU Cores", 1, 12, 4)
    
    st.divider()
    st.header("⚖️ Energy Thresholds (kcal/mol)")
    low_thresh = st.number_input("Threshold 1 (Low)", value=3.0, step=0.5)
    high_thresh = st.number_input("Threshold 2 (High)", value=10.0, step=1.0)
    
    st.divider()
    if st.button("🗑️ Clear All & Reset"):
        reset_all()

# --- メインパネル ---
col1, col2 = st.columns([1, 1])

with col1:
    st.subheader("📥 Export Results")
    uploaded_file = st.file_uploader("Upload XYZ File", type=["xyz"])
    
    # 新しいファイルが来たら古い結果を削除
    if uploaded_file:
        if "current_file" not in st.session_state or st.session_state.current_file != uploaded_file.name:
            if os.path.exists("xtb.trj"):
                os.remove("xtb.trj")
            st.session_state.current_file = uploaded_file.name

    trj_file = "xtb.trj"
    trj_exists = os.path.exists(trj_file)
    is_ready = True if (trj_exists and comp_name) else False

    # --- 解析実行エリア（常に表示） ---
    if uploaded_file:
        st.write("---")
        st.markdown("### 🚀 Calculation Control")
        b_col1, b_col2 = st.columns(2)
        
        with b_col1:
            # trjがあっても「再計算」できるように常に表示
            btn_label = "🔄 Re-run Real Analysis" if trj_exists else "🚀 Run Real CREST/xTB"
            if st.button(btn_label, type="primary"):
                with st.spinner("Calculating..."):
                    # 古いファイルを一旦消す
                    if os.path.exists("xtb.trj"): os.remove("xtb.trj")
                    
                    with open("input.xyz", "wb") as f:
                        f.write(uploaded_file.getbuffer())
                    
                    if calc_mode == "CREST (Conformer Search)":
                        cmd = ["crest", "input.xyz", "--gfn2", "-T", str(cores)]
                        if solvent != "none": cmd.extend(["--alpb", solvent])
                        subprocess.run(cmd)
                        if os.path.exists("crest_conformers.xyz"):
                            os.replace("crest_conformers.xyz", "xtb.trj")
                    else:
                        cmd = ["xtb", "input.xyz", "--opt", "--gfn2"]
                        if solvent != "none": cmd.extend(["--alpb", solvent])
                        subprocess.run(cmd)
                        if os.path.exists("xtbopt.xyz"):
                            os.replace("xtbopt.xyz", "xtb.trj")
                    st.rerun()

        with b_col2:
            if st.button("🛠️ Launch Test (Dummy)"):
                with open("xtb.trj", "wb") as f:
                    f.write(uploaded_file.getbuffer())
                st.rerun()
        st.write("---")

    # データ読み込み・SDF生成ロジック（以下、前回と同じ）
    frames = []
    num_atoms = 0
    min_e = 0.0
    if trj_exists:
        with open(trj_file, 'r') as f:
            lines = f.readlines()
        if lines:
            num_atoms = int(lines[0].strip())
            chunk_size = num_atoms + 2
            for i in range(0, len(lines), chunk_size):
                chunk = lines[i:i+chunk_size]
                if len(chunk) < chunk_size: break
                e_match = re.search(r"energy:\s+([-+]?\d+\.\d+)", chunk[1])
                frames.append({"energy": float(e_match.group(1)) if e_match else 0.0, "coords": chunk})
            frames.sort(key=lambda x: x["energy"])
            if frames: min_e = frames[0]["energy"]

    def make_sdf(threshold):
        if not is_ready: return ""
        sdf = ""
        for i, f in enumerate(frames, start=1):
            rel_e = (f["energy"] - min_e) * 627.509
            if rel_e > threshold: break
            sdf += f"{comp_name}_{i}\nThreshold: {threshold} kcal\n\n"
            sdf += f"{num_atoms:>3}  0  0  0  0  0  0  0  0  0999 V2000\n"
            for line in f["coords"][2:]:
                p = line.split()
                if len(p) >= 4:
                    sdf += f"{float(p[1]):>10.4f}{float(p[2]):>10.4f}{float(p[3]):>10.4f} {p[0]:<3} 0  0  0  0  0\n"
            sdf += "M  END\n> <ENERGY_KCAL>\n{:.4f}\n\n$$$$\n".format(rel_e)
        return sdf

    s1 = make_sdf(low_thresh)
    st.download_button(f"Download {comp_name if comp_name else '---'}_{low_thresh}kcal.sdf", data=s1, file_name=f"{comp_name}_{low_thresh}kcal.sdf", disabled=not is_ready)
    
    s2 = make_sdf(high_thresh)
    st.download_button(f"Download {comp_name if comp_name else '---'}_{high_thresh}kcal.sdf", data=s2, file_name=f"{comp_name}_{high_thresh}kcal.sdf", disabled=not is_ready)

with col2:
    st.subheader("🔍 3D Viewer")
    if is_ready and frames:
        if len(frames) > 1:
            rank = st.slider("Select Conformer", 1, len(frames), 1)
        else:
            st.info("Single Result Displayed")
            rank = 1
        
        f_view = frames[rank-1]
        xyz_data = "".join(f_view["coords"])
        view = py3Dmol.view(width=450, height=450)
        view.addModel(xyz_data, 'xyz')
        view.setStyle({'stick': {'radius': 0.15}, 'sphere': {'scale': 0.25}})
        view.zoomTo()
        components.html(view._make_html(), height=460)
        st.write(f"Relative Energy: **{(f_view['energy'] - min_e) * 627.509:.4f} kcal/mol**")
    else:
        st.markdown('<div style="width:450px; height:450px; background-color:#262730; border-radius:10px; display:flex; align-items:center; justify-content:center; color:#555; border: 1px dashed #444;">Viewer Locked</div>', unsafe_allow_html=True)