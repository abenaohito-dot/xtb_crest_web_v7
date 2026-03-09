import streamlit as st
import pandas as pd
import numpy as np
import re
import os
import py3Dmol
import streamlit.components.v1 as components

# --- ページ設定 ---
st.set_page_config(page_title="XTB-CREST Web Console v7.4", layout="wide")

st.title("🧪 XTB-CREST Web Console v7.4")
st.write("Natural Products Conformer Ensemble Analyzer")

# --- サイドバー：設定 ---
with st.sidebar:
    st.header("⚙️ Settings")
    comp_name = st.text_input("Compound Name", value="", placeholder="Enter name to unlock downloads")
    
    mode = st.radio("Execution Mode", ["xTB-Opt", "xTB-MD", "CREST-MTD"])
    solvent = st.selectbox("Solvent (ALPB)", ["methanol", "water", "chcl3", "benzene", "none"])
    
    st.divider()
    st.header("⚖️ Energy Thresholds (kcal/mol)")
    low_thresh = st.number_input("Threshold 1 (Low)", value=3.0, step=0.5)
    high_thresh = st.number_input("Threshold 2 (High)", value=10.0, step=1.0)

# --- メインパネルのレイアウト ---
col1, col2 = st.columns([1, 1])

with col1:
    st.subheader("📥 Export Results")
    uploaded_file = st.file_uploader("Upload XYZ File", type=["xyz"])
    
    # セッション状態の管理：新しいファイルが来たら古い結果をリセット
    if uploaded_file:
        # アップロードされたファイル名が以前と違う場合、古いtrjを削除
        if "last_uploaded" not in st.session_state or st.session_state.last_uploaded != uploaded_file.name:
            if os.path.exists("xtb.trj"):
                os.remove("xtb.trj")
            st.session_state.last_uploaded = uploaded_file.name

    trj_exists = os.path.exists("xtb.trj")
    is_ready = True if (trj_exists and comp_name) else False

    frames = []
    num_atoms = 0
    min_e = 0.0

    if trj_exists:
        with open("xtb.trj", 'r') as f:
            lines = f.readlines()
        if len(lines) > 0:
            try:
                num_atoms = int(lines[0].strip())
                chunk_size = num_atoms + 2
                for i in range(0, len(lines), chunk_size):
                    chunk = lines[i:i+chunk_size]
                    if len(chunk) < chunk_size: break
                    e_match = re.search(r"energy:\s+([-+]?\d+\.\d+)", chunk[1])
                    frames.append({"energy": float(e_match.group(1)) if e_match else 0.0, "coords": chunk})
                frames.sort(key=lambda x: x["energy"])
                if frames: min_e = frames[0]["energy"]
            except:
                st.error("Error parsing the result file.")

    # SDF生成
    def get_sdf(threshold):
        if not is_ready or not frames: return ""
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

    # ボタン表示（常に表示、条件未達でグレーアウト）
    s1_data = get_sdf(low_thresh)
    st.download_button(
        label=f"Download {comp_name if comp_name else '---'}_{low_thresh}kcal.sdf",
        data=s1_data,
        file_name=f"{comp_name}_{low_thresh}kcal.sdf",
        disabled=not is_ready,
        key="btn1"
    )
    
    s2_data = get_sdf(high_thresh)
    st.download_button(
        label=f"Download {comp_name if comp_name else '---'}_{high_thresh}kcal.sdf",
        data=s2_data,
        file_name=f"{comp_name}_{high_thresh}kcal.sdf",
        disabled=not is_ready,
        key="btn2"
    )

    if not trj_exists:
        st.info("💡 1. Upload XYZ  2. Run xTB/Launch Test")
        if uploaded_file and st.button("🚀 Launch Test Mode"):
            # アップロードされたファイルを元にダミーのtrjを作成
            dummy_content = f"{num_atoms if num_atoms > 0 else '...'}\n energy: -1.0\n" + uploaded_file.getvalue().decode()
            with open("xtb.trj", "w") as f:
                f.write(uploaded_file.getvalue().decode())
            st.rerun()

with col2:
    st.subheader("🔍 3D Viewer")
    if is_ready and frames:
        # スライダーの最小/最大値不整合を回避
        if len(frames) > 1:
            rank = st.slider("Select Conformer", 1, len(frames), 1)
        else:
            st.text("Conformer: 1 / 1 (Single structure)")
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
        st.markdown("""
            <div style="width:450px; height:450px; background-color:#262730; border-radius:10px; display:flex; align-items:center; justify-content:center; color:#555; border: 1px dashed #444;">
                Viewer Locked (Check Input & Name)
            </div>
        """, unsafe_allow_html=True)