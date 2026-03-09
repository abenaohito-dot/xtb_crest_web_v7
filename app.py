import streamlit as st
import pandas as pd
import numpy as np
import re
import os
import subprocess
import py3Dmol
import streamlit.components.v1 as components

# --- ページ設定 ---
st.set_page_config(page_title="XTB-CREST Web Console v8.5", layout="wide")

st.title("🧪 XTB-CREST Web Console v8.5")
st.write("Robust Reliability Mode")

# --- 状態リセット関数 ---
def reset_all():
    files = ["xtb.trj", "input.xyz", "crest_conformers.xyz", "xtbopt.xyz", "crestopt.log", "raw_input.*"]
    for f in files:
        if os.path.exists(f):
            try: os.remove(f)
            except: pass
    st.rerun()

# --- サイドバー ---
with st.sidebar:
    st.header("⚙️ Settings")
    comp_name = st.text_input("Compound Name", value="", placeholder="Required for download")
    calc_mode = st.radio("Method", ["CREST (Conformer Search)", "xTB (Optimization)"])
    solvent = st.selectbox("Solvent", ["methanol", "water", "chcl3", "benzene", "none"])
    cores = st.slider("CPU Cores", 1, 12, 4)
    
    st.divider()
    st.header("⚖️ Thresholds (kcal/mol)")
    low_thresh = st.number_input("Threshold 1", value=3.0, step=0.5)
    high_thresh = st.number_input("Threshold 2", value=10.0, step=1.0)
    
    if st.button("🗑️ Reset Application"):
        reset_all()

# --- メインパネル ---
col1, col2 = st.columns([1, 1])

with col1:
    st.subheader("📥 Input & Control")
    uploaded_file = st.file_uploader("Upload File (XYZ/MOL/PDB)", type=["xyz", "mol", "pdb", "sdf"])
    
    if uploaded_file:
        # ファイル保存とXYZ変換
        ext = uploaded_file.name.split('.')[-1].lower()
        with open(f"raw_input.{ext}", "wb") as f:
            f.write(uploaded_file.getbuffer())
        
        # フォーマット変換
        if ext != "xyz":
            subprocess.run(["xtb", f"raw_input.{ext}", "--steps", "0"], capture_output=True)
            if os.path.exists("xtbopt.xyz"): os.replace("xtbopt.xyz", "input.xyz")
        else:
            os.replace(f"raw_input.{ext}", "input.xyz")

        # 実行ボタン
        if st.button("🚀 Start Calculation", type="primary"):
            # 以前の結果を削除
            if os.path.exists("xtb.trj"): os.remove("xtb.trj")
            
            with st.status("Executing calculation...", expanded=True) as status:
                st.write("Step 1: GFN-FF Pre-optimization...")
                subprocess.run(["xtb", "input.xyz", "--gfnff", "--opt", "-T", str(cores)])
                
                target = "xtbopt.xyz" if os.path.exists("xtbopt.xyz") else "input.xyz"
                
                st.write(f"Step 2: Running {calc_mode}...")
                if "CREST" in calc_mode:
                    cmd = ["crest", target, "--gfn2", "-T", str(cores), "--quick"]
                    if solvent != "none": cmd.extend(["--alpb", solvent])
                    res = subprocess.run(cmd, capture_output=True, text=True)
                    if os.path.exists("crest_conformers.xyz"):
                        os.replace("crest_conformers.xyz", "xtb.trj")
                else:
                    cmd = ["xtb", target, "--opt", "--gfn2", "-T", str(cores)]
                    if solvent != "none": cmd.extend(["--alpb", solvent])
                    res = subprocess.run(cmd, capture_output=True, text=True)
                    if os.path.exists("xtbopt.xyz"):
                        os.replace("xtbopt.xyz", "xtb.trj")
                
                if os.path.exists("xtb.trj"):
                    status.update(label="Calculation Finished!", state="complete")
                    st.rerun()
                else:
                    status.update(label="Calculation Failed", state="error")
                    st.error("Error Log:")
                    st.code(res.stderr)

    # SDF出力
    trj_exists = os.path.exists("xtb.trj")
    is_ready = True if (trj_exists and comp_name) else False
    
    # 簡易的なデータパース
    frames = []
    if trj_exists:
        try:
            with open("xtb.trj", 'r') as f:
                lines = f.readlines()
            if lines:
                num_atoms = int(lines[0].strip())
                for i in range(0, len(lines), num_atoms + 2):
                    frames.append(lines[i:i + num_atoms + 2])
        except: pass

    st.download_button("Download Result (SDF)", data="", disabled=not is_ready)

with col2:
    st.subheader("🔍 3D Viewer")
    if trj_exists and frames:
        idx = st.slider("Conformer", 1, len(frames), 1) if len(frames) > 1 else 1
        xyz_data = "".join(frames[idx-1])
        view = py3Dmol.view(width=450, height=450)
        view.addModel(xyz_data, 'xyz')
        view.setStyle({'stick': {'radius': 0.15}, 'sphere': {'scale': 0.25}})
        view.zoomTo()
        components.html(view._make_html(), height=460)
    else:
        st.info("Ready for input.")