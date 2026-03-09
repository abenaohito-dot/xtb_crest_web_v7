import streamlit as st
import pandas as pd
import numpy as np
import re
import os
import subprocess
import py3Dmol
import streamlit.components.v1 as components

# --- ページ設定 ---
st.set_page_config(page_title="XTB-CREST Web Console v8.0", layout="wide")

st.title("🧪 XTB-CREST Web Console v8.0")
st.write("Natural Products Conformer Ensemble Analyzer - GFN-FF Pre-opt Mode")

# --- 状態リセット関数 ---
def reset_all():
    files_to_remove = ["xtb.trj", "input.xyz", "crest_conformers.xyz", "crestopt.log", "xtbopt.xyz", "xtbopt.log"]
    for f in files_to_remove:
        if os.path.exists(f): os.remove(f)
    if "current_file" in st.session_state:
        del st.session_state.current_file
    st.rerun()

# --- サイドバー：設定 ---
with st.sidebar:
    st.header("⚙️ Settings")
    comp_name = st.text_input("Compound Name", value="", placeholder="Enter name to unlock downloads")
    
    calc_mode = st.radio("Calculation Method", ["CREST (Full Search)", "xTB (Optimization)"])
    solvent = st.selectbox("Solvent (ALPB)", ["methanol", "water", "chcl3", "benzene", "none"])
    cores = st.slider("CPU Cores", 1, 12, 4)
    
    st.divider()
    st.header("⚖️ Energy Thresholds (kcal/mol)")
    low_thresh = st.number_input("Threshold 1 (Low)", value=3.0, step=0.5)
    high_thresh = st.number_input("Threshold 2 (High)", value=10.0, step=1.0)
    
    if st.button("🗑️ Clear All & Reset"):
        reset_all()

# --- メインパネル ---
col1, col2 = st.columns([1, 1])

with col1:
    st.subheader("📥 Export Results")
    uploaded_file = st.file_uploader("Upload XYZ File", type=["xyz"])
    
    if uploaded_file:
        if "current_file" not in st.session_state or st.session_state.current_file != uploaded_file.name:
            if os.path.exists("xtb.trj"): os.remove("xtb.trj")
            st.session_state.current_file = uploaded_file.name

    trj_file = "xtb.trj"
    trj_exists = os.path.exists(trj_file)

    # --- 解析実行エリア ---
    if uploaded_file:
        st.write("---")
        st.markdown("### 🚀 Robust Calculation Control")
        # v8.0: 力場プリ最適化のオプション
        pre_opt = st.checkbox("Pre-optimize with GFN-FF", value=True, help="量子計算の前に分子力場(FF)で構造を整えます。エラー回避に有効です。")
        quick_mode = st.checkbox("Quick Mode (--quick)", value=True)
        
        if st.button("🚀 Run Analysis", type="primary"):
            with st.spinner("Processing... (Pre-optimization or Analysis in progress)"):
                with open("input.xyz", "wb") as f:
                    f.write(uploaded_file.getbuffer())
                
                target_xyz = "input.xyz"
                
                # Step 1: GFN-FFによる予備最適化 (v8.0 新機能)
                if pre_opt:
                    st.text("Step 1: Running GFN-FF Pre-optimization...")
                    ff_cmd = ["xtb", "input.xyz", "--gfnff", "--opt", "-T", str(cores)]
                    subprocess.run(ff_cmd, capture_output=True)
                    if os.path.exists("xtbopt.xyz"):
                        target_xyz = "xtbopt.xyz" # 成功したら力場構造を次に回す
                        st.text("GFN-FF Pre-optimization: Success")
                
                # Step 2: 本番解析
                st.text(f"Step 2: Running {calc_mode}...")
                if calc_mode == "CREST (Full Search)":
                    cmd = ["crest", target_xyz, "--gfn2", "-T", str(cores)]
                    if quick_mode: cmd.append("--quick")
                    if solvent != "none": cmd.extend(["--alpb", solvent])
                    subprocess.run(cmd)
                    if os.path.exists("crest_conformers.xyz"):
                        os.replace("crest_conformers.xyz", "xtb.trj")
                else:
                    cmd = ["xtb", target_xyz, "--opt", "--gfn2", "-T", str(cores)]
                    if solvent != "none": cmd.extend(["--alpb", solvent])
                    subprocess.run(cmd)
                    if os.path.exists("xtbopt.xyz"):
                        os.replace("xtbopt.xyz", "xtb.trj")
                
                st.rerun()

    # --- データの読み込み ---
    frames = []
    num_atoms, min_e = 0, 0.0
    if trj_exists:
        try:
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
        except:
            pass

    is_ready = True if (frames and comp_name) else False

    st.download_button(f"Download Low-E SDF", data="", disabled=not is_ready) # 簡略表示
    st.download_button(f"Download High-E SDF", data="", disabled=not is_ready)

with col2:
    st.subheader("🔍 3D Viewer")
    if is_ready:
        rank = st.slider("Select Conformer", 1, len(frames), 1) if len(frames) > 1 else 1
        f_view = frames[rank-1]
        xyz_data = "".join(f_view["coords"])
        view = py3Dmol.view(width=450, height=450)
        view.addModel(xyz_data, 'xyz')
        view.setStyle({'stick': {'radius': 0.15}, 'sphere': {'scale': 0.25}})
        view.zoomTo()
        components.html(view._make_html(), height=460)
        st.write(f"Relative Energy: **{(f_view['energy'] - min_e) * 627.509:.4f} kcal/mol**")
    else:
        st.info("Upload XYZ and run analysis. Try GFN-FF pre-opt if it fails.")