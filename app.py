import streamlit as st
import pandas as pd
import numpy as np
import re
import os
import subprocess  # コマンド実行用
import py3Dmol
import streamlit.components.v1 as components

# --- ページ設定 ---
st.set_page_config(page_title="XTB-CREST Web Console v7.5", layout="wide")

st.title("🧪 XTB-CREST Web Console v7.5")
st.write("Natural Products Conformer Ensemble Analyzer - Real Engine Mode")

# --- サイドバー：設定 ---
with st.sidebar:
    st.header("⚙️ Settings")
    comp_name = st.text_input("Compound Name", value="", placeholder="Enter name to unlock downloads")
    
    # 計算オプション
    mode = st.radio("Execution Mode", ["iMTD-GC (Full Search)", "Optimization Only"])
    solvent = st.selectbox("Solvent (ALPB)", ["methanol", "water", "chcl3", "benzene", "none"])
    cores = st.slider("CPU Cores", 1, 12, 4) # M4 Proに合わせて調整
    
    st.divider()
    st.header("⚖️ Energy Thresholds (kcal/mol)")
    low_thresh = st.number_input("Threshold 1 (Low)", value=3.0, step=0.5)
    high_thresh = st.number_input("Threshold 2 (High)", value=10.0, step=1.0)

# --- メインパネル ---
col1, col2 = st.columns([1, 1])

with col1:
    st.subheader("📥 Export Results")
    uploaded_file = st.file_uploader("Upload XYZ File", type=["xyz"])
    
    if uploaded_file:
        with open("input.xyz", "wb") as f:
            f.write(uploaded_file.getbuffer())
        
        # 解析実行ボタン（本物）
        if st.button("🚀 Run CREST Analysis (Real calculation)"):
            with st.spinner("CREST is running... Please wait (Check terminal for progress)"):
                # コマンド構築
                cmd = ["crest", "input.xyz", "--gfn2", "-T", str(cores)]
                if solvent != "none":
                    cmd.extend(["--alpb", solvent])
                
                # 実行
                result = subprocess.run(cmd, capture_output=True, text=True)
                
                if os.path.exists("crest_conformers.xyz"):
                    # 結果ファイルをアプリが認識する名前に変更
                    os.replace("crest_conformers.xyz", "xtb.trj")
                    st.success("Analysis Complete!")
                    st.rerun()
                else:
                    st.error("CREST failed to produce results. Check terminal.")
                    st.code(result.stderr)

    # 以下、解析・表示ロジックは継続
    trj_file = "xtb.trj"
    trj_exists = os.path.exists(trj_file)
    is_ready = True if (trj_exists and comp_name) else False

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

    # SDF生成とボタン（グレーアウト制御）
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

    s1_data = get_sdf(low_thresh)
    st.download_button(label=f"Download {comp_name if comp_name else '---'}_{low_thresh}kcal.sdf ({len(frames)} total)", 
                       data=s1_data, file_name=f"{comp_name}_{low_thresh}kcal.sdf", disabled=not is_ready)
    
    s2_data = get_sdf(high_thresh)
    st.download_button(label=f"Download {comp_name if comp_name else '---'}_{high_thresh}kcal.sdf", 
                       data=s2_data, file_name=f"{comp_name}_{high_thresh}kcal.sdf", disabled=not is_ready)

with col2:
    st.subheader("🔍 3D Viewer")
    if is_ready and frames:
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
        st.info("Upload and Run analysis to see conformers.")