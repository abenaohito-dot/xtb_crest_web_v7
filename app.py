import streamlit as st
import pandas as pd
import numpy as np
import re
import os
import subprocess
import py3Dmol
import streamlit.components.v1 as components

# --- ページ設定 ---
st.set_page_config(page_title="XTB-CREST Web Console v8.3", layout="wide")

st.title("🧪 XTB-CREST Web Console v8.3")
st.write("Natural Products Conformer Ensemble Analyzer - Professional Edition")

# --- 状態リセット関数 ---
def reset_all():
    files_to_remove = ["xtb.trj", "input.xyz", "crest_conformers.xyz", "crestopt.log", "xtbopt.xyz", "xtbopt.log", "raw_input.*"]
    for f in files_to_remove:
        if os.path.exists(f): 
            try: os.remove(f)
            except: pass
    if "current_file" in st.session_state:
        del st.session_state.current_file
    st.rerun()

# --- サイドバー：設定（v8.0の全機能を維持） ---
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
    
    st.divider()
    if st.button("🗑️ Clear All & Reset"):
        reset_all()

# --- メインパネル ---
col1, col2 = st.columns([1, 1])

with col1:
    st.subheader("📥 Export Results")
    # v8.3: 対応形式を拡大 (MOL, PDBに対応)
    uploaded_file = st.file_uploader("Upload Structure File", type=["xyz", "mol", "pdb", "sdf"])
    
    if uploaded_file:
        file_ext = uploaded_file.name.split('.')[-1].lower()
        if "current_file" not in st.session_state or st.session_state.current_file != uploaded_file.name:
            if os.path.exists("xtb.trj"): os.remove("xtb.trj")
            st.session_state.current_file = uploaded_file.name

        # 一時保存とフォーマット変換
        temp_raw = f"raw_input.{file_ext}"
        with open(temp_raw, "wb") as f:
            f.write(uploaded_file.getbuffer())
        
        # XYZへの変換・洗浄
        if file_ext != "xyz":
            subprocess.run(["xtb", temp_raw, "--steps", "0"], capture_output=True)
            if os.path.exists("xtbopt.xyz"):
                os.replace("xtbopt.xyz", "input.xyz")
        else:
            os.replace(temp_raw, "input.xyz")

    trj_file = "xtb.trj"
    trj_exists = os.path.exists(trj_file)
    is_ready = True if (trj_exists and comp_name) else False

    # --- 解析実行エリア（v8.0のUIを復元） ---
    if uploaded_file and os.path.exists("input.xyz"):
        st.write("---")
        st.markdown("### 🚀 Calculation Control")
        pre_opt = st.checkbox("Pre-optimize with GFN-FF", value=True)
        quick_mode = st.checkbox("Quick Mode (--quick)", value=True)
        
        b_col1, b_col2 = st.columns(2)
        with b_col1:
            btn_label = "🔄 Re-run Analysis" if trj_exists else "🚀 Run Analysis"
            if st.button(btn_label, type="primary"):
                with st.spinner("Calculating..."):
                    target = "input.xyz"
                    if pre_opt:
                        subprocess.run(["xtb", "input.xyz", "--gfnff", "--opt", "-T", str(cores)])
                        if os.path.exists("xtbopt.xyz"): target = "xtbopt.xyz"
                    
                    if calc_mode == "CREST (Full Search)":
                        cmd = ["crest", target, "--gfn2", "-T", str(cores)]
                        if quick_mode: cmd.append("--quick")
                        if solvent != "none": cmd.extend(["--alpb", solvent])
                        subprocess.run(cmd)
                        if os.path.exists("crest_conformers.xyz"):
                            os.replace("crest_conformers.xyz", "xtb.trj")
                    else:
                        cmd = ["xtb", target, "--opt", "--gfn2", "-T", str(cores)]
                        if solvent != "none": cmd.extend(["--alpb", solvent])
                        subprocess.run(cmd)
                        if os.path.exists("xtbopt.xyz"):
                            os.replace("xtbopt.xyz", "xtb.trj")
                    st.rerun()

        with b_col2:
            if st.button("🛠️ Launch Test (Dummy)"):
                with open("xtb.trj", "wb") as f:
                    f.write(open("input.xyz", "rb").read())
                st.rerun()

    # --- データ解析とSDF生成（v8.0のロジックを復元） ---
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
        except: pass

    def get_sdf(threshold):
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

    st.download_button(f"Download {comp_name if comp_name else '---'}_{low_thresh}kcal.sdf", 
                       data=get_sdf(low_thresh), file_name=f"{comp_name}_{low_thresh}kcal.sdf", disabled=not is_ready)
    st.download_button(f"Download {comp_name if comp_name else '---'}_{high_thresh}kcal.sdf", 
                       data=get_sdf(high_thresh), file_name=f"{comp_name}_{high_thresh}kcal.sdf", disabled=not is_ready)

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
        st.info("Upload XYZ/MOL/PDB and run analysis.")