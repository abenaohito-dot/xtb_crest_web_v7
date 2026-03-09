import streamlit as st
import pandas as pd
import numpy as np
import re
import os
import subprocess
import py3Dmol
import streamlit.components.v1 as components

# --- ページ設定 ---
st.set_page_config(page_title="XTB-CREST Web Console v9.1", layout="wide")

# M4 Pro クラッシュ対策（スタックメモリ拡張） [cite: 3]
env = os.environ.copy()
env["OMP_STACKSIZE"] = "2G"
env["OMP_NUM_THREADS"] = "4"

st.title("🧪 XTB-CREST Web Console v9.1")
st.write("Full Integrated Suite - Dual Search Engine Mode")

# --- 座標抽出関数 (PDB/MOL/SDF) ---
def parse_to_xyz(file_obj, ext):
    lines = file_obj.getvalue().decode().splitlines()
    atoms = []
    if ext == "pdb":
        for l in lines:
            if l.startswith("ATOM") or l.startswith("HETATM"):
                elem = l[76:78].strip() or l[12:14].strip().rstrip('0123456789')
                atoms.append(f"{elem} {l[30:38]} {l[38:46]} {l[46:54]}")
    elif ext == "mol" or ext == "sdf":
        for l in lines[4:]:
            p = l.split()
            if len(p) >= 10 and p[3].isalpha():
                atoms.append(f"{p[3]} {p[0]} {p[1]} {p[2]}")
            if l.startswith("M  END"): break
    if atoms:
        return f"{len(atoms)}\nconverted from {ext}\n" + "\n".join(atoms) + "\n"
    return file_obj.getvalue().decode()

# --- 状態リセット ---
def reset_all():
    files = ["xtb.trj", "input.xyz", "crest_conformers.xyz", "xtbopt.xyz", "crestopt.log", "xtbopt.log", "xtb_confs.xyz"]
    for f in files:
        if os.path.exists(f): os.remove(f)
    if "current_file" in st.session_state:
        del st.session_state.current_file
    st.rerun()

# --- サイドバー：設定 ---
with st.sidebar:
    st.header("⚙️ Settings")
    comp_name = st.text_input("Compound Name", value="", placeholder="e.g. Stilbene_Oligomer")
    
    # v9.1: 手法の選択肢を3つに拡充
    calc_mode = st.radio("Method", [
        "CREST (Conformer Search)", 
        "xTB (Simple Search)", 
        "xTB (Optimization Only)"
    ])
    
    solvent = st.selectbox("Solvent", ["methanol", "water", "chcl3", "benzene", "none"])
    cores = st.slider("Cores", 1, 12, 4)
    
    st.divider()
    st.header("⚖️ Thresholds (kcal/mol)")
    low_thresh = st.number_input("Threshold 1 (Low)", value=3.0, step=0.5)
    high_thresh = st.number_input("Threshold 2 (High)", value=10.0, step=1.0)
    
    if st.button("🗑️ Full Reset"):
        reset_all()

# --- メインパネル ---
col1, col2 = st.columns([1, 1])

with col1:
    st.subheader("📥 Input & Control")
    uploaded_file = st.file_uploader("Upload Structure (PDB/MOL/XYZ)", type=["xyz", "mol", "pdb", "sdf"])
    
    if uploaded_file:
        ext = uploaded_file.name.split('.')[-1].lower()
        if "current_file" not in st.session_state or st.session_state.current_file != uploaded_file.name:
            if os.path.exists("xtb.trj"): os.remove("xtb.trj")
            st.session_state.current_file = uploaded_file.name
        
        with open("input.xyz", "w") as f:
            f.write(parse_to_xyz(uploaded_file, ext))

        if st.button("🚀 Start Calculation", type="primary"):
            with st.status("Executing Analysis...", expanded=True) as status:
                env["OMP_NUM_THREADS"] = str(cores)
                
                # Step 1: 共通のプリ最適化
                st.write("Step 1: GFN-FF Pre-opt...")
                subprocess.run(["xtb", "input.xyz", "--gfnff", "--opt", "-T", str(cores)], env=env)
                target = "xtbopt.xyz" if os.path.exists("xtbopt.xyz") else "input.xyz"
                
                # Step 2: 選択手法による計算
                st.write(f"Step 2: Running {calc_mode}...")
                
                if "CREST" in calc_mode:
                    cmd = ["crest", target, "--gfn2", "-T", str(cores), "--quick"]
                    if solvent != "none": cmd.extend(["--alpb", solvent])
                    subprocess.run(cmd, env=env)
                    if os.path.exists("crest_conformers.xyz"):
                        os.replace("crest_conformers.xyz", "xtb.trj")
                
                elif "Simple Search" in calc_mode:
                    # xTB自体の簡易探索機能
                    cmd = ["xtb", target, "--conf", "--gfn2", "-T", str(cores)]
                    if solvent != "none": cmd.extend(["--alpb", solvent])
                    subprocess.run(cmd, env=env)
                    if os.path.exists("xtb_confs.xyz"):
                        os.replace("xtb_confs.xyz", "xtb.trj")
                
                else:
                    # 通常の1構造最適化
                    cmd = ["xtb", target, "--opt", "--gfn2", "-T", str(cores)]
                    if solvent != "none": cmd.extend(["--alpb", solvent])
                    subprocess.run(cmd, env=env)
                    if os.path.exists("xtbopt.xyz"):
                        os.replace("xtbopt.xyz", "xtb.trj")
                
                if os.path.exists("xtb.trj"):
                    status.update(label="Success!", state="complete")
                    st.rerun()

    # --- 解析とSDF生成ロジック ---
    frames = []
    num_atoms, min_e = 0, 0.0
    if os.path.exists("xtb.trj"):
        with open("xtb.trj", 'r') as f:
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
        if not frames or not comp_name: return ""
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

    st.download_button(f"Download {comp_name}_{low_thresh}kcal.sdf", data=make_sdf(low_thresh), file_name=f"{comp_name}_{low_thresh}kcal.sdf", disabled=not (frames and comp_name))
    st.download_button(f"Download {comp_name}_{high_thresh}kcal.sdf", data=make_sdf(high_thresh), file_name=f"{comp_name}_{high_thresh}kcal.sdf", disabled=not (frames and comp_name))

with col2:
    st.subheader("🔍 3D Viewer")
    if frames:
        if len(frames) > 1:
            rank = st.slider("Select Structure", 1, len(frames), 1)
        else:
            st.info("Single structure optimized.")
            rank = 1
            
        f_view = frames[rank-1]
        xyz_data = "".join(f_view["coords"])
        view = py3Dmol.view(width=450, height=450)
        view.addModel(xyz_data, 'xyz')
        view.setStyle({'stick': {}, 'sphere': {'scale': 0.3}})
        view.zoomTo()
        components.html(view._make_html(), height=460)
        st.write(f"Relative Energy: **{(f_view['energy'] - min_e) * 627.509:.4f} kcal/mol**")
    else:
        st.info("Waiting for calculation results...")