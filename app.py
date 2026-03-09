import streamlit as st
import os
import subprocess
import time
import numpy as np
import py3Dmol
import streamlit.components.v1 as components

st.set_page_config(page_title="XTB Web Console v9.7", layout="wide")

# --- 物理定数・変換 ---
AU_TO_KCAL = 627.509

def reset_all():
    for f in os.listdir("."):
        if f.endswith((".xyz", ".trj", ".pdb", ".mol", ".log", ".inp", ".out", ".sdf")):
            try: os.remove(f)
            except: pass
    st.session_state.clear()
    st.rerun()

# --- サイドバー ---
with st.sidebar:
    st.header("⚙️ Settings")
    comp_name = st.text_input("Compound Name", value="structure")
    calc_mode = st.radio("Method", ["xTB Opt", "xTB MD"])
    cores = st.slider("CPU Cores", 1, 12, 10)
    
    st.divider()
    st.header("🎯 Filtering & Export")
    # 先生のリクエスト：0-10 kcal/mol の選別
    e_window = st.slider("Energy Window (kcal/mol)", 0.0, 10.0, 3.0)
    rmsd_thr = st.slider("RMSD Clustering Threshold (Å)", 0.0, 2.0, 0.5)
    
    if st.button("🗑️ Full Reset"):
        reset_all()

col1, col2 = st.columns([1, 1])

with col1:
    st.subheader("📥 Input Control")
    uploaded_file = st.file_uploader("Upload Structure", type=["pdb", "mol", "xyz"])
    
    if uploaded_file and comp_name:
        in_ext = uploaded_file.name.split('.')[-1].lower()
        input_file = f"{comp_name}_input.{in_ext}"
        with open(input_file, "wb") as f:
            f.write(uploaded_file.getbuffer())
        
        if st.button("🚀 Start Calculation", type="primary"):
            with st.status("Calculating...", expanded=True) as status:
                # 実行前処理
                for f in ["xtbopt.xyz", "xtbopt.mol", "xtb.trj"]:
                    if os.path.exists(f): os.remove(f)

                if "Opt" in calc_mode:
                    cmd = ["xtb", input_file, "--opt", "--gfn", "2", "-P", str(cores)]
                    check_targets = ["xtbopt.mol", "xtbopt.xyz"]
                else:
                    with open("md.inp", "w") as f:
                        f.write("$md\n  temp=500\n  time=50.0\n  dump=100\n$end\n")
                    cmd = ["xtb", input_file, "--md", "--input", "md.inp", "--gfn", "2", "-P", str(cores)]
                    check_targets = ["xtb.trj"]
                
                result = subprocess.run(cmd, capture_output=True, text=True)
                
                # ファイル検出
                found = False
                for i in range(15):
                    for t in check_targets:
                        if os.path.exists(t):
                            time.sleep(1)
                            final_name = f"{comp_name}_raw.trj" if "MD" in calc_mode else f"{comp_name}_opt.{in_ext}"
                            os.replace(t, final_name)
                            st.session_state.current_raw = final_name
                            found = True; break
                    if found: break
                    time.sleep(1)
                
                if found:
                    status.update(label="Calculation Successful!", state="complete")
                    st.rerun()

    # --- 配座選別 & SDF変換ロジック ---
    if "current_raw" in st.session_state and st.session_state.current_raw.endswith(".trj"):
        st.divider()
        st.subheader("📦 Conformer Processing")
        if st.button("Generate Clean SDF (Filter & Cluster)"):
            with st.spinner("Processing trajectory..."):
                # 1. Trajectory読み込み & エネルギー抽出
                with open(st.session_state.current_raw, "r") as f:
                    lines = f.readlines()
                
                num_atoms = int(lines[0].strip())
                frame_size = num_atoms + 2
                frames = [lines[i:i+frame_size] for i in range(0, len(lines), frame_size)]
                
                data = []
                for frame in frames:
                    e_match = [l for l in frame if "energy:" in l]
                    if e_match:
                        e_val = float(e_match[0].split()[2])
                        coords = np.array([[float(x) for x in l.split()[1:4]] for l in frame[2:2+num_atoms]])
                        data.append({"energy": e_val, "coords": coords, "frame": frame})
                
                # 2. エネルギー選別
                min_e = min(d["energy"] for d in data)
                valid_data = [d for d in data if (d["energy"] - min_e) * AU_TO_KCAL <= e_window]
                valid_data.sort(key=lambda x: x["energy"]) # エネルギー順
                
                # 3. RMSDクラスタリング (簡易実装)
                final_conformers = []
                if valid_data:
                    final_conformers.append(valid_data[0])
                    for d in valid_data[1:]:
                        is_unique = True
                        for ref in final_conformers:
                            diff = d["coords"] - ref["coords"]
                            rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
                            if rmsd < rmsd_thr:
                                is_unique = False; break
                        if is_unique:
                            final_conformers.append(d)
                
                # 4. SDF書き出し (簡易変換)
                sdf_name = f"{comp_name}_conformers.sdf"
                with open(sdf_name, "w") as f:
                    for i, d in enumerate(final_conformers):
                        rel_e = (d["energy"] - min_e) * AU_TO_KCAL
                        f.write(f"{comp_name}_conf_{i}\n")
                        f.write(f" xTB-MD / RelE: {rel_e:.3f} kcal/mol\n\n")
                        f.write(f"{num_atoms:>3}{0:>3}  0  0  0  0  0  0  0  0999 V2000\n")
                        for j, l in enumerate(d["frame"][2:2+num_atoms]):
                            parts = l.split()
                            f.write(f"{float(parts[1]):>10.4f}{float(parts[2]):>10.4f}{float(parts[3]):>10.4f} {parts[0]:<3} 0  0  0  0  0  0  0  0  0  0  0  0\n")
                        f.write("M  END\n$$$$\n")
                
                st.session_state.current_sdf = sdf_name
                st.success(f"Extracted {len(final_conformers)} unique conformers within {e_window} kcal/mol.")

    # --- ダウンロードボタン ---
    if "current_sdf" in st.session_state:
        with open(st.session_state.current_sdf, "r") as f:
            st.download_button(f"💾 Download {st.session_state.current_sdf}", data=f.read(), file_name=st.session_state.current_sdf)

with col2:
    st.subheader("🔍 3D Viewer")
    # (Viewer部分は前述のTrajectory再生モードを維持)