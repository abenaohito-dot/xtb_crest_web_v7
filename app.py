import streamlit as st
import os
import subprocess
import time
import numpy as np
import py3Dmol
import streamlit.components.v1 as components

# --- Page Configuration ---
st.set_page_config(page_title="XTB Web Console", layout="wide")

# M4 Pro Environment Settings
env = os.environ.copy()
env["OMP_STACKSIZE"] = "2G"

# --- Constants ---
AU_TO_KCAL = 627.509

def reset_all():
    exts = (".xyz", ".trj", ".pdb", ".mol", ".log", ".inp", ".out", ".sdf", "xtbrestart", "wbo")
    for f in os.listdir("."):
        if f.endswith(exts) or "xtb" in f:
            try: os.remove(f)
            except: pass
    st.session_state.clear()
    st.rerun()

# --- Title ---
st.title("🧪 XTB Web Console v9.8")

# --- Sidebar (Settings) ---
with st.sidebar:
    st.header("Settings")
    comp_name = st.text_input("Compound Name", value="structure")
    calc_mode = st.radio("Method", ["xTB Opt", "xTB MD"])
    cores = st.slider("CPU Cores", 1, 12, 11)
    
    st.divider()
    st.header("Filter & Export")
    e_window = st.slider("Energy Window (kcal/mol)", 0.0, 10.0, 3.0)
    rmsd_thr = st.slider("RMSD Threshold (Å)", 0.0, 2.0, 0.5)
    
    if st.button("Full Reset"):
        reset_all()

col1, col2 = st.columns([1, 1])

with col1:
    st.subheader("Input & Control")
    uploaded_file = st.file_uploader("Upload Structure (PDB/MOL/XYZ)", type=["pdb", "mol", "xyz"])
    
    if uploaded_file and comp_name:
        in_ext = uploaded_file.name.split('.')[-1].lower()
        input_file = f"{comp_name}_input.{in_ext}"
        with open(input_file, "wb") as f:
            f.write(uploaded_file.getbuffer())
        
        if st.button("Start Calculation", type="primary"):
            with st.status("Calculating...", expanded=True) as status:
                # Clean up old files
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
                
                result = subprocess.run(cmd, env=env, capture_output=True, text=True)
                
                # File Detection
                found = False
                for i in range(15):
                    for t in check_targets:
                        if os.path.exists(t):
                            time.sleep(1)
                            final_name = f"{comp_name}_raw.trj" if "MD" in calc_mode else f"{comp_name}_opt.{in_ext}"
                            os.replace(t, final_name)
                            st.session_state.current_raw = final_name
                            st.session_state.output_ext = t.split('.')[-1]
                            found = True; break
                    if found: break
                    time.sleep(1)
                
                if found:
                    status.update(label="Success!", state="complete")
                    st.rerun()
                else:
                    st.error("Output Not Found")
                    st.code(result.stdout[-1000:])

    # --- Processing & SDF Export ---
    if "current_raw" in st.session_state and st.session_state.current_raw.endswith(".trj"):
        st.divider()
        st.subheader("Conformer Processing")
        if st.button("Generate SDF (Filter & Cluster)"):
            with st.spinner("Processing..."):
                with open(st.session_state.current_raw, "r") as f:
                    lines = f.readlines()
                
                num_atoms = int(lines[0].strip())
                f_size = num_atoms + 2
                frames = [lines[i:i+f_size] for i in range(0, len(lines), f_size)]
                
                data = []
                for f in frames:
                    e_m = [l for l in f if "energy:" in l]
                    if e_m:
                        e_val = float(e_m[0].split()[2])
                        c = np.array([[float(x) for x in l.split()[1:4]] for l in f[2:2+num_atoms]])
                        data.append({"energy": e_val, "coords": c, "frame": f})
                
                min_e = min(d["energy"] for d in data)
                valid = [d for d in data if (d["energy"] - min_e) * AU_TO_KCAL <= e_window]
                valid.sort(key=lambda x: x["energy"])
                
                final = []
                if valid:
                    final.append(valid[0])
                    for d in valid[1:]:
                        unique = True
                        for ref in final:
                            diff = d["coords"] - ref["coords"]
                            if np.sqrt(np.mean(np.sum(diff**2, axis=1))) < rmsd_thr:
                                unique = False; break
                        if unique: final.append(d)
                
                sdf_name = f"{comp_name}_conformers.sdf"
                with open(sdf_name, "w") as f:
                    for i, d in enumerate(final):
                        rel_e = (d["energy"] - min_e) * AU_TO_KCAL
                        f.write(f"{comp_name}_conf_{i}\n xTB-MD / RelE: {rel_e:.3f} kcal/mol\n\n")
                        f.write(f"{num_atoms:>3}  0  0  0  0  0  0  0  0  0999 V2000\n")
                        for l in d["frame"][2:2+num_atoms]:
                            p = l.split()
                            f.write(f"{float(p[1]):>10.4f}{float(p[2]):>10.4f}{float(p[3]):>10.4f} {p[0]:<3} 0  0  0  0  0  0  0  0  0  0  0  0\n")
                        f.write("M  END\n$$$$\n")
                
                st.session_state.current_sdf = sdf_name
                st.success(f"Extracted {len(final)} conformers.")

    if "current_sdf" in st.session_state:
        st.download_button(f"Download SDF", data=open(st.session_state.current_sdf, "rb"), file_name=st.session_state.current_sdf)

with col2:
    st.subheader("3D Viewer")
    if "current_raw" in st.session_state:
        with open(st.session_state.current_raw, "r") as f:
            data = f.read()
            if data:
                fmt = st.session_state.get("output_ext", "xyz")
                view = py3Dmol.view(width=500, height=500)
                if ".trj" in st.session_state.current_raw:
                    view.addModelsAsFrames(data, 'xyz')
                    view.animate({'loop': 'backward', 'step': 5})
                else:
                    view.addModel(data, fmt)
                view.setStyle({'stick': {'radius': 0.15}, 'sphere': {'scale': 0.25}})
                view.zoomTo()
                components.html(view._make_html(), height=510)