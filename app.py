import streamlit as st
import os
import subprocess
import re
import py3Dmol
import streamlit.components.v1 as components

st.set_page_config(page_title="XTB-CREST Web Console v8.7", layout="wide")
st.title("🧪 XTB-CREST Web Console v8.7")

# 阿部先生の Mac mini (M4 Pro) で Segmentation fault を防ぐための環境変数
env = os.environ.copy()
env["OMP_STACKSIZE"] = "1G"

# --- サイドバー ---
with st.sidebar:
    st.header("⚙️ Settings")
    comp_name = st.text_input("Name", value="", placeholder="Required")
    calc_mode = st.radio("Method", ["CREST (Search)", "xTB (Opt)"])
    cores = st.slider("Cores", 1, 12, 4)
    if st.button("🗑️ Reset All"):
        for f in ["xtb.trj", "input.xyz", "crest_conformers.xyz", "xtbopt.xyz"]:
            if os.path.exists(f): os.remove(f)
        st.rerun()

col1, col2 = st.columns([1, 1])

with col1:
    st.subheader("📥 Input")
    uploaded_file = st.file_uploader("Upload XYZ", type=["xyz"])
    
    if uploaded_file:
        # 1. 複雑な変換をせず、シンプルに保存
        with open("input.xyz", "wb") as f:
            f.write(uploaded_file.getbuffer())

        # 2. 実行ボタン
        if st.button("🚀 Run Calculation", type="primary"):
            # 以前の残骸を確実に消去
            for f in ["xtb.trj", "crest_conformers.xyz", "xtbopt.xyz"]:
                if os.path.exists(f): os.remove(f)
            
            with st.status("Calculating...", expanded=True):
                # OMP環境変数を適用して実行
                env["OMP_NUM_THREADS"] = str(cores)
                
                if "CREST" in calc_mode:
                    # 阿部先生の環境で最も安定していた --quick オプション
                    cmd = ["crest", "input.xyz", "--gfn2", "-T", str(cores), "--quick"]
                    subprocess.run(cmd, env=env)
                    if os.path.exists("crest_conformers.xyz"):
                        os.replace("crest_conformers.xyz", "xtb.trj")
                else:
                    cmd = ["xtb", "input.xyz", "--opt", "--gfn2", "-T", str(cores)]
                    subprocess.run(cmd, env=env)
                    if os.path.exists("xtbopt.xyz"):
                        os.replace("xtbopt.xyz", "xtb.trj")
            
            if os.path.exists("xtb.trj"):
                st.success("Analysis Complete!")
                st.rerun()
            else:
                st.error("Calculation failed. Segmentation fault may have occurred.")

    # 3. ダウンロードボタン（以前の動作を復元）
    if os.path.exists("xtb.trj") and comp_name:
        with open("xtb.trj", "r") as f:
            st.download_button(f"Download {comp_name}.sdf", data=f.read(), file_name=f"{comp_name}.sdf")

with col2:
    st.subheader("🔍 3D Viewer")
    if os.path.exists("xtb.trj"):
        with open("xtb.trj", "r") as f:
            lines = f.readlines()
            num_atoms = int(lines[0].strip())
            xyz_data = "".join(lines[:num_atoms+2])
            view = py3Dmol.view(width=450, height=450)
            view.addModel(xyz_data, 'xyz')
            view.setStyle({'stick': {}, 'sphere': {'scale': 0.3}})
            view.zoomTo()
            components.html(view._make_html(), height=460)