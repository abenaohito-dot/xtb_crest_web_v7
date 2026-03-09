import streamlit as st
import os
import subprocess
import time
import re
import py3Dmol
import streamlit.components.v1 as components

# --- ページ設定 ---
st.set_page_config(page_title="XTB-CREST Web Console v9.4", layout="wide")

# M4 Pro 安定化設定
env = os.environ.copy()
env["OMP_STACKSIZE"] = "2G"
env["OMP_NUM_THREADS"] = "4"

st.title("🧪 XTB-CREST Web Console v9.4")
st.write("Direct Format Processing & File Sync Mode")

# --- 状態クリア ---
def reset_all():
    # 関連ファイルを一掃
    files = ["xtb.trj", "input.pdb", "input.mol", "input.xyz", "crest_conformers.xyz", "xtbopt.xyz", "xtbtopo.mol", "xtbopt.log"]
    for f in files:
        if os.path.exists(f):
            try: os.remove(f)
            except: pass
    if "current_file" in st.session_state:
        del st.session_state.current_file
    st.rerun()

# --- サイドバー ---
with st.sidebar:
    st.header("⚙️ Settings")
    comp_name = st.text_input("Compound Name", value="", placeholder="Required for download")
    calc_mode = st.radio("Method", ["xTB (Optimization Only)", "CREST (Conformer Search)"])
    cores = st.slider("CPU Cores", 1, 12, 4)
    
    st.divider()
    if st.button("🗑️ Full Reset"):
        reset_all()

col1, col2 = st.columns([1, 1])

with col1:
    st.subheader("📥 Input Control")
    # PDBとMOLを直接受け付ける
    uploaded_file = st.file_uploader("Upload (PDB or MOL)", type=["pdb", "mol"])
    
    if uploaded_file:
        ext = uploaded_file.name.split('.')[-1].lower()
        target_file = f"input.{ext}"
        
        # ファイルをそのまま保存（中身を一切いじらない）
        with open(target_file, "wb") as f:
            f.write(uploaded_file.getbuffer())
        
        if st.button("🚀 Start Calculation", type="primary"):
            with st.status("Calculating...", expanded=True) as status:
                # 実行前に古い出力を削除
                for f in ["xtbopt.xyz", "crest_conformers.xyz", "xtb.trj"]:
                    if os.path.exists(f): os.remove(f)

                # コマンド構築
                if "xTB" in calc_mode:
                    # --opt で最適化、--gfn 2 で計算
                    cmd = ["xtb", target_file, "--opt", "--gfn", "2", "-T", str(cores)]
                else:
                    # CREST実行
                    cmd = ["crest", target_file, "--gfn2", "-T", str(cores), "--quick"]
                
                st.write(f"Running: {' '.join(cmd)}")
                result = subprocess.run(cmd, env=env, capture_output=True, text=True)
                
                # ファイル生成の同期（M4 Proの書き込み待ち）
                found = False
                check_target = "xtbopt.xyz" if "xTB" in calc_mode else "crest_conformers.xyz"
                
                for i in range(10): # 最大10秒待機
                    if os.path.exists(check_target):
                        # 読み込みエラーを防ぐため、コピーが完了するまで微小待機
                        time.sleep(0.5)
                        os.replace(check_target, "xtb.trj")
                        found = True
                        break
                    time.sleep(1)
                    st.write(f"Waiting for output... ({i+1}s)")

                if found:
                    status.update(label="Calculation Successful!", state="complete")
                    st.rerun()
                else:
                    status.update(label="File Not Found", state="error")
                    st.error("The output file was not generated. Check the log below:")
                    st.code(result.stderr)

    # --- ダウンロードボタン ---
    if os.path.exists("xtb.trj") and comp_name:
        with open("xtb.trj", "r") as f:
            st.download_button(f"Download {comp_name}.sdf", data=f.read(), file_name=f"{comp_name}.sdf")

with col2:
    st.subheader("🔍 3D Viewer")
    if os.path.exists("xtb.trj"):
        with open("xtb.trj", "r") as f:
            lines = f.readlines()
            if lines:
                num = int(lines[0].strip())
                xyz_data = "".join(lines[:num+2])
                view = py3Dmol.view(width=450, height=450)
                view.addModel(xyz_data, 'xyz')
                view.setStyle({'stick': {}, 'sphere': {'scale': 0.3}})
                view.zoomTo()
                components.html(view._make_html(), height=460)
                
    elif uploaded_file:
        st.info("Calculation is required to view optimized structure.")