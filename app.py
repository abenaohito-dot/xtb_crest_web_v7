import streamlit as st
import pandas as pd
import numpy as np
import re
import os
import py3Dmol
import streamlit.components.v1 as components

# --- ページ設定 ---
st.set_page_config(page_title="XTB-CREST Web Console", layout="wide")

# カスタムCSSでデザインを調整
st.markdown("""
    <style>
    .main { background-color: #0f172a; }
    .stButton>button { width: 100%; background-color: #38bdf8; color: white; border-radius: 5px; }
    </style>
    """, unsafe_allow_html=True)

st.title("XTB-CREST Web Console v7.0")
st.write("Natural Products Conformer Ensemble Analyzer")

# --- サイドバー：設定 ---
with st.sidebar:
    st.header("⚙️ Settings")
    # デフォルト名を空に変更
    comp_name = st.text_input("Compound Name", value="", placeholder="e.g. Stilbene_Oligomer")
    
    mode = st.radio("Execution Mode", ["xTB-Opt", "xTB-MD", "CREST-MTD"])
    solvent = st.selectbox("Solvent (ALPB)", ["methanol", "water", "chcl3", "benzene", "none"])
    cores = st.slider("CPU Cores", 1, 16, 4)
    
    st.divider()
    st.header("⚖️ Energy Thresholds (kcal/mol)")
    # 可変しきい値の設定
    low_thresh = st.number_input("Threshold 1 (Low)", value=3.0, step=0.5)
    high_thresh = st.number_input("Threshold 2 (High)", value=10.0, step=1.0)

# --- メインパネル ---
uploaded_file = st.file_uploader("Upload XYZ File", type=["xyz"])

if uploaded_file:
    # ファイル名からデフォルト名を推測（入力がない場合）
    display_name = comp_name if comp_name else uploaded_file.name.rsplit('.', 1)[0]
    
    # ファイルの一時保存
    with open("input.xyz", "wb") as f:
        f.write(uploaded_file.getbuffer())
    
    col1, col2 = st.columns([1, 1])

    # 計算結果ファイル（xtb.trj）の解析ロジック
    trj_file = "xtb.trj" 
    if os.path.exists(trj_file):
        with open(trj_file, 'r') as f:
            lines = f.readlines()
        
        # 構造のパース
        num_atoms = int(lines[0].strip())
        frames = []
        chunk_size = num_atoms + 2
        for i in range(0, len(lines), chunk_size):
            chunk = lines[i:i+chunk_size]
            if len(chunk) < chunk_size: break
            # エネルギーの抽出
            e_match = re.search(r"energy:\s+([-+]?\d+\.\d+)", chunk[1])
            frames.append({
                "energy": float(e_match.group(1)) if e_match else 0.0, 
                "coords": chunk
            })
        
        # エネルギー順にソート
        frames.sort(key=lambda x: x["energy"])
        min_e = frames[0]["energy"]

        with col1:
            st.subheader("📥 Export Results")
            
            # SDF生成関数
            def make_sdf(threshold):
                sdf = ""
                count = 0
                for i, f in enumerate(frames, start=1):
                    rel_e = (f["energy"] - min_e) * 627.509
                    if rel_e > threshold: break
                    sdf += f"{display_name}_{i}\nThreshold: {threshold} kcal\n\n"
                    sdf += f"{num_atoms:>3}  0  0  0  0  0  0  0  0  0999 V2000\n"
                    for line in f["coords"][2:]:
                        p = line.split()
                        if len(p) >= 4:
                            sdf += f"{float(p[1]):>10.4f}{float(p[2]):>10.4f}{float(p[3]):>10.4f} {p[0]:<3} 0  0  0  0  0\n"
                    sdf += "M  END\n> <ENERGY_KCAL>\n{:.4f}\n\n$$$$\n".format(rel_e)
                    count += 1
                return sdf, count

            # ダウンロードボタン：Threshold 1
            s1, c1 = make_sdf(low_thresh)
            st.download_button(
                label=f"Download {display_name}_{low_thresh}kcal.sdf ({c1} structs)",
                data=s1,
                file_name=f"{display_name}_{low_thresh}kcal.sdf",
                mime="chemical/x-mdlsdf"
            )
            
            # ダウンロードボタン：Threshold 2
            s2, c2 = make_sdf(high_thresh)
            st.download_button(
                label=f"Download {display_name}_{high_thresh}kcal.sdf ({c2} structs)",
                data=s2,
                file_name=f"{display_name}_{high_thresh}kcal.sdf",
                mime="chemical/x-mdlsdf"
            )
            
            # エネルギー統計の表示
            st.info(f"Total conformers found: {len(frames)}")

        with col2:
            st.subheader("🔍 3D Viewer")
            rank = st.slider("Select Conformer (by Energy Rank)", 1, len(frames), 1)
            f_view = frames[rank-1]
            e_rel_now = (f_view["energy"] - min_e) * 627.509
            
            # 3Dビューワーの表示 (py3Dmol HTML埋め込み)
            xyz_data = "".join(f_view["coords"])
            view = py3Dmol.view(width=450, height=450)
            view.addModel(xyz_data, 'xyz')
            view.setStyle({'stick': {'radius': 0.15}, 'sphere': {'scale': 0.25}})
            view.zoomTo()
            
            # コンポーネントとして出力
            components.html(view._make_html(), height=460)
            st.write(f"Relative Energy: **{e_rel_now:.4f} kcal/mol**")

    else:
        st.warning("No calculation results (xtb.trj) detected. Please ensure xtb has been run in the working directory.")
        if st.button("🚀 Launch Test (Create dummy results)"):
            # テスト用に現在の入力をダミーのtrjとして保存する機能（デバッグ用）
            with open("xtb.trj", "w") as f:
                f.write(f"{uploaded_file.getvalue().decode()}\n")
            st.rerun()