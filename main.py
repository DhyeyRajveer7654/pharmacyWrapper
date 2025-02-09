import streamlit as st
import streamlit.components.v1 as components
import prompts
import chat_with_gpt
from string import Template
from rdkit import Chem
from rdkit.Chem import Draw
import requests

size = (250, 250)

# Set Page Configuration
st.set_page_config(page_title="QAI Model", layout="wide", page_icon="🧪")

# Apply Custom Styles
st.markdown("""
    <style>
        body { background-color: #0e1117; color: white; font-family: 'Arial', sans-serif; }
        .stTextInput>div>div>input, .stSelectbox>div>div>select, .stTextArea>div>textarea { 
            border-radius: 5px !important; padding: 10px;
        }
        div[data-testid="stTextInput"] input,
        div[data-testid="stTextArea"] textarea,
        div[data-testid="stNumberInput"] input,
        div[data-testid="stSelectbox"] > div[data-baseweb="select"],
        div[data-testid="stSlider"] > div[data-baseweb="slider"] {
            background-color: white !important;
            color: black !important;
            border: 1px solid #ccc !important;
            border-radius: 10px !important;
            box-shadow: none !important;
        }
        div[data-testid="stTextInput"] input:hover,
        div[data-testid="stTextArea"] textarea:hover,
        div[data-testid="stNumberInput"] input:hover,
        div[data-testid="stSelectbox"] > div[role="combobox"]:hover {
            border: 1px solid #888 !important;
        }
        div[data-testid="stTextInput"] input:focus,
        div[data-testid="stTextArea"] textarea:focus,
        div[data-testid="stNumberInput"] input:focus,
        div[data-testid="stSelectbox"] > div[role="combobox"]:focus {
            border: 1px solid #007BFF !important;
            outline: none !important;
        }
        .stButton>button:hover { 
            background: linear-gradient(90deg, #00D4FF, #007BFF);
            transform: scale(1.05);
        }
        .title { color: #00D4FF; text-align: center; font-size: 30px; font-weight: bold; }
        .subtitle { color: #cccccc; text-align: center; font-size: 18px; }
        .card {
            background-color: #1e222a; padding: 20px; border-radius: 12px; 
            box-shadow: 0px 4px 10px rgba(255, 255, 255, 0.1); margin: 20px;
        }
    </style>
""", unsafe_allow_html=True)

# Page Navigation
if "page" not in st.session_state:
    st.session_state.page = "form"
if "api_response" not in st.session_state:
    st.session_state.api_response = None

options = dict()

# 📌 RESULT PAGE
elif st.session_state.page == "result":

    if st.button("🔙 Go Back to Form"):
        st.session_state.page = "form"
        st.experimental_rerun()

    st.markdown("""
        <style>
            body { background-color: black !important; color: white !important; }
        </style>
    """, unsafe_allow_html=True)

    st.markdown('<div style="text-align:center; color:#007BFF; font-size:30px; font-weight:bold;">📑 Submission Summary</div>', unsafe_allow_html=True)

    st.markdown(f"**💊 Product Name:** {st.session_state.product_name}")
    st.markdown(f"**📦 Quantity of Medicine:** {st.session_state.quanOfMed}")
    st.markdown(f"**⚡ Power of Drug:** {st.session_state.powerOfDrug}")

    st.markdown("### 📋 Generated Report")
    if st.session_state.api_response:
        response_table = f"""
        <style>
            table {{ border-collapse: collapse; width: 100%; }}
            th, td {{ border: 1px solid black; padding: 8px; text-align: left; }}
            th {{ background-color: #007BFF; color: white; }}
        </style>
        {st.session_state.api_response.replace('<html>', '').replace('</html>', '')}
        """
        st.markdown(response_table, unsafe_allow_html=True)
    else:
        st.warning("⚠️ No response received from API.")

    if st.session_state.ftir_required:
        with st.spinner("📡 Fetching FTIR Data..."):
            ftir_data = chat_with_gpt.get_ftir_from_gpt(st.session_state.product_name)
            st.markdown("### 🔬 FTIR Data")
            st.write(ftir_data)
