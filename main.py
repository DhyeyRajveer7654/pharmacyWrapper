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

# Apply Custom Styles for Professional & Attractive Look
st.markdown("""
    <style>
        body { background: linear-gradient(135deg, #1a2a6c, #b21f1f, #fdbb2d); color: white; font-family: 'Arial', sans-serif; }
        .stTextInput>div>div>input, .stSelectbox>div>div>select, .stTextArea>div>textarea { 
            border-radius: 8px !important; padding: 12px;
            background-color: #ffffff; color: #000000; border: 2px solid #007BFF;
        }
        .stButton>button { 
            background: linear-gradient(90deg, #007BFF, #00D4FF);
            color: white; font-weight: bold; padding: 10px 20px; border-radius: 8px;
        }
        .stButton>button:hover { 
            background: linear-gradient(90deg, #00D4FF, #007BFF);
            transform: scale(1.08);
        }
        .title { color: #ffffff; text-align: center; font-size: 32px; font-weight: bold; text-shadow: 2px 2px 8px #000; }
        .subtitle { color: #eeeeee; text-align: center; font-size: 20px; text-shadow: 1px 1px 6px #000; }
        .card {
            background-color: rgba(255, 255, 255, 0.1); padding: 20px; border-radius: 15px; 
            box-shadow: 0px 4px 15px rgba(255, 255, 255, 0.2); margin: 20px;
        }
        table { width: 100%; border-collapse: collapse; color: white; background-color: rgba(255,255,255,0.1); }
        th, td { border: 1px solid #ddd; padding: 10px; text-align: center; }
        th { background-color: #007BFF; color: white; }
        tr:nth-child(even) { background-color: rgba(255,255,255,0.2); }
        tr:hover { background-color: #007BFF; color: white; }
    </style>
""", unsafe_allow_html=True)

# Page Navigation
if "page" not in st.session_state:
    st.session_state.page = "form"
if "api_response" not in st.session_state:
    st.session_state.api_response = None

options = dict()

def showStructure(product_name):
    product_code_prompt = prompts.STRUCTURE_PROMPT.substitute(product_name=product_name)
    product_code = chat_with_gpt.chatWithGpt(product_code_prompt)
    if product_code == "NO DRUG FOUND":
        return ""
    m = Chem.MolFromSmiles(product_code)
    fig = Draw.MolToImage(m, size=size)
    return fig
    
# 📌 FORM PAGE
if st.session_state.page == "form":
    st.markdown('<div class="title">🧪 QAI Model - AI-Powered Quality Assurance</div>', unsafe_allow_html=True)
    st.markdown('<div class="subtitle">🔍 Enter details below to generate a pharmaceutical quality report.</div>', unsafe_allow_html=True)

    col1, col2 = st.columns(2)
    with col1:
        options["product_name"] = st.text_input("💊 Product Name", placeholder="e.g., Paracetamol")
    with col2:
        options["quanOfMed"] = st.text_input("📦 Quantity of Medicine", placeholder="e.g., 1000 tablets")
    options["jurisdiction"] = st.selectbox("🌎 Select Jurisdiction", 
            ["INDIAN PHARMACOPIEA", "BRITISH PHARMACOPIEA", "UNITED STATES PHARMACOPOEIA", "MARTINDALE-EXTRA PHARMACOPIEA", "COMPARE WITH ALL"])
    options["powerOfDrug"] = st.text_input("⚡ Power of Drug", placeholder="e.g., 500 mg")
    options["typeOfInfo"] = st.selectbox("📊 Select Information Required:", 
            ["METHOD OF PREPARATION", "CHARACTARIZATION/EVALUATION", "Both of above", "CHECK RESULTS"])
    options["ftir_required"] = st.checkbox("📱 Retrieve FTIR Data")

    if st.button("🚀 Submit & Generate Report"):
        if not all([options["product_name"], options["quanOfMed"], options["powerOfDrug"]]):
            st.error("⚠️ Please fill in all required fields!")
        else:
            prompt = prompts.getPromptForOptions(options)
            with st.spinner("🛠️ Processing... Please wait"):
                api_response = chat_with_gpt.chatWithGpt(prompt)
                st.session_state.api_response = api_response

            st.session_state.update(options)
            st.session_state.page = "result"
            st.experimental_rerun()

# 📌 RESULT PAGE
elif st.session_state.page == "result":
    if st.button("🔙 Go Back to Form"):
        st.session_state.page = "form"
        st.experimental_rerun()
    st.markdown("### 📋 Generated Report")
    if st.session_state.api_response:
        components.html(st.session_state.api_response, height=800, scrolling=True)
    else:
        st.warning("⚠️ No response received from API.")
