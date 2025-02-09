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
        if st.button("Get structure"):
            if not options["product_name"]:
                st.error("⚠️ Please write product name!")
            else:
                with st.spinner("🛠️ Processing... Please wait"):
                    fig = showStructure(options["product_name"])
                if not fig:
                    st.error("⚠️ Drug not found, please input a valid drug name")
                else:
                    st.image(fig, caption=f"{options['product_name']} Molecule")

    with col2:
        options["quanOfMed"] = st.text_input("📦 Quantity of Medicine", placeholder="e.g., 1000 tablets")
        options["jurisdiction"] = st.selectbox("🌎 Select Jurisdiction", 
            ["INDIAN PHARMACOPIEA", "BRITISH PHARMACOPIEA", "UNITED STATES PHARMACOPOEIA", "MARTINDALE-EXTRA PHARMACOPIEA", "COMPARE WITH ALL"])
        options["powerOfDrug"] = st.text_input("⚡ Power of Drug", placeholder="e.g., 500 mg")

    options["typeOfInfo"] = st.selectbox("📊 Select Information Required:", 
            ["METHOD OF PREPARATION", "CHARACTARIZATION/EVALUATION", "Both of above", "CHECK RESULTS"])

    if options["typeOfInfo"] == "CHECK RESULTS":
        options["resultsToCheck"] = st.text_area("🔍 Enter Your Results:", height=200, placeholder="Paste lab results here...", key="checkResults")

    options["ftir_required"] = st.checkbox("📱 Retrieve FTIR Data")

    submit_button = st.button("🚀 Submit & Generate Report")
    if submit_button:
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

    st.markdown('<div style="text-align:center; color:#007BFF; font-size:30px; font-weight:bold;">📝 Submission Summary</div>', unsafe_allow_html=True)
    st.markdown(f"**💊 Product Name:** {st.session_state.product_name}")
    st.markdown(f"**📦 Quantity of Medicine:** {st.session_state.quanOfMed}")
    st.markdown(f"**⚡ Power of Drug:** {st.session_state.powerOfDrug}")

    st.markdown("### 📋 Generated Report")
    if st.session_state.api_response:
        components.html(st.session_state.api_response, height=800, scrolling=True)
    else:
        st.warning("⚠️ No response received from API.")

    if st.session_state.ftir_required:
        with st.spinner("📱 Fetching FTIR Data..."):
            ftir_data = chat_with_gpt.get_ftir_from_gpt(st.session_state.product_name)
            st.markdown("### 🔬 FTIR Data")
            st.write(ftir_data)
