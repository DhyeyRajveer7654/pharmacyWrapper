import streamlit as st
import streamlit.components.v1 as components
import prompts
import chat_with_gpt
from string import Template
from rdkit import Chem
from rdkit.Chem import Draw
import requests

# Set Page Configuration
st.set_page_config(page_title="QAI Model", layout="wide", page_icon="🧪")

# Apply Pharma-Based Styles
st.markdown("""
    <style>
        body { background-color: #87CEEB; color: #0b3d91; font-family: 'Arial', sans-serif; }
        .stApp { background-color: #F8F9FA; }
        .input-box {
            padding: 10px;
            border-radius: 10px;
            background: white;
            box-shadow: 2px 2px 10px rgba(0, 0, 0, 0.1);
            margin-bottom: 10px;
            border: 1px solid #0078D4;
        }
        .report-container {
            background: white;
            padding: 15px;
            border-radius: 10px;
            box-shadow: 0px 4px 10px rgba(0,0,0,0.1);
        }
        .title {
            font-size: 24px;
            font-weight: bold;
            color: #0078D4;
            margin-bottom: 20px;
        }
        .stButton>button:hover {
            background: linear-gradient(90deg, #00D4FF, #007BFF);
            transform: scale(1.05);
        }
    </style>
""", unsafe_allow_html=True)

# Page Navigation
if "page" not in st.session_state:
    st.session_state.page = "form"
if "api_response" not in st.session_state:
    st.session_state.api_response = None

options = dict()

def get_pubchem_structure(product_name):
    """Fetches molecular structure from PubChem."""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{product_name}/property/CanonicalSMILES/JSON"
    response = requests.get(url)
    if response.status_code == 200:
        try:
            smiles = response.json()["PropertyTable"]["Properties"][0]["CanonicalSMILES"]
            mol = Chem.MolFromSmiles(smiles)
            return Draw.MolToImage(mol, size=(250, 250))
        except (KeyError, IndexError):
            return None
    return None

# 📌 FORM PAGE
if st.session_state.page == "form":
    st.markdown('<div class="title">🧪 QAI Model - AI-Powered Quality Assurance</div>', unsafe_allow_html=True)

    # Input Form
    col1, col2 = st.columns(2)

    with col1:
        options["product_name"] = st.text_input("💊 Product Name", placeholder="e.g., Paracetamol")
        if st.button("Get Structure"):
            if not options["product_name"]:
                st.error("⚠️ Please enter a product name!")
            else:
                with st.spinner("Fetching structure..."):
                    structure_img = get_pubchem_structure(options["product_name"])
                if structure_img:
                    st.image(structure_img, caption=f"{options['product_name']} Molecule")
                else:
                    st.error("⚠️ Drug structure not found. Please enter a valid drug name.")

    with col2:
        options["quanOfMed"] = st.text_input("📦 Quantity of Medicine", placeholder="e.g., 1000 tablets")
        options["jurisdiction"] = st.selectbox("🌎 Select Jurisdiction", 
            ["INDIAN PHARMACOPIEA", "BRITISH PHARMACOPIEA", "UNITED STATES PHARMACOPOEIA", "MARTINDALE-EXTRA PHARMACOPIEA", "COMPARE WITH ALL"])
        options["powerOfDrug"] = st.text_input("⚡ Power of Drug", placeholder="e.g., 500 mg")

    options["typeOfInfo"] = st.selectbox("📊 Select Information Required:", 
            ["METHOD OF PREPARATION", "CHARACTARIZATION/EVALUATION", "Both of above", "CHECK RESULTS"])

    if options["typeOfInfo"] == "CHECK RESULTS":
        options["resultsToCheck"] = st.text_area("🔍 Enter Your Results:", height=200, placeholder="Paste lab results here...", key="checkResults")

    options["ftir_required"] = st.checkbox("📡 Retrieve FTIR Data")

    submit_button = st.button("🚀 Submit & Generate Report")
    if submit_button:
        if not all([options["product_name"], options["quanOfMed"], options["powerOfDrug"]]):
            st.error("⚠️ Please fill in all required fields!")
        else:
            prompt = prompts.getPromptForOptions(options)
            with st.spinner("🛠️ Processing... Please wait"):
                api_response = chat_with_gpt.chatWithGpt(prompt)
                st.session_state.api_response = api_response
            st.session_state.page = "result"
            st.experimental_rerun()

# 📌 RESULT PAGE
elif st.session_state.page == "result":
    if st.button("🔙 Go Back to Form"):
        st.session_state.page = "form"
        st.experimental_rerun()

    # Display submission summary
    st.markdown('<div class="title">📑 Submission Summary</div>', unsafe_allow_html=True)
    st.markdown(f"**💊 Product Name:** {st.session_state.product_name}")
    st.markdown(f"**📦 Quantity of Medicine:** {st.session_state.quanOfMed}")
    st.markdown(f"**⚡ Power of Drug:** {st.session_state.powerOfDrug}")

    # Display report results in a table
    st.markdown("### 📋 Generated Report")
    if st.session_state.api_response:
        st.markdown(f"<div class='report-container'>{st.session_state.api_response}</div>", unsafe_allow_html=True)
    else:
        st.warning("⚠️ No response received from API.")

    # FTIR Data Retrieval
    if st.session_state.ftir_required:
        with st.spinner("📡 Fetching FTIR Data..."):
            ftir_data = chat_with_gpt.get_ftir_from_gpt(st.session_state.product_name)
            st.markdown("### 🔬 FTIR Data")
            st.write(ftir_data)
