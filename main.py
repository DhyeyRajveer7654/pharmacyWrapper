import streamlit as st
import streamlit.components.v1 as components
import prompts
import chat_with_gpt
from rdkit import Chem
from rdkit.Chem import Draw
import requests

st.set_page_config(page_title="QAI Model", layout="wide", page_icon="🧪")

# **Apply Pharma UI**
st.markdown("""
    <style>
        body { background-color: #f8fbff; color: #0b3d91; font-family: 'Arial', sans-serif; }
        .stButton>button { background: linear-gradient(90deg, #004080, #007BFF); color: white; border-radius: 10px; }
        .title { color: #004080; text-align: center; font-size: 32px; font-weight: bold; }
        .subtitle { color: #0b3d91; text-align: center; font-size: 18px; }
    </style>
""", unsafe_allow_html=True)

if "page" not in st.session_state:
    st.session_state.page = "form"

options = {}

# **STRUCTURE FEATURE**
def get_structure(product_name):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{product_name}/property/CanonicalSMILES/JSON"
    response = requests.get(url)
    if response.status_code == 200:
        smiles = response.json().get("PropertyTable", {}).get("Properties", [{}])[0].get("CanonicalSMILES", "NO DATA")
        if smiles != "NO DATA":
            mol = Chem.MolFromSmiles(smiles)
            return Draw.MolToImage(mol, size=(250, 250))
    return None

# **FORM PAGE**
if st.session_state.page == "form":
    st.markdown('<div class="title">🧪 QAI Model - AI-Powered Pharma QA</div>', unsafe_allow_html=True)

    options["product_name"] = st.text_input("💊 Product Name", placeholder="e.g., Paracetamol")
    if st.button("Get Structure"):
        if options["product_name"]:
            with st.spinner("Fetching Structure..."):
                structure_img = get_structure(options["product_name"])
                if structure_img:
                    st.image(structure_img, caption="Molecular Structure")
                else:
                    st.error("Structure not found.")

    options["quanOfMed"] = st.text_input("📦 Quantity", placeholder="e.g., 1000 tablets")
    options["powerOfDrug"] = st.text_input("⚡ Strength", placeholder="e.g., 500 mg")
    options["jurisdiction"] = st.selectbox("🌍 Jurisdiction", ["INDIAN PHARMACOPIEA", "BP", "USP", "COMPARE ALL"])
    options["typeOfInfo"] = st.radio("📊 Select Report Type:", ["METHOD OF PREPARATION", "CHARACTERIZATION/EVALUATION", "CHECK RESULTS"])
    options["ftir_required"] = st.checkbox("🔬 Retrieve FTIR Data")

    if st.button("🚀 Generate Report"):
        if not all([options["product_name"], options["quanOfMed"], options["powerOfDrug"]]):
            st.error("⚠️ Fill all required fields!")
        else:
            prompt = prompts.getPromptForOptions(options)
            with st.spinner("Generating report..."):
                st.session_state.api_response = chat_with_gpt.chatWithGpt(prompt)
            st.session_state.page = "result"
            st.experimental_rerun()

# **RESULT PAGE**
elif st.session_state.page == "result":
    st.markdown('<div class="title">📑 Generated Report</div>', unsafe_allow_html=True)
    
    if st.session_state.api_response:
        components.html(st.session_state.api_response, height=700, scrolling=True)
    else:
        st.warning("⚠️ No response received.")

    if st.session_state.ftir_required:
        with st.spinner("Fetching FTIR Data..."):
            ftir_data = chat_with_gpt.get_ftir_from_gpt(st.session_state.product_name)
            st.markdown("### 🔬 FTIR Data")
            st.write(ftir_data)

    if st.button("🔙 Back"):
        st.session_state.page = "form"
        st.experimental_rerun()
