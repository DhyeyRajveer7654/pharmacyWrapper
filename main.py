import streamlit as st
import streamlit.components.v1 as components
import prompts
import chat_with_gpt
from string import Template
from rdkit import Chem
from rdkit.Chem import Draw
import requests
import os
from PIL import Image

# Initialize session state
if "page" not in st.session_state:
    st.session_state.page = "form"
if "api_response" not in st.session_state:
    st.session_state.api_response = None
if "structure" not in st.session_state:
    st.session_state.structure = None
if "ftir_graph" not in st.session_state:
    st.session_state.ftir_graph = None

# Ensure set_page_config is at the very beginning
if not hasattr(st, "_is_configured"):
    st.set_page_config(
        page_title="Pharma QA Assistant",
        layout="wide",
        page_icon="🧪",
        initial_sidebar_state="collapsed"
    )
    st._is_configured = True

# Load custom CSS
st.markdown(prompts.TABLE_STYLE, unsafe_allow_html=True)

# Function to get molecular structure
def show_structure(product_name):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{product_name}/property/IsomericSMILES/JSON"
    response = requests.get(url)
    if response.status_code == 200:
        smiles = response.json()['PropertyTable']['Properties'][0]['IsomericSMILES']
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return Draw.MolToImage(mol, size=(400, 400))
    return None

# Function to get FTIR graph
def show_ftir(product_name):
    product_name = product_name.lower()
    image_path = f"attached_assets/{product_name}.png"
    if os.path.exists(image_path):
        return Image.open(image_path)
    return None

# Main Form Page
if st.session_state.page == "form":
    st.markdown("""
        <div class="header-container">
            <h1 class="main-title">🧪 AI-Powered Pharmaceutical QA</h1>
            <p class="subtitle">Advanced Analysis & Reporting System</p>
        </div>
    """, unsafe_allow_html=True)

    col1, col2 = st.columns([3, 2])

    with col1:
        product_name = st.text_input("💊 Product Name", placeholder="e.g., Paracetamol")
        quan_of_med = st.text_input("📦 Quantity", placeholder="e.g., 1000 tablets")
        power_of_drug = st.text_input("⚡ Strength", placeholder="e.g., 500 mg")

        jurisdiction = st.selectbox("🌎 Jurisdiction", ["INDIAN PHARMACOPIEA", "BRITISH PHARMACOPIEA", "UNITED STATES PHARMACOPOEIA", "COMPARE WITH ALL"])
        type_of_info = st.selectbox("📊 Information Required", ["METHOD OF PREPARATION", "CHARACTERIZATION/EVALUATION", "Both", "CHECK RESULTS"])
        ftir_required = st.checkbox("📡 Include FTIR Analysis")

        results_to_check = ""
        if type_of_info == "CHECK RESULTS":
            results_to_check = st.text_area("🔍 Laboratory Results", height=150, placeholder="Enter your lab results here...")

        col_submit, col_ftir, col_struct = st.columns([2, 1, 1])
        with col_submit:
            submit_button = st.button("🚀 Generate Report")
        with col_ftir:
            show_ftir_button = st.button("📊 FTIR Graph")
        with col_struct:
            get_structure = st.button("🔬 Show Structure")

    with col2:
        if st.session_state.ftir_graph is not None:
            st.image(st.session_state.ftir_graph, caption=f"FTIR Graph - {product_name}", use_container_width=True)
        elif st.session_state.structure is not None:
            st.image(st.session_state.structure, caption="Molecular Structure", use_container_width=True)
        else:
            st.info("Enter a product name and click 'Show Structure' or 'FTIR Graph' to view")

    if show_ftir_button and product_name:
        with st.spinner("📊 Loading FTIR graph..."):
            ftir_image = show_ftir(product_name)
            if ftir_image:
                st.session_state.ftir_graph = ftir_image
                st.session_state.structure = None
                st.rerun()
            else:
                st.warning("⚠️ FTIR graph not available for this compound")

    if get_structure and product_name:
        with st.spinner("🔍 Generating structure..."):
            structure = show_structure(product_name)
            if structure:
                st.session_state.structure = structure
                st.session_state.ftir_graph = None
                st.rerun()
            else:
                st.warning("⚠️ Could not generate structure for this compound")

    if submit_button:
        if not all([product_name, quan_of_med, power_of_drug]):
            st.error("⚠️ Please fill in all required fields!")
        else:
            options = {
                "product_name": product_name,
                "quanOfMed": quan_of_med,
                "powerOfDrug": power_of_drug,
                "jurisdiction": jurisdiction,
                "typeOfInfo": type_of_info,
                "ftir_required": ftir_required
            }
            if type_of_info == "CHECK RESULTS":
                options["resultsToCheck"] = results_to_check

            with st.spinner("🔬 Analyzing and generating report..."):
                prompt = prompts.getPromptForOptions(options)
                api_response = chat_with_gpt.chatWithGpt(prompt)
                if api_response:
                    st.session_state.api_response = api_response
                    st.session_state.update(options)
                    st.session_state.page = "result"
                    st.experimental_rerun()
                else:
                    st.error("⚠️ Failed to generate report. Please check API settings.")

# Result Page
elif st.session_state.page == "result":
    if st.button("🔙 Back to Form"):
        st.session_state.page = "form"
        st.experimental_rerun()

    if st.session_state.api_response:
        st.markdown("### 📊 Report Summary")
        st.markdown(f"**💊 Product:** {st.session_state.product_name}")
        st.markdown(f"**📦 Quantity:** {st.session_state.quanOfMed}")
        st.markdown(f"**⚡ Strength:** {st.session_state.powerOfDrug}")

        st.markdown("### 📑 Generated Report")
        st.markdown(prompts.TABLE_STYLE + st.session_state.api_response, unsafe_allow_html=True)

        if st.session_state.get('ftir_required', False):
            with st.spinner("📡 Fetching FTIR Data..."):
                ftir_data = chat_with_gpt.get_ftir_from_gpt(st.session_state.product_name)
                if ftir_data:
                    st.markdown("### 🔬 FTIR Analysis")
                    st.markdown(prompts.TABLE_STYLE + ftir_data, unsafe_allow_html=True)

        if st.session_state.structure is not None:
            st.markdown("### 🧪 Molecular Structure")
            st.image(st.session_state.structure, caption="Molecular Structure", width=400)