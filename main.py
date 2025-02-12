import streamlit as st
import rdkit
from rdkit import Chem
from rdkit.Chem import Draw
import requests
import prompts
import chat_with_gpt

# Initialize session state
if "page" not in st.session_state:
    st.session_state.page = "form"
if "api_response" not in st.session_state:
    st.session_state.api_response = None
if "structure" not in st.session_state:
    st.session_state.structure = None

# Set Page Configuration (MUST BE FIRST)
st.set_page_config(
    page_title="Pharma QA Assistant",
    layout="wide",
    page_icon="🧪",
    initial_sidebar_state="collapsed"
)

# Load custom CSS
try:
    with open('static/custom.css') as f:
        st.markdown(f'<style>{f.read()}</style>', unsafe_allow_html=True)
except Exception as e:
    st.warning("Custom CSS file not found. Using default styling.")

# Additional styling for FTIR checkbox
st.markdown("""
    <style>
    div[data-testid="stCheckbox"] {
        background-color: rgba(255, 255, 255, 0.1);
        padding: 10px;
        border-radius: 5px;
        border: 1px solid rgba(255, 255, 255, 0.2);
    }
    div[data-testid="stCheckbox"] label {
        color: #ffffff !important;
        font-weight: 500;
    }
    </style>
""", unsafe_allow_html=True)

def show_structure(compound_name):
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{compound_name}/property/IsomericSMILES/JSON"
        response = requests.get(url)

        if response.status_code == 200:
            smiles = response.json()['PropertyTable']['Properties'][0]['IsomericSMILES']
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                img = Draw.MolToImage(mol, size=(400, 400))
                return img
            return None
        return None
    except Exception as e:
        st.error(f"Error generating structure: {str(e)}")
        return None

# Main Form Page
if st.session_state.page == "form":
    st.markdown('''
        <div class="header-container">
            <h1 class="main-title">🧪 AI-Powered Pharmaceutical QA</h1>
            <p class="subtitle">Advanced Analysis & Reporting System</p>
        </div>
    ''', unsafe_allow_html=True)

    col1, col2 = st.columns([3, 2])

    with col1:
        with st.form("input_form"):
            product_name = st.text_input("💊 Product Name", placeholder="e.g., Paracetamol")

            col_a, col_b = st.columns(2)
            with col_a:
                quan_of_med = st.text_input("📦 Quantity", placeholder="e.g., 1000 tablets")
            with col_b:
                power_of_drug = st.text_input("⚡ Strength", placeholder="e.g., 500 mg")

            jurisdiction = st.selectbox(
                "🌎 Jurisdiction",
                ["INDIAN PHARMACOPIEA", "BRITISH PHARMACOPIEA", "UNITED STATES PHARMACOPOEIA", 
                 "MARTINDALE-EXTRA PHARMACOPIEA", "COMPARE WITH ALL"]
            )

            type_of_info = st.selectbox(
                "📊 Information Required",
                ["METHOD OF PREPARATION", "CHARACTARIZATION/EVALUATION", 
                 "Both of above", "CHECK RESULTS"]
            )

            if type_of_info == "CHECK RESULTS":
                results_to_check = st.text_area(
                    "🔍 Laboratory Results",
                    height=150,
                    placeholder="Enter your lab results here..."
                )

            ftir_required = st.checkbox("📡 Include FTIR Analysis", help="Check to include FTIR spectral analysis in the report")

            col_submit, col_struct = st.columns([3, 1])
            with col_submit:
                submit_button = st.form_submit_button("🚀 Generate Report")
            with col_struct:
                get_structure = st.form_submit_button("🔬 Show Structure")

            if get_structure and product_name:
                with st.spinner("🔍 Generating structure..."):
                    structure = show_structure(product_name)
                    if structure:
                        st.session_state.structure = structure
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
                            st.rerun()
                        else:
                            st.error("⚠️ Failed to generate report. Please check API configuration.")

    with col2:
        if st.session_state.structure is not None:
            st.image(st.session_state.structure, caption="Molecular Structure", use_column_width=True)
        else:
            st.info("Enter a product name and click 'Show Structure' to view the molecular structure")

# Result Page
elif st.session_state.page == "result":
    if st.button("🔙 Back to Form"):
        st.session_state.page = "form"
        st.rerun()

    if st.session_state.api_response:
        st.markdown("### 📊 Report Summary")
        st.markdown(f"""
        <div class="report-summary">
            <p><strong>💊 Product:</strong> {st.session_state.product_name}</p>
            <p><strong>📦 Quantity:</strong> {st.session_state.quanOfMed}</p>
            <p><strong>⚡ Strength:</strong> {st.session_state.powerOfDrug}</p>
        </div>
        """, unsafe_allow_html=True)

        st.markdown("### 📑 Generated Report")
        st.markdown(st.session_state.api_response, unsafe_allow_html=True)

        if st.session_state.get('ftir_required', False):
            with st.spinner("📡 Fetching FTIR Data..."):
                ftir_data = chat_with_gpt.get_ftir_from_gpt(st.session_state.product_name)
                if ftir_data:
                    st.markdown("### 🔬 FTIR Analysis")
                    st.markdown(ftir_data, unsafe_allow_html=True)

        if st.session_state.structure is not None:
            st.markdown("### 🧪 Molecular Structure")
            st.image(st.session_state.structure, caption="Molecular Structure", width=400)