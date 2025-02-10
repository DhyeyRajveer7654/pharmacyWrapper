import streamlit as st
import streamlit.components.v1 as components
import prompts
import chat_with_gpt
from string import Template
from rdkit import Chem
from rdkit.Chem import Draw
import requests

# Set Page Configuration
st.set_page_config(
    page_title="QAI Model - Pharmaceutical Quality Assurance",
    layout="wide",
    page_icon="🧪",
    initial_sidebar_state="expanded"
)

# Apply Professional Styling
st.markdown("""
<style>
    /* Global Styles */
    body {
        background-color: #F8F9FA;
        color: #0B3D91;
        font-family: 'Helvetica Neue', Arial, sans-serif;
    }
    
    /* Header Styling */
    .main-header {
        background: linear-gradient(135deg, #0B3D91, #1E88E5);
        color: white;
        padding: 2rem;
        border-radius: 15px;
        margin-bottom: 2rem;
        text-align: center;
        box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
    }
    
    .main-title {
        font-size: 2.5rem;
        font-weight: 700;
        margin-bottom: 1rem;
        text-transform: uppercase;
        letter-spacing: 2px;
    }
    
    .subtitle {
        font-size: 1.2rem;
        opacity: 0.9;
        font-weight: 300;
    }
    
    /* Form Elements */
    .stTextInput>div>div>input,
    .stSelectbox>div>div>select,
    .stTextArea>div>textarea {
        border-radius: 10px !important;
        border: 2px solid #E3F2FD !important;
        padding: 1rem !important;
        background-color: white !important;
        transition: all 0.3s ease !important;
    }
    
    .stTextInput>div>div>input:focus,
    .stSelectbox>div>div>select:focus,
    .stTextArea>div>textarea:focus {
        border-color: #1E88E5 !important;
        box-shadow: 0 0 0 2px rgba(30, 136, 229, 0.2) !important;
    }
    
    /* Button Styling */
    .stButton>button {
        background: linear-gradient(135deg, #0B3D91, #1E88E5);
        color: white;
        border: none;
        padding: 0.75rem 2rem;
        border-radius: 10px;
        font-weight: 600;
        transition: all 0.3s ease;
        width: 100%;
        text-transform: uppercase;
        letter-spacing: 1px;
    }
    
    .stButton>button:hover {
        transform: translateY(-2px);
        box-shadow: 0 4px 12px rgba(30, 136, 229, 0.3);
    }
    
    /* Card Styling */
    .card {
        background: white;
        border-radius: 15px;
        padding: 2rem;
        margin: 1rem 0;
        box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
        border-left: 5px solid #0B3D91;
    }
    
    /* Results Table */
    .results-table {
        width: 100%;
        border-collapse: collapse;
        margin: 1rem 0;
        background: white;
        border-radius: 10px;
        overflow: hidden;
        box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
    }
    
    .results-table th {
        background: #0B3D91;
        color: white;
        padding: 1rem;
        text-align: left;
    }
    
    .results-table td {
        padding: 1rem;
        border-bottom: 1px solid #E3F2FD;
    }
    
    .results-table tr:hover {
        background-color: #F8F9FA;
    }
    
    /* Loading Spinner */
    .stSpinner {
        text-align: center;
        padding: 2rem;
    }
    
    /* Error Messages */
    .error-message {
        background-color: #FFE0E0;
        color: #D32F2F;
        padding: 1rem;
        border-radius: 10px;
        margin: 1rem 0;
        border-left: 5px solid #D32F2F;
    }
</style>
""", unsafe_allow_html=True)

# Initialize Session State
if "page" not in st.session_state:
    st.session_state.page = "form"
if "api_response" not in st.session_state:
    st.session_state.api_response = None

# Helper Functions
def get_cid_from_name(drug_name):
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{drug_name}/cids/JSON"
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            return response.json()["IdentifierList"]["CID"][0]
        return None
    except Exception:
        return None

def get_pubchem_product_code(product_name):
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{product_name}/property/CanonicalSMILES/JSON"
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            return response.json()["PropertyTable"]["Properties"][0]["CanonicalSMILES"]
        return None
    except Exception:
        return None

def showStructure(product_name):
    if not product_name:
        st.error("Please enter a product name first.")
        return None
        
    with st.spinner("Fetching molecular structure..."):
        product_code = get_pubchem_product_code(product_name)
        if not product_code:
            product_code_prompt = prompts.STRUCTURE_PROMPT.substitute(product_name=product_name)
            product_code = chat_with_gpt.chatWithGpt(product_code_prompt)
            
            if product_code == "NO DRUG FOUND":
                st.error("No molecular structure found for this drug.")
                return None

        try:
            m = Chem.MolFromSmiles(product_code)
            if m is not None:
                return Draw.MolToImage(m, size=(300, 300))
        except Exception as e:
            st.error(f"Error generating structure: {str(e)}")
        
        return None

# Main Application Logic
if st.session_state.page == "form":
    # Header
    st.markdown("""
        <div class="main-header">
            <h1 class="main-title">🧪 QAI Model</h1>
            <p class="subtitle">AI-Powered Pharmaceutical Quality Assurance</p>
        </div>
    """, unsafe_allow_html=True)

    # Form Layout
    col1, col2 = st.columns([2, 1])
    
    with col1:
        with st.container():
            st.markdown('<div class="card">', unsafe_allow_html=True)
            options = {}
            options["product_name"] = st.text_input("💊 Product Name", 
                placeholder="e.g., Paracetamol",
                help="Enter the name of the pharmaceutical product")
                
            if st.button("View Molecular Structure", key="structure_button"):
                fig = showStructure(options["product_name"])
                if fig:
                    st.image(fig, caption=f"{options['product_name']} Molecular Structure")
            
            options["quanOfMed"] = st.text_input("📦 Quantity of Medicine",
                placeholder="e.g., 1000 tablets",
                help="Specify the total quantity to be manufactured")
                
            options["powerOfDrug"] = st.text_input("⚡ Power of Drug",
                placeholder="e.g., 500 mg",
                help="Enter the strength of each unit")
            st.markdown('</div>', unsafe_allow_html=True)

    with col2:
        with st.container():
            st.markdown('<div class="card">', unsafe_allow_html=True)
            options["jurisdiction"] = st.selectbox("🌎 Select Jurisdiction",
                ["INDIAN PHARMACOPIEA",
                 "BRITISH PHARMACOPIEA",
                 "UNITED STATES PHARMACOPOEIA",
                 "MARTINDALE-EXTRA PHARMACOPIEA",
                 "COMPARE WITH ALL"])
                 
            options["typeOfInfo"] = st.selectbox("📊 Select Information Required:",
                ["METHOD OF PREPARATION",
                 "CHARACTARIZATION/EVALUATION",
                 "Both of above",
                 "CHECK RESULTS"])

            if options["typeOfInfo"] == "CHECK RESULTS":
                options["resultsToCheck"] = st.text_area("🔍 Enter Your Results:",
                    height=200,
                    placeholder="Paste lab results here...",
                    key="checkResults")

            options["ftir_required"] = st.checkbox("📡 Include FTIR Data Analysis")
            st.markdown('</div>', unsafe_allow_html=True)

    # Submit Button
    if st.button("🚀 Generate Quality Report", type="primary"):
        if not all([options["product_name"], options["quanOfMed"], options["powerOfDrug"]]):
            st.error("⚠️ Please fill in all required fields!")
        else:
            with st.spinner("🔄 Generating comprehensive report..."):
                prompt = prompts.getPromptForOptions(options)
                api_response = chat_with_gpt.chatWithGpt(prompt)
                
                if "Error:" in api_response:
                    st.error(api_response)
                else:
                    st.session_state.api_response = api_response
                    st.session_state.update(options)
                    st.session_state.page = "result"
                    st.experimental_rerun()

elif st.session_state.page == "result":
    # Results Page Header
    st.markdown("""
        <div class="main-header">
            <h1 class="main-title">📊 Quality Analysis Report</h1>
            <p class="subtitle">Comprehensive Evaluation Results</p>
        </div>
    """, unsafe_allow_html=True)

    if st.button("← Back to Form"):
        st.session_state.page = "form"
        st.experimental_rerun()

    # Display Results
    with st.container():
        st.markdown('<div class="card">', unsafe_allow_html=True)
        st.markdown(f"**Product:** {st.session_state.product_name}")
        st.markdown(f"**Quantity:** {st.session_state.quanOfMed}")
        st.markdown(f"**Strength:** {st.session_state.powerOfDrug}")
        st.markdown('</div>', unsafe_allow_html=True)

        if st.session_state.api_response:
            st.markdown(st.session_state.api_response, unsafe_allow_html=True)

            if st.session_state.ftir_required:
                with st.spinner("📡 Analyzing FTIR Data..."):
                    ftir_data = chat_with_gpt.get_ftir_from_gpt(st.session_state.product_name)
                    if "Error:" not in ftir_data:
                        st.markdown("### 🔬 FTIR Analysis")
                        st.markdown(ftir_data, unsafe_allow_html=True)
                    else:
                        st.error(ftir_data)
        else:
            st.error("⚠️ No response received. Please try again.")

## file_path: .streamlit/config.toml
