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
    page_title="PharmQA™ - Quality Assurance Platform",
    layout="wide",
    page_icon="🧬"
)

# Apply Custom Styles
st.markdown("""
    <style>
        @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap');
        
        /* Global Styles */
        body {
            font-family: 'Inter', sans-serif;
            background: linear-gradient(135deg, #f8f9fa 0%, #e9ecef 100%);
            color: #1E293B;
        }
        
        /* Header Styles */
        .main-header {
            background: linear-gradient(90deg, #0B3D91 0%, #1E4D9E 100%);
            padding: 2rem;
            border-radius: 0 0 20px 20px;
            box-shadow: 0 4px 12px rgba(0,0,0,0.1);
            margin-bottom: 2rem;
            color: white;
            text-align: center;
        }
        
        .main-title {
            font-size: 2.5rem;
            font-weight: 700;
            margin-bottom: 0.5rem;
            background: linear-gradient(90deg, #FFFFFF 0%, #E3F2FD 100%);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
        }
        
        .subtitle {
            font-size: 1.1rem;
            opacity: 0.9;
            font-weight: 500;
        }

        /* Form Card Styles */
        .form-card {
            background: white;
            padding: 2rem;
            border-radius: 16px;
            box-shadow: 0 4px 20px rgba(0,0,0,0.08);
            border: 1px solid #E2E8F0;
            margin-bottom: 1.5rem;
            transition: transform 0.2s ease;
        }
        
        .form-card:hover {
            transform: translateY(-2px);
        }

        /* Input Styling */
        .stTextInput > div > div > input,
        .stSelectbox > div > div > select,
        .stTextArea > div > textarea {
            background-color: #F8FAFC;
            border: 2px solid #E2E8F0;
            border-radius: 10px;
            padding: 0.75rem 1rem;
            font-size: 1rem;
            transition: all 0.2s ease;
            color: #1E293B;
        }

        .stTextInput > div > div > input:focus,
        .stSelectbox > div > div > select:focus,
        .stTextArea > div > textarea:focus {
            border-color: #0B3D91;
            box-shadow: 0 0 0 2px rgba(11,61,145,0.1);
        }

        /* Button Styling */
        .stButton > button {
            background: linear-gradient(90deg, #0B3D91 0%, #1E4D9E 100%);
            color: white;
            padding: 0.75rem 2rem;
            border-radius: 10px;
            border: none;
            font-weight: 600;
            transition: all 0.3s ease;
            width: 100%;
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }

        .stButton > button:hover {
            transform: translateY(-2px);
            box-shadow: 0 4px 12px rgba(11,61,145,0.2);
        }

        /* Results Table Styling */
        .results-table {
            border-collapse: separate;
            border-spacing: 0;
            width: 100%;
            border-radius: 12px;
            overflow: hidden;
            box-shadow: 0 4px 12px rgba(0,0,0,0.05);
        }

        .results-table th {
            background: #0B3D91;
            color: white;
            padding: 1rem;
            font-weight: 600;
            text-align: left;
        }

        .results-table td {
            padding: 1rem;
            border-bottom: 1px solid #E2E8F0;
            background: white;
        }

        .results-table tr:last-child td {
            border-bottom: none;
        }

        /* Loading Animation */
        .stSpinner {
            border: 4px solid #f3f3f3;
            border-top: 4px solid #0B3D91;
            border-radius: 50%;
            width: 40px;
            height: 40px;
            animation: spin 1s linear infinite;
        }

        @keyframes spin {
            0% { transform: rotate(0deg); }
            100% { transform: rotate(360deg); }
        }

        /* Alert Styling */
        .stAlert {
            border-radius: 10px;
            border: none;
            padding: 1rem;
            margin: 1rem 0;
        }

        /* Molecule Visualization */
        .molecule-card {
            background: white;
            padding: 1.5rem;
            border-radius: 12px;
            box-shadow: 0 4px 12px rgba(0,0,0,0.05);
            text-align: center;
        }

        .molecule-title {
            font-weight: 600;
            margin-bottom: 1rem;
            color: #0B3D91;
        }
    </style>
""", unsafe_allow_html=True)

# Initialize session state
if "page" not in st.session_state:
    st.session_state.page = "form"
if "api_response" not in st.session_state:
    st.session_state.api_response = None

options = dict()

# Utility Functions
def get_pubchem_product_code(product_name):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{product_name}/property/CanonicalSMILES/JSON"
    response = requests.get(url)
    
    if response.status_code == 200:
        try:
            return response.json()["PropertyTable"]["Properties"][0]["CanonicalSMILES"]
        except (KeyError, IndexError):
            return None
    return None

def showStructure(product_name):
    product_code = get_pubchem_product_code(product_name)
    
    if not product_code:
        product_code_prompt = prompts.STRUCTURE_PROMPT.substitute(product_name=product_name)
        product_code = chat_with_gpt.chatWithGpt(product_code_prompt)
        
        if product_code == "NO DRUG FOUND":
            return None

    try:
        m = Chem.MolFromSmiles(product_code)
        if m is not None:
            return Draw.MolToImage(m, size=(300, 300))
    except Exception as e:
        st.error(f"Error generating structure: {str(e)}")
    return None

# Main UI
if st.session_state.page == "form":
    # Header
    st.markdown('''
        <div class="main-header">
            <h1 class="main-title">🧬 PharmQA™ Platform</h1>
            <p class="subtitle">Advanced Pharmaceutical Quality Assurance Analysis</p>
        </div>
    ''', unsafe_allow_html=True)

    # Form Section
    with st.container():
        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown('<div class="form-card">', unsafe_allow_html=True)
            options["product_name"] = st.text_input("Product Name", placeholder="e.g., Paracetamol")
            if st.button("Generate Structure"):
                if options["product_name"]:
                    with st.spinner("Generating molecular structure..."):
                        fig = showStructure(options["product_name"])
                        if fig:
                            st.markdown('<div class="molecule-card">', unsafe_allow_html=True)
                            st.markdown(f'<p class="molecule-title">{options["product_name"]} Structure</p>', unsafe_allow_html=True)
                            st.image(fig)
                            st.markdown('</div>', unsafe_allow_html=True)
                        else:
                            st.error("Unable to generate structure. Please verify the product name.")
                else:
                    st.warning("Please enter a product name first.")
            st.markdown('</div>', unsafe_allow_html=True)

        with col2:
            st.markdown('<div class="form-card">', unsafe_allow_html=True)
            options["quanOfMed"] = st.text_input("Quantity", placeholder="e.g., 1000 tablets")
            options["powerOfDrug"] = st.text_input("Strength", placeholder="e.g., 500 mg")
            options["jurisdiction"] = st.selectbox("Pharmacopoeia Reference",
                ["INDIAN PHARMACOPIEA", "BRITISH PHARMACOPIEA", "UNITED STATES PHARMACOPOEIA",
                 "MARTINDALE-EXTRA PHARMACOPIEA", "COMPARE WITH ALL"])
            st.markdown('</div>', unsafe_allow_html=True)

    # Analysis Options
    st.markdown('<div class="form-card">', unsafe_allow_html=True)
    options["typeOfInfo"] = st.selectbox("Analysis Type",
        ["METHOD OF PREPARATION", "CHARACTARIZATION/EVALUATION", "Both of above", "CHECK RESULTS"])

    if options["typeOfInfo"] == "CHECK RESULTS":
        options["resultsToCheck"] = st.text_area("Laboratory Results", 
            height=200, placeholder="Enter your laboratory results here...")

    options["ftir_required"] = st.checkbox("Include FTIR Analysis")
    
    if st.button("Generate Analysis Report"):
        if all([options["product_name"], options["quanOfMed"], options["powerOfDrug"]]):
            with st.spinner("Analyzing data and generating report..."):
                prompt = prompts.getPromptForOptions(options)
                api_response = chat_with_gpt.chatWithGpt(prompt)
                st.session_state.api_response = api_response
                st.session_state.update(options)
                st.session_state.page = "result"
                st.experimental_rerun()
        else:
            st.error("Please fill in all required fields.")
    st.markdown('</div>', unsafe_allow_html=True)

# Results Page
elif st.session_state.page == "result":
    if st.button("← Back to Analysis"):
        st.session_state.page = "form"
        st.experimental_rerun()

    st.markdown('''
        <div class="main-header">
            <h1 class="main-title">Analysis Report</h1>
            <p class="subtitle">Comprehensive Quality Analysis Results</p>
        </div>
    ''', unsafe_allow_html=True)

    # Display Results
    if st.session_state.api_response:
        st.markdown('<div class="form-card">', unsafe_allow_html=True)
        st.markdown(f"""
            <div class="results-table">
                {st.session_state.api_response}
            </div>
        """, unsafe_allow_html=True)
        st.markdown('</div>', unsafe_allow_html=True)

        if st.session_state.ftir_required:
            with st.spinner("Retrieving FTIR data..."):
                ftir_data = chat_with_gpt.get_ftir_from_gpt(st.session_state.product_name)
                if ftir_data:
                    st.markdown('<div class="form-card">', unsafe_allow_html=True)
                    st.markdown("### FTIR Analysis")
                    st.markdown(ftir_data, unsafe_allow_html=True)
                    st.markdown('</div>', unsafe_allow_html=True)
    else:
        st.error("No analysis results available. Please try again.")
