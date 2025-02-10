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
st.set_page_config(page_title="QAI Model - AI-Powered Quality Assistance", layout="wide", page_icon="🧪")

# Apply Custom Styles
st.markdown("""
    <style>
        @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap');

        body { 
            font-family: 'Inter', sans-serif;
            background-color: #F4F6F9; /* Light professional pharma background */
            color: #1E3A8A;
        }
        
        /* Header */
        .header {
            background: linear-gradient(90deg, #0B3D91 0%, #1E4D9E 100%);
            padding: 2rem;
            border-radius: 12px;
            box-shadow: 0 4px 20px rgba(11,61,145,0.15);
            color: white;
            text-align: center;
            font-size: 26px;
            font-weight: 700;
        }

        /* Form Containers */
        .form-container {
            background: white;
            padding: 1.8rem;
            border-radius: 12px;
            box-shadow: 0 6px 20px rgba(0,0,0,0.08);
            border: 1px solid #E2E8F0;
            margin-bottom: 1.5rem;
        }

        /* Input Fields */
        .stTextInput > div > div > input,
        .stSelectbox > div > div > select,
        .stTextArea > div > textarea {
            background-color: #FFFFFF;
            border: 2px solid #007BFF;
            border-radius: 10px;
            padding: 12px;
            font-size: 16px;
            color: #1E3A8A;
            transition: all 0.2s ease;
        }

        .stTextInput > div > div > input:focus,
        .stSelectbox > div > div > select:focus,
        .stTextArea > div > textarea:focus {
            border-color: #004085;
            box-shadow: 0 0 6px rgba(11,61,145,0.2);
        }

        /* Buttons */
        .stButton > button {
            background: linear-gradient(90deg, #00B4DB, #0083B0);
            color: white;
            padding: 14px;
            border-radius: 10px;
            font-weight: 600;
            transition: all 0.3s ease;
            border: none;
            width: 100%;
            text-transform: uppercase;
            font-size: 16px;
        }

        .stButton > button:hover {
            background: linear-gradient(90deg, #007BFF, #004085);
            transform: scale(1.05);
            box-shadow: 0 8px 20px rgba(11,61,145,0.2);
        }

        /* Results Table */
        .results-table {
            border-collapse: collapse;
            width: 100%;
            border-radius: 10px;
            box-shadow: 0 4px 15px rgba(0,0,0,0.08);
        }

        .results-table th {
            background: #0B3D91;
            color: white;
            padding: 1rem;
            text-align: left;
            font-weight: 600;
        }

        .results-table td {
            padding: 1rem;
            border-bottom: 1px solid #E2E8F0;
            background: white;
        }

        .results-table tr:hover td {
            background-color: #F8FAFC;
        }

    </style>
""", unsafe_allow_html=True)

# Store session state
if "page" not in st.session_state:
    st.session_state.page = "form"
if "api_response" not in st.session_state:
    st.session_state.api_response = None

options = dict()

# Utility Functions
def get_pubchem_product_code(product_name):
    """Fetches the PubChem Canonical SMILES code for a given product name."""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{product_name}/property/CanonicalSMILES/JSON"
    response = requests.get(url)
    
    if response.status_code == 200:
        try:
            return response.json()["PropertyTable"]["Properties"][0]["CanonicalSMILES"]
        except (KeyError, IndexError):
            return None
    return None

def showStructure(product_name):
    """Generates and returns a molecular structure image from PubChem data."""
    product_code = get_pubchem_product_code(product_name)

    if not product_code:
        product_code_prompt = prompts.STRUCTURE_PROMPT.substitute(product_name=product_name)
        product_code = chat_with_gpt.chatWithGpt(product_code_prompt)

        if product_code == "NO DRUG FOUND":
            return None  

    try:
        m = Chem.MolFromSmiles(product_code)
        if m is not None:
            return Draw.MolToImage(m, size=(250, 250))
    except Exception as e:
        st.error(f"Error generating structure: {str(e)}")
    
    return None  

# 📌 FORM PAGE
if st.session_state.page == "form":
    # Header
    st.markdown('<div class="header">🧪 QAI Model - AI-Powered Quality Assurance</div>', unsafe_allow_html=True)

    # Form Section
    with st.container():
        col1, col2 = st.columns(2)

        with col1:
            st.markdown('<div class="form-container">', unsafe_allow_html=True)
            options["product_name"] = st.text_input("💊 Product Name", placeholder="e.g., Paracetamol")
            
            if st.button("🔬 Get Structure"):
                if options["product_name"]:
                    with st.spinner("Generating molecular structure..."):
                        fig = showStructure(options["product_name"])
                        if fig:
                            st.image(fig, caption=f"{options['product_name']} Structure")
                        else:
                            st.error("⚠️ No valid structure found. Please check the product name.")
                else:
                    st.warning("⚠️ Please enter a product name first.")
                    
            st.markdown('</div>', unsafe_allow_html=True)

        with col2:
            st.markdown('<div class="form-container">', unsafe_allow_html=True)
            options["quanOfMed"] = st.text_input("📦 Quantity", placeholder="e.g., 1000 tablets")
            options["powerOfDrug"] = st.text_input("⚡ Strength", placeholder="e.g., 500 mg")
            options["jurisdiction"] = st.selectbox("🌎 Pharmacopoeia Reference",
                ["INDIAN PHARMACOPIEA", "BRITISH PHARMACOPIEA", "UNITED STATES PHARMACOPOEIA",
                 "MARTINDALE-EXTRA PHARMACOPIEA", "COMPARE WITH ALL"])
            st.markdown('</div>', unsafe_allow_html=True)

    # Submit Button
    if st.button("🚀 Generate Report"):
        if all([options["product_name"], options["quanOfMed"], options["powerOfDrug"]]):
            with st.spinner("Analyzing data and generating report..."):
                prompt = prompts.getPromptForOptions(options)
                api_response = chat_with_gpt.chatWithGpt(prompt)
                st.session_state.api_response = api_response
                st.session_state.page = "result"
                st.experimental_rerun()
        else:
            st.error("⚠️ Please fill in all required fields.")

