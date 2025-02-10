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
            background: #F4F6F9;
            color: #1E293B;
        }

        /* Header Styles */
        .main-header {
            background: linear-gradient(90deg, #0B3D91 0%, #1E4D9E 100%);
            padding: 2rem;
            border-radius: 12px;
            box-shadow: 0 4px 20px rgba(11,61,145,0.15);
            margin-bottom: 2rem;
            color: white;
            text-align: center;
        }

        .main-title {
            font-size: 2.5rem;
            font-weight: 700;
        }

        .subtitle {
            font-size: 1.2rem;
            font-weight: 500;
        }

        /* Form Styling */
        .form-container {
            background: white;
            padding: 2rem;
            border-radius: 12px;
            box-shadow: 0 6px 20px rgba(0,0,0,0.08);
            border: 1px solid #E2E8F0;
            margin-bottom: 1.5rem;
        }

        /* Input Fields */
        input, select, textarea {
            background-color: #FFFFFF;
            border: 2px solid #E2E8F0;
            border-radius: 8px;
            padding: 12px;
            font-size: 1rem;
            transition: all 0.2s ease;
            color: #1E293B;
            width: 100%;
        }

        input:focus, select:focus, textarea:focus {
            border-color: #0B3D91;
            box-shadow: 0 0 0 3px rgba(11,61,145,0.1);
            background-color: #FFFFFF;
        }

        /* Buttons */
        .stButton > button {
            background: linear-gradient(90deg, #0B3D91 0%, #1E4D9E 100%);
            color: white;
            padding: 12px;
            border-radius: 8px;
            font-weight: 600;
            transition: all 0.3s ease;
            width: 100%;
        }

        .stButton > button:hover {
            transform: translateY(-2px);
            box-shadow: 0 8px 20px rgba(11,61,145,0.2);
        }

        /* Table Styling */
        .results-table {
            border-collapse: collapse;
            width: 100%;
            border-radius: 8px;
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

# Store user session state
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
    st.markdown("""
        <div class="main-header">
            <h1 class="main-title">🧬 PharmQA™ Platform</h1>
            <p class="subtitle">Advanced Pharmaceutical Quality Assurance Analysis</p>
        </div>
    """, unsafe_allow_html=True)

    # Form Section
    with st.container():
        col1, col2 = st.columns(2)

        with col1:
            st.markdown('<div class="form-container">', unsafe_allow_html=True)
            options["product_name"] = st.text_input("Product Name", placeholder="e.g., Paracetamol")
            if st.button("Generate Structure"):
                if options["product_name"]:
                    with st.spinner("Generating molecular structure..."):
                        fig = showStructure(options["product_name"])
                        if fig:
                            st.image(fig, caption=f"{options['product_name']} Structure")
                        else:
                            st.error("Unable to generate structure. Please verify the product name.")
                else:
                    st.warning("Please enter a product name first.")
            st.markdown('</div>', unsafe_allow_html=True)

        with col2:
            st.markdown('<div class="form-container">', unsafe_allow_html=True)
            options["quanOfMed"] = st.text_input("Quantity", placeholder="e.g., 1000 tablets")
            options["powerOfDrug"] = st.text_input("Strength", placeholder="e.g., 500 mg")
            options["jurisdiction"] = st.selectbox("Pharmacopoeia Reference",
                ["INDIAN PHARMACOPIEA", "BRITISH PHARMACOPIEA", "UNITED STATES PHARMACOPOEIA",
                 "MARTINDALE-EXTRA PHARMACOPIEA", "COMPARE WITH ALL"])
            st.markdown('</div>', unsafe_allow_html=True)

    # Submit Button
    if st.button("Generate Analysis Report"):
        if all([options["product_name"], options["quanOfMed"], options["powerOfDrug"]]):
            with st.spinner("Analyzing data and generating report..."):
                prompt = prompts.getPromptForOptions(options)
                api_response = chat_with_gpt.chatWithGpt(prompt)
                st.session_state.api_response = api_response
                st.session_state.page = "result"
                st.experimental_rerun()
        else:
            st.error("Please fill in all required fields.")
