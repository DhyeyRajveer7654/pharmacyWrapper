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
        @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;600&display=swap');

        /* Global Styles */
        body {
            background: linear-gradient(135deg, #0A192F 0%, #001F3F 100%);
            color: #E3F2FD;
            font-family: 'Inter', sans-serif;
        }

        /* Styled Containers */
        .container {
            max-width: 90%;
            margin: auto;
            padding: 1.5rem;
            border-radius: 15px;
            background: rgba(255, 255, 255, 0.08);
            box-shadow: 0 8px 20px rgba(0, 123, 255, 0.1);
            backdrop-filter: blur(10px);
            transition: all 0.3s ease-in-out;
        }

        .container:hover {
            transform: scale(1.02);
        }

        /* Input Fields */
        div[data-testid="stTextInput"] input,
        div[data-testid="stTextArea"] textarea,
        div[data-testid="stSelectbox"] > div[data-baseweb="select"],
        div[data-testid="stNumberInput"] input {
            background: #001F3F;
            color: white;
            border: 2px solid #007BFF;
            border-radius: 8px;
            padding: 12px;
            font-size: 16px;
            box-shadow: 0 4px 10px rgba(0, 123, 255, 0.2);
            transition: all 0.3s ease;
        }

        /* Focus effect for inputs */
        div[data-testid="stTextInput"] input:focus,
        div[data-testid="stTextArea"] textarea:focus,
        div[data-testid="stNumberInput"] input:focus {
            border-color: #00D4FF;
            box-shadow: 0 0 10px rgba(0, 212, 255, 0.5);
            outline: none;
        }

        /* Buttons */
        .stButton>button {
            background: linear-gradient(90deg, #00D4FF, #007BFF);
            color: white;
            border-radius: 10px;
            padding: 12px;
            font-weight: bold;
            text-transform: uppercase;
            transition: all 0.3s ease-in-out;
            border: none;
            width: 100%;
            font-size: 16px;
        }

        .stButton>button:hover {
            background: linear-gradient(90deg, #007BFF, #004085);
            transform: scale(1.07);
            box-shadow: 0 8px 20px rgba(0, 123, 255, 0.3);
        }

        /* Titles & Headers */
        .title {
            color: #00D4FF;
            text-align: center;
            font-size: 32px;
            font-weight: bold;
            text-transform: uppercase;
            margin-bottom: 20px;
        }

        .subtitle {
            color: #E3F2FD;
            text-align: center;
            font-size: 20px;
            font-weight: 600;
            margin-bottom: 30px;
        }

        /* Card-Style Sections */
        .card {
            background: rgba(255, 255, 255, 0.08);
            padding: 20px;
            border-radius: 12px;
            box-shadow: 2px 2px 15px rgba(0, 123, 255, 0.2);
            margin: 20px;
            border-left: 5px solid #007BFF;
            transition: all 0.3s ease-in-out;
        }

        .card:hover {
            transform: translateY(-5px);
            box-shadow: 0 12px 25px rgba(0, 123, 255, 0.3);
        }

        /* Tables */
        .table-container {
            overflow-x: auto;
            margin-top: 20px;
        }

        .results-table {
            border-collapse: collapse;
            width: 100%;
            border-radius: 10px;
            box-shadow: 0 4px 15px rgba(0, 123, 255, 0.3);
            background: rgba(255, 255, 255, 0.1);
        }

        .results-table th {
            background: #007BFF;
            color: white;
            padding: 14px;
            text-align: left;
            font-weight: 600;
            font-size: 16px;
        }

        .results-table td {
            padding: 12px;
            border-bottom: 1px solid rgba(255, 255, 255, 0.1);
            color: white;
        }

        .results-table tr:hover td {
            background-color: rgba(0, 212, 255, 0.1);
        }

        /* Loading Spinner */
        .stSpinner {
            border: 4px solid rgba(0, 123, 255, 0.1);
            border-top: 4px solid #00D4FF;
            border-radius: 50%;
            width: 40px;
            height: 40px;
            animation: spin 1s linear infinite;
        }

        @keyframes spin {
            0% { transform: rotate(0deg); }
            100% { transform: rotate(360deg); }
        }

        /* Alerts */
        .stAlert {
            border-radius: 10px;
            padding: 1rem;
            margin: 1rem 0;
            box-shadow: 0 4px 12px rgba(0, 123, 255, 0.1);
        }

        /* Molecule Visualization */
        .molecule-card {
            background: rgba(255, 255, 255, 0.08);
            padding: 20px;
            border-radius: 12px;
            box-shadow: 2px 2px 15px rgba(0, 123, 255, 0.2);
            text-align: center;
            transition: all 0.3s ease-in-out;
        }

        .molecule-card:hover {
            transform: translateY(-5px);
            box-shadow: 0 12px 25px rgba(0, 123, 255, 0.3);
        }

        .molecule-title {
            font-weight: 600;
            margin-bottom: 1rem;
            color: #00D4FF;
            font-size: 18px;
        }

        /* Checkbox Styling */
        .stCheckbox > label > div {
            background: rgba(255, 255, 255, 0.1);
            border: 2px solid #007BFF;
            border-radius: 8px;
            padding: 0.75rem;
            transition: all 0.2s ease;
        }

        .stCheckbox > label > div:hover {
            border-color: #00D4FF;
            background-color: rgba(0, 212, 255, 0.1);
        }
    </style>
""", unsafe_allow_html=True)

# Page Navigation
if "page" not in st.session_state:
    st.session_state.page = "form"
if "api_response" not in st.session_state:
    st.session_state.api_response = None

options = dict()

def get_cid_from_name(drug_name):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{drug_name}/cids/JSON"
    response = requests.get(url)
    
    if response.status_code == 200:
        try:
            cids = response.json()["IdentifierList"]["CID"]
            return cids[0]  # Return the first matching CID
        except (KeyError, IndexError):
            return None
    else:
        return None
    
def get_pubchem_product_code(product_name):
    product_code_from_pubchem = ""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{product_name}/property/CanonicalSMILES/JSON"
    response = requests.get(url)
    if response.status_code == 200:
        try:
            smiles = response.json()["PropertyTable"]["Properties"][0]["CanonicalSMILES"]
            product_code_from_pubchem=smiles
        except (KeyError, IndexError):
            product_code_from_pubchem = "NO DRUG FOUND"
    else:
        product_code_from_pubchem="NO DRUG FOUND"
    if product_code_from_pubchem=="NO DRUG FOUND":
        return ""
    else:
        return product_code_from_pubchem

def showStructure(product_name):
    product_code = ""
    product_code_from_pubchem = get_pubchem_product_code(product_name)
    if product_code_from_pubchem=="":
        product_code_prompt = prompts.STRUCTURE_PROMPT.substitute(product_name=product_name)
        print("Prompt is: "+product_code_prompt)
        product_code = chat_with_gpt.chatWithGpt(product_code_prompt)
        if product_code == "NO DRUG FOUND":
            return ""
    else:
        product_code = product_code_from_pubchem
    
    print("product code is: "+product_code)
    print("product code from pubchem: "+product_code_from_pubchem)
    m = Chem.MolFromSmiles(product_code)
    return fig
    
# 📌 FORM PAGE
if st.session_state.page == "form":

    st.markdown('<div class="title">🧪 QAI Model - AI-Powered Quality Assurance</div>', unsafe_allow_html=True)
    st.markdown('<div class="subtitle">🔍 Enter details below to generate a pharmaceutical quality report.</div>', unsafe_allow_html=True)

    # User Input Form
    # with st.form("input_form"):
    col1, col2 = st.columns(2)

    with col1:
        options["product_name"] = st.text_input("💊 Product Name", placeholder="e.g., Paracetamol")
        if st.button("Get structure"):
            if ("product_name" not in options) or ("product_name" in options and options["product_name"]==""):
                st.error("⚠️ Please write product name!")
            else:
                with st.spinner("🛠️ Processing... Please wait"):
                    fig = showStructure(options["product_name"])
                if fig=="":
                    st.error("⚠️ Drug not found, please input a valid drug name")
                else:
                    st.image(fig, caption=f"{options["product_name"]} Molecule")
                

    with col2:
        options["quanOfMed"] = st.text_input("📦 Quantity of Medicine", placeholder="e.g., 1000 tablets")
        options["jurisdiction"] = st.selectbox("🌎 Select Jurisdiction", 
            ["INDIAN PHARMACOPIEA", "BRITISH PHARMACOPIEA", "UNITED STATES PHARMACOPOEIA", "MARTINDALE-EXTRA PHARMACOPIEA", "COMPARE WITH ALL"])
        options["powerOfDrug"] = st.text_input("⚡ Power of Drug", placeholder="e.g., 500 mg")

    # st.markdown('<div class="card">', unsafe_allow_html=True)
    options["typeOfInfo"] = st.selectbox("📊 Select Information Required:", 
            ["METHOD OF PREPARATION", "CHARACTARIZATION/EVALUATION", "Both of above", "CHECK RESULTS"])

    if options["typeOfInfo"] == "CHECK RESULTS":
        options["resultsToCheck"] = st.text_area("🔍 Enter Your Results:", height=200, placeholder="Paste lab results here...",key="checkResults")

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

            st.session_state.update(options)
            st.session_state.page = "result"
            st.experimental_rerun()

# 📌 RESULT PAGE
elif st.session_state.page == "result":

    if st.button("🔙 Go Back to Form"):
        st.session_state.page = "form"
        st.experimental_rerun()
    # Apply White Background for Result Page
    st.markdown("""
        <style>
            body { background-color: black !important; color: white !important; }
        </style>
    """, unsafe_allow_html=True)

    st.markdown('<div style="text-align:center; color:#007BFF; font-size:30px; font-weight:bold;">📑 Submission Summary</div>', unsafe_allow_html=True)

    st.markdown(f"**💊 Product Name:** {st.session_state.product_name}")
    st.markdown(f"**📦 Quantity of Medicine:** {st.session_state.quanOfMed}")
    st.markdown(f"**⚡ Power of Drug:** {st.session_state.powerOfDrug}")

    st.markdown("### 📋 Generated Report")
    if st.session_state.api_response:
        st.markdown(st.session_state.api_response)
        # components.html(st.session_state.api_response, height=1000, width=1000, scrolling=True)
    else:
        st.warning("⚠️ No response received from API.")

    if st.session_state.ftir_required:
        with st.spinner("📡 Fetching FTIR Data..."):
            ftir_data = chat_with_gpt.get_ftir_from_gpt(st.session_state.product_name)
            st.markdown("### 🔬 FTIR Data")
            st.write(ftir_data)

