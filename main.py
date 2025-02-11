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
        @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap');

        /* Global Styles */
        body {
            font-family: 'Inter', sans-serif;
            color: #1a1a1a;
            background-color: #f8fafc;
        }

        /* Header Styles */
        .title {
            background: linear-gradient(135deg, #2563eb, #1e40af);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            font-size: 2.5rem;
            font-weight: 700;
            text-align: center;
            margin: 2rem 0;
            padding: 1rem;
        }

        .subtitle {
            color: #4b5563;
            text-align: center;
            font-size: 1.2rem;
            font-weight: 400;
            margin-bottom: 2.5rem;
        }

        /* Card Styles */
        .card {
            background: white;
            border-radius: 12px;
            padding: 1.5rem;
            box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.1), 0 2px 4px -1px rgba(0, 0, 0, 0.06);
            margin: 1rem 0;
            border: 1px solid #e5e7eb;
            transition: transform 0.2s ease, box-shadow 0.2s ease;
        }

        .card:hover {
            transform: translateY(-2px);
            box-shadow: 0 10px 15px -3px rgba(0, 0, 0, 0.1), 0 4px 6px -2px rgba(0, 0, 0, 0.05);
        }

        /* Input Field Styles */
        div[data-testid="stTextInput"] input,
        div[data-testid="stTextArea"] textarea,
        div[data-testid="stSelectbox"] > div[data-baseweb="select"] {
            background-color: #ffffff !important;
            border: 2px solid #e2e8f0 !important;
            border-radius: 8px !important;
            padding: 0.75rem !important;
            font-size: 1rem !important;
            transition: all 0.2s ease !important;
            color: #1a1a1a !important;
        }

        div[data-testid="stTextInput"] input:focus,
        div[data-testid="stTextArea"] textarea:focus,
        div[data-testid="stSelectbox"] > div[data-baseweb="select"]:focus {
            border-color: #2563eb !important;
            box-shadow: 0 0 0 3px rgba(37, 99, 235, 0.1) !important;
            outline: none !important;
        }

        /* Button Styles */
        .stButton > button {
            background: linear-gradient(135deg, #2563eb, #1e40af) !important;
            color: white !important;
            border: none !important;
            padding: 0.75rem 1.5rem !important;
            border-radius: 8px !important;
            font-weight: 600 !important;
            font-size: 1rem !important;
            transition: all 0.2s ease !important;
            width: 100%;
            text-transform: none !important;
            letter-spacing: 0.5px !important;
        }

        .stButton > button:hover {
            transform: translateY(-2px) !important;
            box-shadow: 0 4px 12px rgba(37, 99, 235, 0.2) !important;
        }

        /* Section Headers */
        h1, h2, h3 {
            color: #1e40af;
            margin: 1.5rem 0 1rem 0;
        }

        /* Results Page Styling */
        .results-content {
            background: white;
            padding: 2rem;
            border-radius: 12px;
            box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
            margin-top: 2rem;
        }

        /* Loading Spinner */
        .stSpinner > div {
            border-top-color: #2563eb !important;
        }

        /* Alert and Info Messages */
        .stAlert {
            border-radius: 8px !important;
            padding: 1rem !important;
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
    fig = Draw.MolToImage(m, size=size)
    return fig

# Form Page
if st.session_state.page == "form":
    st.markdown('<div class="title">🧪 QAI Model - AI-Powered Quality Assurance</div>', unsafe_allow_html=True)
    st.markdown('<div class="subtitle">Advanced pharmaceutical quality analysis powered by artificial intelligence</div>', unsafe_allow_html=True)

    # User Input Form
    col1, col2 = st.columns(2)

    with col1:
        with st.container():
            st.markdown('<div class="card">', unsafe_allow_html=True)
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
                        st.image(fig, caption=f"{options['product_name']} Molecule")
            st.markdown('</div>', unsafe_allow_html=True)

    with col2:
        with st.container():
            st.markdown('<div class="card">', unsafe_allow_html=True)
            options["quanOfMed"] = st.text_input("📦 Quantity of Medicine", placeholder="e.g., 1000 tablets")
            options["jurisdiction"] = st.selectbox("🌎 Select Jurisdiction", 
                ["INDIAN PHARMACOPIEA", "BRITISH PHARMACOPIEA", "UNITED STATES PHARMACOPOEIA", "MARTINDALE-EXTRA PHARMACOPIEA", "COMPARE WITH ALL"])
            options["powerOfDrug"] = st.text_input("⚡ Power of Drug", placeholder="e.g., 500 mg")
            st.markdown('</div>', unsafe_allow_html=True)

    st.markdown('<div class="card">', unsafe_allow_html=True)
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
    st.markdown('</div>', unsafe_allow_html=True)

# Result Page
elif st.session_state.page == "result":
    if st.button("🔙 Go Back to Form"):
        st.session_state.page = "form"
        st.experimental_rerun()

    st.markdown('<div class="results-content">', unsafe_allow_html=True)
    st.markdown('<h1 style="text-align: center; color: #1e40af;">📑 Submission Summary</h1>', unsafe_allow_html=True)

    st.markdown(f"**💊 Product Name:** {st.session_state.product_name}")
    st.markdown(f"**📦 Quantity of Medicine:** {st.session_state.quanOfMed}")
    st.markdown(f"**⚡ Power of Drug:** {st.session_state.powerOfDrug}")

    st.markdown('<h2 style="color: #1e40af;">📋 Generated Report</h2>', unsafe_allow_html=True)
    if st.session_state.api_response:
        st.markdown(st.session_state.api_response)
    else:
        st.warning("⚠️ No response received from API.")

    if st.session_state.ftir_required:
        with st.spinner("📡 Fetching FTIR Data..."):
            ftir_data = chat_with_gpt.get_ftir_from_gpt(st.session_state.product_name)
            st.markdown("### 🔬 FTIR Data")
            st.write(ftir_data)
    st.markdown('</div>', unsafe_allow_html=True)