import streamlit as st

# Add a button with dropdown options in the upper right corner
import streamlit as st

# Create a button in the top right corner using Streamlit columns
col1, col2 = st.columns([8, 1])  # Adjust column widths to position the button

with col2:
    selected_option = st.selectbox("Options", ["Select", "Pharmacopoeial Compliance", "Regulatory Compliance"])

# Handle the selected option
if selected_option == "Pharmacopoeial Compliance":
    st.write("You selected Pharmacopoeial Compliance.")
elif selected_option == "Regulatory Compliance":
    st.write("You selected Regulatory Compliance.")
import streamlit.components.v1 as components
import prompts
import chat_with_gpt
from string import Template
from rdkit import Chem
from rdkit.Chem import Draw
import requests
import os
import streamlit as st

size = (250, 250)

# Set Page Configuration
st.set_page_config(page_title="QAI Model AI-Powered Quality Assistance", layout="wide", page_icon="🧪")

# Enhanced Professional UI Styling
st.markdown("""
    <style>
        /* Global Styles */
        body {
            background-color: #f0f2f6;
            color: #1e293b;
            font-family: 'Inter', 'sans serif';
        }

        /* Header Styling */
        .main-header {
            background: linear-gradient(135deg, #0052cc, #00a3bf);
            color: white;
            padding: 2rem;
            border-radius: 10px;
            margin-bottom: 2rem;
            text-align: center;
        }

        /* Form Elements */
        div[data-testid="stTextInput"] input,
        div[data-testid="stTextArea"] textarea,
        div[data-testid="stSelectbox"] > div[data-baseweb="select"] {
            background-color: white;
            border: 1px solid #e2e8f0;
            border-radius: 8px;
            padding: 0.75rem;
            font-size: 1rem;
            transition: all 0.3s ease;
            box-shadow: 0 2px 4px rgba(0,0,0,0.05);
        }

        div[data-testid="stTextInput"] input:focus,
        div[data-testid="stTextArea"] textarea:focus {
            border-color: #0052cc;
            box-shadow: 0 0 0 2px rgba(0,82,204,0.2);
        }

        /* Button Styling */
        .stButton > button {
            width: 100%;
            background: linear-gradient(135deg, #0052cc, #00a3bf);
            color: white;
            border: none;
            padding: 0.75rem 1.5rem;
            border-radius: 8px;
            font-weight: 600;
            transition: all 0.3s ease;
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }

        .stButton > button:hover {
            transform: translateY(-2px);
            box-shadow: 0 4px 12px rgba(0,82,204,0.2);
        }

        /* Card Styling */
        .card div[data-testid="stSelectbox"]{
            background: white;
            border-radius: 10px;
            padding: 1.5rem;
            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
            margin-bottom: 1.5rem;
        }

        /* Results Table */
        .table-container {
            background: white;
            border-radius: 10px;
            padding: 1rem;
            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
            border-collapse: collapse;
        }

        table {
            width: 100%;
            border-collapse: collapse;
            border-spacing: 0;
            margin: 1rem 0;
        }

        th {
            background: #0052cc;
            color: white;
            padding: 1rem;
            text-align: left;
            font-weight: 600;
        }

        td {
            padding: 1rem;
            border-bottom: 1px solid #e2e8f0;
        }

        tr:hover {
            background: #f8fafc;
        }

        /* Loading Spinner */
        .stSpinner > div {
            border-color: #0052cc !important;
        }

        /* Success Message */
        .success-message {
            background: #dcfce7;
            color: #166534;
            padding: 1rem;
            border-radius: 8px;
            margin: 1rem 0;
        }

        /* Error Message */
        .error-message {
            background: #fee2e2;
            color: #991b1b;
            padding: 1rem;
            border-radius: 8px;
            margin: 1rem 0;
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
    if m:
        return Draw.MolToImage(m, size=(400, 400))
    return None

# Directory where FTIR images are stored
FTIR_IMAGE_DIR = "./"

def get_ftir_image(product_name):
    """Fetches the corresponding FTIR image for the given product name."""
    image_filename = f"{product_name.lower()}.png"
    image_path = os.path.join(FTIR_IMAGE_DIR, image_filename)
    if os.path.exists(image_path):
        return image_path
    return None

# 📌 FORM PAGE
if st.session_state.page == "form":
    st.markdown('<div class="main-header"><h1>🧪 QAI Model AI-Powered Quality Assistance</h1><p> CREATED BY :- MEERA ACHARYA & RAJ PATEL</P><p>Enter details below to generate a comprehensive quality report</p></div>', unsafe_allow_html=True)

    # User Input Form in a card layout
    # st.markdown('<div class="card">', unsafe_allow_html=True)
    col1, col2 = st.columns(2)

    with col1:
        options["product_name"] = st.text_input("💊 Product Name", placeholder="e.g., Paracetamol")
        if st.button("🔬 Get Structure"):
            if not options["product_name"]:
                st.error("⚠️ Please write product name!")
            else:
                with st.spinner("🛠️ Processing... Please wait"):
                    fig = showStructure(options["product_name"])
                if fig == "":
                    st.error("⚠️ Drug not found, please input a valid drug name")
                else:
                    st.image(fig, caption=f"{options['product_name']} Molecule")

        if st.button("📊 Show FTIR Graph"):
            if options.get("product_name"):  # Ensure product name exists
                ftir_image = get_ftir_image(options["product_name"])
                if ftir_image:
                    st.image(ftir_image, caption=f"FTIR Graph for {options['product_name']}", use_column_width=True)
                else:
                    st.error(f"⚠️ No FTIR data available for {options['product_name']}.")
            else:
                st.error("⚠️ Please enter a product name.")            

    with col2:
        options["quanOfMed"] = st.text_input("📦 Quantity of Medicine", placeholder="e.g., 1000 tablets")
        options["jurisdiction"] = st.selectbox("🌎 Select Jurisdiction", 
            ["INDIAN PHARMACOPIEA", "BRITISH PHARMACOPIEA", "UNITED STATES PHARMACOPOEIA", "MARTINDALE-EXTRA PHARMACOPIEA", "COMPARE WITH ALL"])
        options["powerOfDrug"] = st.text_input("⚡ Power of Drug", placeholder="e.g., 500 mg")

    # st.markdown('</div>', unsafe_allow_html=True)

    # Analysis Options in a separate card
    # st.markdown('<div class="card">', unsafe_allow_html=True)
    options["typeOfInfo"] = st.selectbox("📊 Select Analysis Type:", 
            ["METHOD OF PREPARATION", "CHARACTARIZATION/EVALUATION", "Both of above", "CHECK RESULTS"])

    if options["typeOfInfo"] == "CHECK RESULTS":
        options["resultsToCheck"] = st.text_area("🔍 Enter Your Results:", height=200, placeholder="Paste lab results here...", key="checkResults")

    options["ftir_required"] = st.checkbox("📡 Include FTIR Analysis")
    # st.markdown('</div>', unsafe_allow_html=True)

    # Submit button with enhanced styling
    submit_button = st.button("🚀 Generate Report")
    if submit_button:
        if not all([options.get("product_name"), options.get("quanOfMed"), options.get("powerOfDrug")]):
            st.error("⚠️ Please fill in all required fields!")
        else:
            prompt = prompts.getPromptForOptions(options)
            with st.spinner("🛠️ Generating comprehensive report... Please wait"):
                api_response = chat_with_gpt.chatWithGpt(prompt)
                st.session_state.api_response = api_response

            st.session_state.update(options)
            st.session_state.page = "result"
            st.experimental_rerun()

# 📌 RESULT PAGE
elif st.session_state.page == "result":
    st.markdown('<div class="main-header"><h1>📑 Quality Analysis Report</h1></div>', unsafe_allow_html=True)
    
    if st.button("🔙 Return to Form", key="back_button"):
        st.session_state.page = "form"
        st.experimental_rerun()

    # st.markdown('<div class="card">', unsafe_allow_html=True)
    st.markdown("### 📋 Analysis Details")
    st.markdown(f"**💊 Product:** {st.session_state.product_name}",)
    st.markdown(f"**📦 Quantity:** {st.session_state.quanOfMed}")
    st.markdown(f"**⚡ Strength:** {st.session_state.powerOfDrug}")
    # st.markdown('</div>', unsafe_allow_ht ml=True)

    if st.session_state.get("ftir_required"):
        with st.spinner("📡 Analyzing FTIR Data..."):
            ftir_data = chat_with_gpt.get_ftir_from_gpt(st.session_state.product_name)
            components.html("### 🔬 FTIR Analysis")
            st.markdown(ftir_data, unsafe_allow_html=True)
            # components.html(ftir_data)

    if st.session_state.api_response:
        components.html("<div class='table-container'>"+st.session_state.api_response+"</div>",height=800,width=1000,scrolling=True)
    else:
        st.warning("⚠️ No response received. Please try again.")
