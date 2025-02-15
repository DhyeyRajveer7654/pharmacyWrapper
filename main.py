import streamlit as st
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
st.set_page_config(page_title="QAI Model", layout="wide", page_icon="🧪")

# Apply Custom Styles
st.markdown("""
    <style>
        body { background-color: #0e1117; color: white; font-family: 'Arial', sans-serif; }
        .stTextInput>div>div>input, .stSelectbox>div>div>select, .stTextArea>div>textarea { 
            border-radius: 5px !important; padding: 10px;
        }
        div[data-testid="stTextInput"] input,
        div[data-testid="stTextArea"] textarea,
        div[data-testid="stNumberInput"] input,
        div[data-testid="stSelectbox"] > div[data-baseweb="select"],
        div[data-testid="stSlider"] > div[data-baseweb="slider"] {
            background-color: white !important;
            color: black !important;
            border: 1px solid #ccc !important;
            border-radius: 10px !important;
            box-shadow: none !important; /* Remove focus glow */
        }

        /* Add hover effect for better UX */
        div[data-testid="stTextInput"] input:hover,
        div[data-testid="stTextArea"] textarea:hover,
        div[data-testid="stNumberInput"] input:hover,
        div[data-testid="stSelectbox"] > div[role="combobox"]:hover {
            border: 1px solid #888 !important; /* Darker border on hover */
        }

        /* Ensure focus border stands out */
        div[data-testid="stTextInput"] input:focus,
        div[data-testid="stTextArea"] textarea:focus,
        div[data-testid="stNumberInput"] input:focus,
        div[data-testid="stSelectbox"] > div[role="combobox"]:focus {
            border: 1px solid #007BFF !important; /* Blue border on focus */
            outline: none !important;
        }
        .stButton>button:hover { 
            background: linear-gradient(90deg, #00D4FF, #007BFF);
            transform: scale(1.05);
        }
        .title { color: #00D4FF; text-align: center; font-size: 30px; font-weight: bold; }
        .subtitle { color: #cccccc; text-align: center; font-size: 18px; }
        .card {
            background-color: #1e222a; padding: 20px; border-radius: 12px; 
            box-shadow: 0px 4px 10px rgba(255, 255, 255, 0.1); margin: 20px;
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

# Directory where FTIR images are stored (Make sure to update this if images are in a different folder)
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
        
        if st.button("📊 Show FTIR Graph"):
            if "product_name" in options and options["product_name"]:  # Ensure product name exists
                ftir_image = get_ftir_image(options["product_name"])  # Fetch correct image

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

## Result Page
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
