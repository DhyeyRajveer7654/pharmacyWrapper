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
        body { 
            background-color: #E3F2FD; /* Light sky blue */
            color: #0B3D91; /* Deep pharma blue */
            font-family: 'Arial', sans-serif;
        }
        
        /* Styled input fields */
        .stTextInput>div>div>input, 
        .stSelectbox>div>div>select, 
        .stTextArea>div>textarea { 
            border-radius: 8px !important; 
            padding: 12px;
        }

        /* Input field styling */
        div[data-testid="stTextInput"] input,
        div[data-testid="stTextArea"] textarea,
        div[data-testid="stNumberInput"] input,
        div[data-testid="stSelectbox"] > div[data-baseweb="select"],
        div[data-testid="stSlider"] > div[data-baseweb="slider"] {
            background-color: #FFFFFF !important;
            color: #0B3D91 !important;
            border: 2px solid #007BFF !important;
            border-radius: 12px !important;
            box-shadow: 2px 2px 8px rgba(0, 123, 255, 0.2) !important;
        }

        /* Hover effect for inputs */
        div[data-testid="stTextInput"] input:hover,
        div[data-testid="stTextArea"] textarea:hover,
        div[data-testid="stNumberInput"] input:hover,
        div[data-testid="stSelectbox"] > div[role="combobox"]:hover {
            border: 2px solid #0056b3 !important;
        }

        /* Focus effect for inputs */
        div[data-testid="stTextInput"] input:focus,
        div[data-testid="stTextArea"] textarea:focus,
        div[data-testid="stNumberInput"] input:focus,
        div[data-testid="stSelectbox"] > div[role="combobox"]:focus {
            border: 2px solid #004085 !important;
            outline: none !important;
        }

        /* Button Styling */
        .stButton>button {
            background: linear-gradient(90deg, #00B4DB, #0083B0);
            color: white;
            border-radius: 10px;
            padding: 10px 20px;
            font-weight: bold;
            transition: all 0.3s ease-in-out;
            border: none;
        }

        .stButton>button:hover {
            background: linear-gradient(90deg, #007BFF, #004085);
            transform: scale(1.07);
        }
        
        /* Title and subtitle */
        .title {
            color: #004085;
            text-align: center;
            font-size: 32px;
            font-weight: bold;
            text-transform: uppercase;
            margin-bottom: 20px;
        }
        
        .subtitle {
            color: #0B3D91;
            text-align: center;
            font-size: 20px;
            font-weight: 600;
        }
        
        /* Card-style sections */
        .card {
            background-color: #FFFFFF;
            padding: 20px;
            border-radius: 12px;
            box-shadow: 2px 2px 15px rgba(0, 0, 0, 0.1);
            margin: 20px;
            border-left: 5px solid #007BFF;
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
    """Fetches the PubChem CID for a given drug name."""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{drug_name}/cids/JSON"
    response = requests.get(url)
    
    if response.status_code == 200:
        try:
            cids = response.json()["IdentifierList"]["CID"]
            return cids[0]  # Return the first matching CID
        except (KeyError, IndexError):
            return None
    return None

def get_pubchem_product_code(product_name):
    """Fetches the PubChem Canonical SMILES code for a given product name."""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{product_name}/property/CanonicalSMILES/JSON"
    response = requests.get(url)
    
    if response.status_code == 200:
        try:
            smiles = response.json()["PropertyTable"]["Properties"][0]["CanonicalSMILES"]
            return smiles
        except (KeyError, IndexError):
            return None
    return None

def showStructure(product_name):
    """Generates and returns a molecular structure image from PubChem data."""
    product_code = get_pubchem_product_code(product_name)
    
    if not product_code:
        product_code_prompt = prompts.STRUCTURE_PROMPT.substitute(product_name=product_name)
        print("Prompt is: " + product_code_prompt)
        product_code = chat_with_gpt.chatWithGpt(product_code_prompt)
        
        if product_code == "NO DRUG FOUND":
            return None  # Return None if no valid structure

    print("Product code is: " + product_code)

    try:
        m = Chem.MolFromSmiles(product_code)
        if m is not None:
            fig = Draw.MolToImage(m, size=(250, 250))
            return fig  # Return the generated molecular image
    except Exception as e:
        print(f"Error generating molecular structure: {e}")
    
    return None  # If structure generation fails, return None
    
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
    response_table = f"""
    <style>
        table {{
            width: 100%;
            border-collapse: collapse;
            background-color: #FFFFFF; /* White background for clarity */
            color: #0B3D91; /* Professional deep blue text */
            font-family: 'Arial', sans-serif;
            font-size: 16px;
            border-radius: 8px;
            box-shadow: 2px 2px 10px rgba(0, 0, 0, 0.1);
        }}
        th, td {{
            border: 1px solid #007BFF; /* Blue border */
            padding: 12px;
            text-align: left;
        }}
        th {{
            background-color: #007BFF; /* White header background */
            color: white;
        }}
        tr:nth-child(even) {{
            background-color: #007BFF; /* Blue for alternate rows */
        }}
        tr:hover {{
            background-color: #CFE2FF; /* Slightly darker blue hover effect */
        }}
    </style>
    {st.session_state.api_response}
    """
    st.markdown(response_table, unsafe_allow_html=True)
else:
    st.warning("⚠️ No response received from API.")

# 📌 FTIR Data Section
if st.session_state.ftir_required:
    with st.spinner("📡 Fetching FTIR Data..."):
        ftir_data = chat_with_gpt.get_ftir_from_gpt(st.session_state.product_name)
        st.markdown("### 🔬 FTIR Data")
        
        if ftir_data:
            ftir_table = f"""
            <table>
                <tr><th>Wavenumber (cm⁻¹)</th><th>Functional Group</th><th>Peak Description</th></tr>
                {ftir_data}
            </table>
            """
            st.markdown(ftir_table, unsafe_allow_html=True)
        else:
            st.warning("⚠️ No FTIR data available.")

