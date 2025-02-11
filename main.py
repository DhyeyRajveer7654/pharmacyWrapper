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

import streamlit as st

import streamlit as st

import streamlit as st
import time

# Custom CSS for UI Enhancements
st.markdown(
    """
    <style>
        /* Boxed layout with blue border */
        .boxed-container {
            border: 3px solid #007BFF;
            padding: 20px;
            border-radius: 10px;
            box-shadow: 0px 0px 15px rgba(0, 123, 255, 0.2);
            animation: fadeIn 0.5s ease-in-out;
        }

        /* Smooth Opening Animation */
        @keyframes fadeIn {
            from {
                opacity: 0;
                transform: scale(0.9);
            }
            to {
                opacity: 1;
                transform: scale(1);
            }
        }

        /* Input fields with blue border */
        .stTextInput, .stSelectbox {
            border: 2px solid #007BFF !important;
            border-radius: 5px;
        }

        /* Search button beside input fields */
        .search-container {
            display: flex;
            align-items: center;
            gap: 10px;
        }

        /* Style for the search button */
        .search-btn {
            background-color: #007BFF;
            color: white;
            border: none;
            padding: 5px 12px;
            border-radius: 5px;
            cursor: pointer;
            transition: 0.3s;
        }

        .search-btn:hover {
            background-color: #0056b3;
        }
    </style>
    """,
    unsafe_allow_html=True,
)

# Model UI inside the blue-bordered box
with st.container():
    st.markdown('<div class="boxed-container">', unsafe_allow_html=True)

    st.markdown("### 🔬 QAI Model - AI-Powered Quality Assurance")
    st.write("Advanced pharmaceutical quality analysis powered by artificial intelligence")

    # Product Name Input
    product_name = st.text_input("📝 Product Name", placeholder="e.g., Paracetamol")

    # Indian Pharmacopoeia with Search Button
    col1, col2 = st.columns([3, 1])
    with col1:
        jurisdiction = st.text_input("🌍 Select Jurisdiction", placeholder="Enter Jurisdiction")
    with col2:
        st.markdown('<button class="search-btn">🔍</button>', unsafe_allow_html=True)

    # Method of Preparation with Search Button
    col3, col4 = st.columns([3, 1])
    with col3:
        method = st.text_input("📖 Method of Preparation", placeholder="Enter Method")
    with col4:
        st.markdown('<button class="search-btn">🔍</button>', unsafe_allow_html=True)

    # Closing the container
    st.markdown('</div>', unsafe_allow_html=True)

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