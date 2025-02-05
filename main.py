import streamlit as st
import streamlit.components.v1 as components
import openai
import prompts
from openai import OpenAI
import chat_with_gpt
from rdkit import Chem
from rdkit.Chem import Draw
from io import BytesIO
from PIL import Image

# Set page configuration
st.set_page_config(page_title="QAI Model", layout="centered")

# Initialize session state
if "page" not in st.session_state:
    st.session_state.page = "form"
if "api_response" not in st.session_state:
    st.session_state.api_response = None
if "structure_image" not in st.session_state:
    st.session_state.structure_image = None

# Navigation functions
def show_form():
    st.session_state.page = "form"

def show_result():
    st.session_state.page = "result"

# Function to generate structure using RDKit
def generate_structure():
    smiles = prompts.get_smiles_for_product(st.session_state.product_name)
    if smiles:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            fig = Draw.MolToImage(mol)
            st.session_state.structure_image = fig

# Page logic
if st.session_state.page == "form":
    st.title("QAI Model")
    st.write("Please fill out the form below and submit.")

    col1, col2 = st.columns([3, 1])
    with col1:
        product_name = st.text_input("Product Name", placeholder="For example: Paracetamol")
    with col2:
        if st.button("Get Structure"):
            if product_name.strip():
                st.session_state.product_name = product_name
                generate_structure()
            else:
                st.warning("Please enter a product name first.")

    quanOfMed = st.text_input("Quantity of medicine", placeholder="For example: 1000 capsules, 1000 ml")
    powerOfDrug = st.text_input("Power of drug", placeholder="For example: 10 mg")
    
    typeOfInfo = st.selectbox("Select information required", [
        "METHOD OF PREPARATION",
        "CHARACTERIZATION/EVALUATION",
        "Both of above",
        "CHECK RESULTS"
    ])

    jurisdiction = st.selectbox("Select jurisdiction", [
        "INDIAN PHARMACOPIEA",
        "BRITISH PHARMACOPIEA",
        "UNITED STATES PHARMACOPOEIA",
        "COMPARE WITH ALL OF THEM"
    ])

    if typeOfInfo == "CHECK RESULTS":
        resultsToCheck = st.text_area("Write your results", placeholder="Enter evaluation results here...", height=250)

    if st.button("Submit"):
        if not product_name.strip() or not quanOfMed.strip() or not powerOfDrug.strip():
            st.error("Please fill in all the required fields!")
        else:
            st.session_state.product_name = product_name
            st.session_state.quanOfMed = quanOfMed
            st.session_state.powerOfDrug = powerOfDrug
            st.session_state.typeOfInfo = typeOfInfo
            st.session_state.jurisdiction = jurisdiction

            prompt = prompts.getPromptForOptions({
                "product_name": product_name,
                "quanOfMed": quanOfMed,
                "powerOfDrug": powerOfDrug,
                "typeOfInfo": typeOfInfo,
                "jurisdiction": jurisdiction,
                "resultsToCheck": resultsToCheck if typeOfInfo == "CHECK RESULTS" else None
            })

            with st.spinner("Generating report..."):
                st.session_state.api_response = chat_with_gpt.chatWithGpt(prompt)

            st.write("Click submit again to see results")
            show_result()

elif st.session_state.page == "result":
    st.title("Submission Summary")
    st.write(f"**Product Name**: {st.session_state.product_name}")
    st.write(f"**Quantity of meds**: {st.session_state.quanOfMed}")
    st.write(f"**Power of drug**: {st.session_state.powerOfDrug}")

    st.write("### Report")
    if st.session_state.api_response:
        components.html(st.session_state.api_response, height=1000, width=1000, scrolling=True)
    else:
        st.warning("No response from ChatGPT API.")

    # Display generated structure if available
    if st.session_state.structure_image:
        st.write("### Molecular Structure")
        st.image(st.session_state.structure_image, caption="Generated Structure", use_column_width=True)

    if st.button("Go Back"):
        st.session_state.clear()
        st.experimental_rerun()
