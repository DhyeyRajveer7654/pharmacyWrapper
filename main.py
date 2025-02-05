import streamlit as st
import streamlit.components.v1 as components
import openai
import prompts
from openai import OpenAI
import chat_with_gpt
import requests

# Set page configuration
st.set_page_config(page_title="QAI Model", layout="centered")

# Initialize session state for navigation and API response
if "page" not in st.session_state:
    st.session_state.page = "form"
if "api_response" not in st.session_state:
    st.session_state.api_response = None
if "structure_image" not in st.session_state:
    st.session_state.structure_image = None

# Define functions to handle navigation
def show_form():
    st.session_state.page = "form"

def show_result():
    st.session_state.page = "result"

# Fetch chemical structure image
def fetch_structure_image(product_name):
    prompt = f"Find and display a trusted chemical structure image of {product_name} from sources like PubChem, PubMed, or Google Scholar."
    image_url = chat_with_gpt.chatWithGpt(prompt)
    return image_url

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
                st.session_state.structure_image = fetch_structure_image(product_name)
            else:
                st.error("Please enter a product name first.")

    if st.session_state.structure_image:
        st.image(st.session_state.structure_image, caption=f"Chemical Structure of {product_name}", use_column_width=True)

    options = dict()
    options["product_name"] = product_name
    options["quanOfMed"] = st.text_input("Quantity of medicine", placeholder="For example: 1000 capsules, 1000 ml")
    options["powerOfDrug"] = st.text_input("Power of drug", placeholder="For example: 10 mg")

    options["typeOfInfo"] = st.selectbox("Select information required", 
                              ["METHOD OF PREPARATION", 
                               "CHARACTARIZATION/EVALUATION", 
                               "Both of above",
                               "CHECK RESULTS"])

    options["jurisdiction"] = st.selectbox("Select jurisdiction", 
                              ["INDIAN PHARMACOPIEA", 
                               "BRITISH PHARMACOPIEA", 
                               "UNITED STATES PHARMACOPOEIA", 
                               "COMPARE WITH ALL OF THEM"])

    if options["typeOfInfo"] == "CHECK RESULTS":
        options["resultsToCheck"] = st.text_area("Write your results", placeholder="Enter evaluation results here...", height=250)

    if st.button("Submit"):
        if product_name.strip() == "" or options["quanOfMed"].strip() == "" or options["powerOfDrug"].strip() == "":
            st.error("Please fill in all the required fields!")
        else:
            prompt = prompts.getPromptForOptions(options)
            with st.spinner("Generating report..."):
                api_response = chat_with_gpt.chatWithGpt(prompt)
                st.session_state.api_response = api_response

            st.session_state.product_name = product_name
            st.session_state.quanOfMed = options["quanOfMed"]
            st.session_state.typeOfInfo = options["typeOfInfo"]
            st.session_state.jurisdiction = options["jurisdiction"]
            st.session_state.powerOfDrug = options["powerOfDrug"]

            show_result()

elif st.session_state.page == "result":
    if st.button("Go Back"):
        st.session_state.clear()
        st.experimental_rerun()

    st.title("Submission Summary")
    st.write("Thank you for your submission! Here are the details:")

    st.write(f"**Product Name**: {st.session_state.product_name}")
    st.write(f"**Quantity of meds**: {st.session_state.quanOfMed}")
    st.write(f"**Power of drug**: {st.session_state.powerOfDrug}")

    st.write("### Report")
    if st.session_state.api_response:
        components.html(st.session_state.api_response, height=1000, width=1000, scrolling=True)
    else:
        st.warning("No response from ChatGPT API.")
