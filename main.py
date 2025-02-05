import streamlit as st
import streamlit.components.v1 as components
import openai
import prompts
from openai import OpenAI
import chat_with_gpt

# Set page configuration
st.set_page_config(page_title="QAI Model", layout="centered")

# Initialize session state for navigation and API response
if "page" not in st.session_state:
    st.session_state.page = "form"  # Default page is the form page
if "api_response" not in st.session_state:
    st.session_state.api_response = None

# Options
options = dict()

# Define functions to handle navigation
def show_form():
    st.session_state.page = "form"

def show_result():
    st.session_state.page = "result"

def display_form():
    st.title("QAI Model")
    st.write("Please fill out the form below and submit.")

    # Input fields
    options["product_name"] = st.text_input("Product Name", placeholder="For example: Paracetamol")
    options["quanOfMed"] = st.text_input("Quantity of medicine", placeholder=" For example: 1000 capsules, 1000 ml")
    options["powerOfDrug"] = st.text_input("Power of drug", placeholder=" For example: 10 mg")

    # Options
    options["typeOfInfo"] = st.selectbox("Select information required", 
                              ["METHOD OF PREPARATION", 
                               "CHARACTARIZATION/EVALUATION", 
                               "Both of above",
                               "CHECK RESULTS" 
                               ])
    options["jurisdiction"] = st.selectbox("Select jurisdiction", 
                              ["INDIAN PHARMACOPIEA", 
                               "BRITISH PHARMACOPIEA", 
                               "UNITED STATES PHARMACOPOEIA", 
                               "COMPARE WITH ALL OF THEM"])

    if options["typeOfInfo"] == "CHECK RESULTS":
        options["resultsToCheck"] = st.text_area("Write your results", placeholder="""For e.g. 
        The tablet has an acceptable appearance with good shape and color.
        The IR spectrum matches the expected profile for Azithromycin.
        The HPLC results are consistent with the standard.
        The weight variation is ±2.8%.
        The tablet hardness is 5 kg.
        The friability is 0.8325%.
        The disintegration time is 23 minutes.
        The dissolution rate is 96.5%.
        The assay of Azithromycin content is 100%.
        """, key="check