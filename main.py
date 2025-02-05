import streamlit as st
import streamlit.components.v1 as components
import openai
import prompts
import chat_with_gpt
import requests
from PIL import Image
from io import BytesIO

# Set page configuration
st.set_page_config(page_title="QAI Model", layout="centered")

# Initialize session state
if "page" not in st.session_state:
    st.session_state.page = "form"
if "api_response" not in st.session_state:
    st.session_state.api_response = None
if "chemical_structure_image" not in st.session_state:
    st.session_state.chemical_structure_image = None

# Navigation functions
def show_form():
    st.session_state.page = "form"

def show_result():
    st.session_state.page = "result"

# Fetch Chemical Structure Image
def fetch_chemical_structure(drug_name):
    prompt = f"Provide a high-quality image (PNG, JPEG, or JPG) of the chemical structure of {drug_name}. The structure must be sourced from trusted scientific databases such as Pharmacopoeias, PubChem, PubMed, NCBI, or peer-reviewed scholarly articles. Include the IUPAC name and molecular formula alongside the structure if available. Ensure the structure matches the standard representation from these authoritative sources."
    
    response = chat_with_gpt.chatWithGpt(prompt)
    
    # Simulate image fetching by interpreting the API response as a URL (adjust based on actual API behavior)
    if response and response.startswith("http"):
        img_response = requests.get(response)
        if img_response.status_code == 200:
            st.session_state.chemical_structure_image = Image.open(BytesIO(img_response.content))
        else:
            st.error("Failed to fetch the chemical structure image.")
    else:
        st.warning("No image URL found in API response.")

# Page Logic
if st.session_state.page == "form":
    st.title("QAI Model")
    st.write("Please fill out the form below and submit.")

    # Input Fields
    product_name = st.text_input("Product Name", placeholder="For example: Paracetamol")
    st.session_state.product_name = product_name  # Save for later use
    
    # Button to fetch chemical structure
    if st.button("Get Chemical Structure"):
        if product_name.strip():
            with st.spinner("Fetching chemical structure image..."):
                fetch_chemical_structure(product_name)
        else:
            st.error("Please enter a product name to fetch the chemical structure.")

    # Display fetched image
    if st.session_state.chemical_structure_image:
        st.image(st.session_state.chemical_structure_image, caption=f"Chemical Structure of {product_name}", use_column_width=True)

    # Other Input Fields
    quan_of_med = st.text_input("Quantity of medicine", placeholder="For example: 1000 capsules, 1000 ml")
    power_of_drug = st.text_input("Power of drug", placeholder="For example: 10 mg")
    
    type_of_info = st.selectbox("Select information required", ["METHOD OF PREPARATION", "CHARACTARIZATION/EVALUATION", "Both of above", "CHECK RESULTS"])
    jurisdiction = st.selectbox("Select jurisdiction", ["INDIAN PHARMACOPIEA", "BRITISH PHARMACOPIEA", "UNITED STATES PHARMACOPOEIA", "COMPARE WITH ALL OF THEM"])

    if type_of_info == "CHECK RESULTS":
        results_to_check = st.text_area("Write your results", placeholder="Enter detailed evaluation results here", height=250)

    # Submit Button
    if st.button("Submit"):
        if not product_name.strip() or not quan_of_med.strip() or not power_of_drug.strip():
            st.error("Please fill in all the required fields!")
        else:
            options = {
                "product_name": product_name,
                "quanOfMed": quan_of_med,
                "powerOfDrug": power_of_drug,
                "typeOfInfo": type_of_info,
                "jurisdiction": jurisdiction
            }
            if type_of_info == "CHECK RESULTS":
                options["resultsToCheck"] = results_to_check

            # Generate report
            with st.spinner("Generating report..."):
                api_response = chat_with_gpt.chatWithGpt(prompts.getPromptForOptions(options))
                st.session_state.api_response = api_response

            # Navigate to result page
            st.write("Click submit again to see results")
            show_result()

elif st.session_state.page == "result":
    if st.button("Go Back"):
        st.session_state.clear()
        st.experimental_rerun()

    st.title("Submission Summary")
    st.write(f"**Product Name**: {st.session_state.product_name}")

    if st.session_state.chemical_structure_image:
        st.image(st.session_state.chemical_structure_image, caption=f"Chemical Structure of {st.session_state.product_name}", use_column_width=True)

    st.write("### Report")
    if st.session_state.api_response:
        components.html(st.session_state.api_response, height=1000, width=1000, scrolling=True)
    else:
        st.warning("No response from ChatGPT API.")
