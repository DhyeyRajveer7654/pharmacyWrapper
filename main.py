import streamlit as st
import openai
import prompts  # Ensure that prompts is imported
import chat_with_gpt  # Assuming this is the module for interacting with GPT
import base64

# Set page configuration
st.set_page_config(page_title="QAI Model", layout="centered")

# Initialize session state for navigation and API response
if "page" not in st.session_state:
    st.session_state.page = "form"  # Default page is the form page
if "api_response" not in st.session_state:
    st.session_state.api_response = None
if "chemical_response" not in st.session_state:
    st.session_state.chemical_response = None

# Options
def show_form():
    st.session_state.page = "form"

def show_result():
    st.session_state.page = "result"

# Page logic
if st.session_state.page == "form":
    # Form page
    st.title("QAI Model")
    st.write("Please fill out the form below and submit.")
    
    # Input fields for full query
    product_name = st.text_input("Product Name", placeholder="For example: Paracetamol")
    quanOfMed = st.text_input("Quantity of medicine", placeholder="For example: 1000 capsules, 1000 ml")
    powerOfDrug = st.text_input("Power of drug", placeholder="For example: 10 mg")
    
    # Options
    typeOfInfo = st.selectbox("Select information required", [
        "METHOD OF PREPARATION",
        "CHARACTARIZATION/EVALUATION",
        "Both of above",
        "CHECK RESULTS"
    ])
    
    jurisdiction = st.selectbox("Select jurisdiction", [
        "INDIAN PHARMACOPIEA",
        "BRITISH PHARMACOPIEA",
        "UNITED STATES PHARMACOPOEIA",
        "COMPARE WITH ALL OF THEM"
    ])
    
    resultsToCheck = ""
    if typeOfInfo == "CHECK RESULTS":
        resultsToCheck = st.text_area("Write your results", placeholder="Enter evaluation results here...", height=250)
    
    # Submit button
    if st.button("Submit"):
        if not all([product_name.strip(), quanOfMed.strip(), powerOfDrug.strip()]):
            st.error("Please fill in all the required fields!")
        else:
            # Prepare data for ChatGPT API request
            options = {
                "product_name": product_name,
                "quanOfMed": quanOfMed,
                "powerOfDrug": powerOfDrug,
                "typeOfInfo": typeOfInfo,
                "jurisdiction": jurisdiction,
                "resultsToCheck": resultsToCheck if typeOfInfo == "CHECK RESULTS" else ""
            }
            prompt = prompts.getPromptForOptions(options)
            
            # Interact with ChatGPT API request
            with st.spinner("Generating report..."):
                st.session_state.api_response = chat_with_gpt.chatWithGpt(prompt)
            
            # Save inputs in session state
            st.session_state.update(options)
            
            # Navigate to result page
            show_result()
    
    st.write("---")
    
    # New section for chemical structure query
    st.subheader("Get Chemical Structure")
    chemical_product_name = st.text_input("Enter only the Product Name", placeholder="For example: Aspirin")
    
    if st.button("Get Chemical Structure"):
        if chemical_product_name.strip() == "":
            st.error("Please enter a product name!")
        else:
            # Generate prompt for chemical structure
            chemical_prompt = prompts.getPromptForChemicalStructure(chemical_product_name)
            
            # Interact with ChatGPT API for chemical structure
            with st.spinner("Fetching chemical structure..."):
                st.session_state.chemical_response = chat_with_gpt.chatWithGpt(chemical_prompt)
            
            # Navigate to result page
            show_result()

elif st.session_state.page == "result":
    # Button to go back to the form page
    if st.button("Go Back"):
        st.session_state.clear()  # Clears session state
        st.experimental_rerun()
    
    # Result page
    st.title("Submission Summary")
    st.write("Thank you for your submission! Here are the details:")
    
    # Display submitted details
    st.write(f"**Product Name**: {st.session_state.get('product_name', 'N/A')}")
    st.write(f"**Quantity of meds**: {st.session_state.get('quanOfMed', 'N/A')}")
    st.write(f"**Power of drug**: {st.session_state.get('powerOfDrug', 'N/A')}")
    
    # Display API response
    st.write("### Report")
    if st.session_state.api_response:
        st.write(st.session_state.api_response)
    else:
        st.warning("No response from ChatGPT API.")
    
    # Display Chemical Structure Result (Image)
    st.write("### Chemical Structure")
    
    if st.session_state.chemical_response:
        # Assuming the AI response is base64-encoded image data
        try:
        import base64
from io import BytesIO
from PIL import Image

if st.session_state.chemical_response:
    # If the response is a base64 string, handle it
    if st.session_state.chemical_response.startswith('data:image/png;base64,'):
        img_data = st.session_state.chemical_response.split(",")[1]
        img_data = base64.b64decode(img_data)
        img = Image.open(BytesIO(img_data))
        st.image(img, caption="Chemical Structure Image", use_column_width=True)
    else:
        st.warning("The provided response is not a valid image.")
