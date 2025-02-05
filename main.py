import streamlit as st
import streamlit.components.v1 as components
import openai
from string import Template
import prompts
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

# Page logic
if st.session_state.page == "form":

    # Form page
    st.title("QAI Model")
    st.write("Please fill out the form below and submit.")

    # Create columns for horizontal layout
    col1, col2 = st.columns([3, 1])  # Adjust the ratio as needed

    with col1:
        # Input fields
        options["product_name"] = st.text_input("Product Name", placeholder="For example: Paracetamol")
        options["quanOfMed"] = st.text_input("Quantity of medicine", placeholder=" For example: 1000 capsules, 1000 ml")
        options["powerOfDrug"] = st.text_input("Power of drug", placeholder="For example: 10 mg")

        # Options for selection
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
            options["resultsToCheck"] = st.text_area("Write your results", placeholder="""For example: 
            The tablet has an acceptable appearance with good shape and color.
            The IR spectrum matches the expected profile for Azithromycin.""", height=250)

    with col2:
        # Button to get structure
        if st.button("Get Structure"):
            # Placeholder for future functionality
            st.write("Functionality to retrieve structure will be added here.")

    # Submit button
    if st.button("Submit"):
        if options["product_name"].strip() == "" or options["quanOfMed"].strip() == "" or options["powerOfDrug"].strip() == "":
            st.error("Please fill in all the required fields!")
        else:
            # Prepare data for ChatGPT API request
            prompt = prompts.getPromptForOptions(options)

            # Interact with ChatGPT API
            with st.spinner("Generating report..."):
                api_response = chat_with_gpt.chatWithGpt(prompt)
                st.session_state.api_response = api_response

            # Save inputs in session state
            st.session_state.product_name = options["product_name"]
            st.session_state.quanOfMed = options["quanOfMed"]
            st.session_state.typeOfInfo = options["typeOfInfo"]
            st.session_state.jurisdiction = options["jurisdiction"]
            st.session_state.powerOfDrug = options["powerOfDrug"]

            # Navigate to result page
            st.write("Click submit again to see results")
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
    st.write(f"**Product Name**: {st.session_state.product_name}")
    st.write(f"**Quantity of meds**: {st.session_state.quanOfMed}")
    st.write(f"**Power of drug**: {st.session_state.powerOfDrug}")

    # Display API response
    st.write("### Report")
    if st.session_state.api_response:
        components.html(st.session_state.api_response, height=1000, width=1000, scrolling=True)
    else:
        st.warning("No response from ChatGPT API.")
