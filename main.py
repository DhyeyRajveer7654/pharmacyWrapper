import streamlit as st
import streamlit.components.v1 as components
import openai
import prompts
import chat_with_gpt
import webbrowser

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

    # Input fields for full query
    options["product_name"] = st.text_input("Product Name", placeholder="For example: Paracetamol")
    options["quanOfMed"] = st.text_input("Quantity of medicine", placeholder="For example: 1000 capsules, 1000 ml")
    options["powerOfDrug"] = st.text_input("Power of drug", placeholder="For example: 10 mg")

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
        options["resultsToCheck"] = st.text_area("Write your results", placeholder="Enter evaluation results here...", key="checkResults", height=250)
    
    # Submit button
    if st.button("Submit"):
        if not all([options["product_name"].strip(), options["quanOfMed"].strip(), options["powerOfDrug"].strip()]):
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
            st.session_state.powerOfDrug = options["powerOfDrug"]
            st.session_state.typeOfInfo = options["typeOfInfo"]
            st.session_state.jurisdiction = options["jurisdiction"]

            # Navigate to result page
            st.write("Click submit again to see results")
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

            # Interact with ChatGPT API
            with st.spinner("Fetching chemical structure..."):
                chemical_response = chat_with_gpt.chatWithGpt(chemical_prompt)

            # Open the response in a new tab
            new_tab_url = f"https://www.google.com/search?q={chemical_product_name}+chemical+structure"
            webbrowser.open_new_tab(new_tab_url)

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
