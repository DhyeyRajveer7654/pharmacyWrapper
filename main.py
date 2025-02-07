import streamlit as st
import streamlit.components.v1 as components
import prompts
import chat_with_gpt

# Set custom page layout
st.set_page_config(page_title="QAI Model", layout="centered", page_icon="🧪")

# Custom styles for a professional look
st.markdown("""
    <style>
        body { background-color: #0e1117; color: white; font-family: 'Arial', sans-serif; }
        .stTextInput>div>div>input, .stSelectbox>div>div>select, .stTextArea>div>textarea { 
            background-color: #1e222a !important; color: white !important; border-radius: 10px !important; }
        .stButton>button { background: linear-gradient(90deg, #007BFF, #00D4FF); color: white; 
            border-radius: 10px; font-size: 16px; padding: 10px; font-weight: bold; }
        .stButton>button:hover { background: linear-gradient(90deg, #00D4FF, #007BFF); }
        .reportview-container { background-color: #0e1117; }
        .stMarkdown h1, .stMarkdown h2, .stMarkdown h3 { color: #00D4FF; text-align: center; }
    </style>
""", unsafe_allow_html=True)

# Navigation setup
if "page" not in st.session_state:
    st.session_state.page = "form"  
if "api_response" not in st.session_state:
    st.session_state.api_response = None

options = dict()

# Define navigation functions
def show_form():
    st.session_state.page = "form"

def show_result():
    st.session_state.page = "result"

# FORM PAGE
if st.session_state.page == "form":

    st.markdown("## 🧪 Welcome to QAI Model - AI-Powered Quality Assurance")
    st.write("### Fill in the details below to generate a professional pharmaceutical quality report.")

    # User Input Form
    with st.form("input_form"):
        col1, col2 = st.columns(2)

        with col1:
            options["product_name"] = st.text_input("📌 Product Name", placeholder="e.g., Paracetamol")
            options["powerOfDrug"] = st.text_input("⚡ Power of Drug", placeholder="e.g., 500 mg")
            
        with col2:
            options["quanOfMed"] = st.text_input("📦 Quantity of Medicine", placeholder="e.g., 1000 tablets")
            options["jurisdiction"] = st.selectbox("🌎 Select Jurisdiction", 
                ["INDIAN PHARMACOPIEA", "BRITISH PHARMACOPIEA", "UNITED STATES PHARMACOPOEIA", "COMPARE WITH ALL"])

        st.divider()  # Horizontal line for better UI separation

        options["typeOfInfo"] = st.radio("📊 Select Information Required:", 
                ["METHOD OF PREPARATION", "CHARACTARIZATION/EVALUATION", "Both of above", "CHECK RESULTS"])

        if options["typeOfInfo"] == "CHECK RESULTS":
            options["resultsToCheck"] = st.text_area("🔍 Enter Your Results:", height=200, placeholder="Paste lab results here...")

        options["ftir_required"] = st.checkbox("📡 Retrieve FTIR Data")

        submit_button = st.form_submit_button("🚀 Submit & Generate Report")

    if submit_button:
        if not all([options["product_name"], options["quanOfMed"], options["powerOfDrug"]]):
            st.error("⚠️ Please fill in all required fields!")
        else:
            # Generate prompt
            prompt = prompts.getPromptForOptions(options)

            # Call GPT API
            with st.spinner("🛠️ Processing... Please wait"):
                api_response = chat_with_gpt.chatWithGpt(prompt)
                st.session_state.api_response = api_response

            # Save user input
            st.session_state.update(options)

            # Navigate to result page
            st.write("✅ Data submitted successfully! Click 'Submit' again to see results.")
            show_result()

# RESULT PAGE
elif st.session_state.page == "result":
    st.markdown("## 📑 Submission Summary")

    # Display input details
    st.markdown(f"**📌 Product Name:** {st.session_state.product_name}")
    st.markdown(f"**📦 Quantity of Medicine:** {st.session_state.quanOfMed}")
    st.markdown(f"**⚡ Power of Drug:** {st.session_state.powerOfDrug}")

    # Display API response
    st.write("### 📋 Generated Report")
    if st.session_state.api_response:
        components.html(st.session_state.api_response, height=1000, width=1000, scrolling=True)
    else:
        st.warning("⚠️ No response received from GPT API.")

    # FTIR Data Retrieval
    if st.session_state.ftir_required:
        with st.spinner("📡 Fetching FTIR Data..."):
            ftir_data = chat_with_gpt.get_ftir_from_gpt(st.session_state.product_name)
            st.write("### 🔬 FTIR Data from GPT-4")
            st.write(ftir_data)

    if st.button("🔙 Go Back to Form"):
        st.session_state.clear()  
        st.experimental_rerun()
