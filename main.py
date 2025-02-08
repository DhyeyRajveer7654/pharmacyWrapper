import streamlit as st
import streamlit.components.v1 as components
import prompts
import chat_with_gpt

# Set Page Configuration
st.set_page_config(page_title="QAI Model", layout="wide", page_icon="🧪")

# Page State Management
if "page" not in st.session_state:
    st.session_state.page = "form"
if "api_response" not in st.session_state:
    st.session_state.api_response = None

options = dict()

# 📌 FORM PAGE
if st.session_state.page == "form":
    st.title("🧪 QAI Model - AI-Powered Quality Assurance")

    with st.form("input_form"):
        options["product_name"] = st.text_input("💊 Product Name", placeholder="e.g., Paracetamol")
        options["powerOfDrug"] = st.text_input("⚡ Power of Drug", placeholder="e.g., 500 mg")
        options["quanOfMed"] = st.text_input("📦 Quantity of Medicine", placeholder="e.g., 1000 tablets")

        options["jurisdiction"] = st.selectbox("🌎 Select Jurisdiction", 
            ["INDIAN PHARMACOPIEA", "BRITISH PHARMACOPIEA", "UNITED STATES PHARMACOPOEIA", "COMPARE WITH ALL"])

        options["typeOfInfo"] = st.radio("📊 Select Information Required:", 
            ["METHOD OF PREPARATION", "CHARACTERIZATION/EVALUATION", "Both of above", "CHECK RESULTS"])

        if options["typeOfInfo"] == "CHECK RESULTS":
            options["resultsToCheck"] = st.text_area("🔍 Enter Your Results:", height=200)

        options["ftir_required"] = st.checkbox("📡 Retrieve FTIR Data")

        submit_button = st.form_submit_button("🚀 Generate Report")

    if submit_button:
        prompt = prompts.getPromptForOptions(options)
        api_response = chat_with_gpt.chatWithGpt(prompt)
        st.session_state.api_response = api_response
        st.session_state.page = "result"
        st.experimental_rerun()

# 📌 RESULT PAGE
elif st.session_state.page == "result":
    st.title("📑 Submission Summary")
    st.write(f"**💊 Product Name:** {st.session_state.product_name}")
    st.write(f"**📦 Quantity:** {st.session_state.quanOfMed}")
    st.write(f"**⚡ Power:** {st.session_state.powerOfDrug}")

    st.markdown("### 📋 Generated Report")
    if st.session_state.api_response:
        components.html(st.session_state.api_response, height=1000, width=1000, scrolling=True)
    else:
        st.warning("⚠️ No response received.")

    if st.button("🔙 Go Back"):
        st.session_state.page = "form"
        st.experimental_rerun()
