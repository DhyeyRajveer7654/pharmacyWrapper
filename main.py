import streamlit as st
import streamlit.components.v1 as components
import prompts
import chat_with_gpt

# Set Page Configuration with custom theme
st.set_page_config(page_title="QAI Model", layout="wide", page_icon="🧪")

# Apply Enhanced Styling
st.markdown("""
    <style>
        /* Background and Main App */
        .stApp { background: linear-gradient(135deg, #1e1e2e, #3a3a55) !important; color: white; }
        body { font-family: 'Arial', sans-serif; }

        /* Form Cards */
        .card { background: rgba(255, 255, 255, 0.1); padding: 20px; border-radius: 12px; box-shadow: 0px 4px 10px rgba(255, 255, 255, 0.2); margin: 20px; }

        /* Input Fields */
        .stTextInput > div > div > input, .stSelectbox > div > div > select, .stTextArea > div > textarea { 
            background-color: #2a2d40 !important; color: white !important; border-radius: 8px !important; padding: 10px;
            border: 1px solid #555;
        }

        /* Buttons */
        .stButton>button { 
            background: linear-gradient(90deg, #ff7eb3, #ff758c); 
            color: white; border-radius: 10px; font-size: 16px; padding: 12px; font-weight: bold; border: none;
            transition: all 0.3s ease-in-out; cursor: pointer;
        }
        .stButton>button:hover { 
            background: linear-gradient(90deg, #ff758c, #ff7eb3);
            transform: scale(1.05);
        }

        /* Titles and Headers */
        .title { color: #ff758c; text-align: center; font-size: 32px; font-weight: bold; text-shadow: 2px 2px 10px rgba(255, 117, 140, 0.6); }
        .subtitle { color: #dddddd; text-align: center; font-size: 18px; margin-bottom: 20px; }

    </style>
""", unsafe_allow_html=True)

# Page Navigation
if "page" not in st.session_state:
    st.session_state.page = "form"
if "api_response" not in st.session_state:
    st.session_state.api_response = None

options = {}

# 📌 FORM PAGE
if st.session_state.page == "form":

    st.markdown('<div class="title">🧪 QAI Model - AI-Powered Quality Assurance</div>', unsafe_allow_html=True)
    st.markdown('<div class="subtitle">🔍 Enter details below to generate a pharmaceutical quality report.</div>', unsafe_allow_html=True)

    with st.form("input_form"):
        col1, col2 = st.columns(2)

        with col1:
            options["product_name"] = st.text_input("💊 Product Name", placeholder="e.g., Paracetamol")
            options["powerOfDrug"] = st.text_input("⚡ Power of Drug", placeholder="e.g., 500 mg")

        with col2:
            options["quanOfMed"] = st.text_input("📦 Quantity of Medicine", placeholder="e.g., 1000 tablets")
            options["jurisdiction"] = st.selectbox("🌎 Select Jurisdiction", 
                ["INDIAN PHARMACOPIEA", "BRITISH PHARMACOPIEA", "UNITED STATES PHARMACOPOEIA", "COMPARE WITH ALL"])

        st.markdown('<div class="card">', unsafe_allow_html=True)
        options["typeOfInfo"] = st.radio("📊 Select Information Required:", 
                ["METHOD OF PREPARATION", "CHARACTERIZATION/EVALUATION", "Both of above", "CHECK RESULTS"])

        if options["typeOfInfo"] == "CHECK RESULTS":
            options["resultsToCheck"] = st.text_area("🔍 Enter Your Results:", height=150, placeholder="Paste lab results here...")

        options["ftir_required"] = st.checkbox("📡 Retrieve FTIR Data")

        submit_button = st.form_submit_button("🚀 Generate Report")

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

# 📌 RESULT PAGE
elif st.session_state.page == "result":
    st.markdown('<div class="title">📑 Submission Summary</div>', unsafe_allow_html=True)

    st.markdown(f"**💊 Product Name:** {st.session_state.product_name}")
    st.markdown(f"**📦 Quantity of Medicine:** {st.session_state.quanOfMed}")
    st.markdown(f"**⚡ Power of Drug:** {st.session_state.powerOfDrug}")

    st.markdown("### 📋 Generated Report")
    if st.session_state.api_response:
        components.html(st.session_state.api_response, height=1000, width=1000, scrolling=True)
    else:
        st.warning("⚠️ No response received from GPT API.")

    if st.session_state.ftir_required:
        with st.spinner("📡 Fetching FTIR Data..."):
            ftir_data = chat_with_gpt.get_ftir_from_gpt(st.session_state.product_name)
            st.markdown("### 🔬 FTIR Data")
            st.write(ftir_data)

    if st.button("🔙 Go Back to Form"):
        st.session_state.page = "form"
        st.experimental_rerun()
