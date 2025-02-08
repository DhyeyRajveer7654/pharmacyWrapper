import streamlit as st
import streamlit.components.v1 as components
import requests
import prompts
import chat_with_gpt

# Set Page Configuration
st.set_page_config(page_title="QAI Model", layout="wide", page_icon="🧪")

# Apply Custom Styles
st.markdown("""
    <style>
        body { background-color: #0e1117; color: white; font-family: 'Arial', sans-serif; }
        .stTextInput>div>div>input, .stSelectbox>div>div>select, .stTextArea>div>textarea { 
            background-color: #1e222a !important; color: white !important; border-radius: 10px !important; padding: 10px;
        }
        .stButton>button { 
            background: linear-gradient(90deg, #007BFF, #00D4FF); 
            color: white; border-radius: 10px; font-size: 16px; padding: 10px; font-weight: bold; border: none;
            transition: 0.3s;
        }
        .stButton>button:hover { 
            background: linear-gradient(90deg, #00D4FF, #007BFF);
            transform: scale(1.05);
        }
        .title { color: #00D4FF; text-align: center; font-size: 30px; font-weight: bold; }
        .subtitle { color: #cccccc; text-align: center; font-size: 18px; }
        .card {
            background-color: #1e222a; padding: 20px; border-radius: 12px; 
            box-shadow: 0px 4px 10px rgba(255, 255, 255, 0.1); margin: 20px;
        }
    </style>
""", unsafe_allow_html=True)

# Function to fetch structure from PubChem
def fetch_pubchem_structure(drug_name):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{drug_name}/PNG"
    response = requests.get(url)
    if response.status_code == 200:
        return url
    return None

# Page Navigation
if "page" not in st.session_state:
    st.session_state.page = "form"
if "api_response" not in st.session_state:
    st.session_state.api_response = None
if "structure_image" not in st.session_state:
    st.session_state.structure_image = None

options = dict()

# 📌 FORM PAGE
if st.session_state.page == "form":

    st.markdown('<div class="title">🧪 QAI Model - AI-Powered Quality Assurance</div>', unsafe_allow_html=True)
    st.markdown('<div class="subtitle">🔍 Enter details below to generate a pharmaceutical quality report.</div>', unsafe_allow_html=True)

    # User Input Form
    with st.form("input_form"):
        col1, col2 = st.columns([3, 1])

        with col1:
            options["product_name"] = st.text_input("💊 Product Name", placeholder="e.g., Paracetamol")

        with col2:
            if st.form_submit_button("🔍 Get Structure"):
                if options["product_name"]:
                    structure_url = fetch_pubchem_structure(options["product_name"])
                    if structure_url:
                        st.session_state.structure_image = structure_url
                    else:
                        st.session_state.structure_image = None
                        st.warning("⚠️ No structure found on PubChem.")

        if st.session_state.structure_image:
            st.image(st.session_state.structure_image, caption="Chemical Structure from PubChem", use_column_width=True)

        col3, col4 = st.columns(2)

        with col3:
            options["powerOfDrug"] = st.text_input("⚡ Power of Drug", placeholder="e.g., 500 mg")

        with col4:
            options["quanOfMed"] = st.text_input("📦 Quantity of Medicine", placeholder="e.g., 1000 tablets")

        options["jurisdiction"] = st.selectbox("🌎 Select Jurisdiction", 
                ["INDIAN PHARMACOPIEA", "BRITISH PHARMACOPIEA", "UNITED STATES PHARMACOPOEIA", "MARTINDALE-EXTRA PHARMACOPIEA", "COMPARE WITH ALL"])

        st.markdown('<div class="card">', unsafe_allow_html=True)
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
            prompt = prompts.getPromptForOptions(options)
            with st.spinner("🛠️ Processing... Please wait"):
                api_response = chat_with_gpt.chatWithGpt(prompt)
                st.session_state.api_response = api_response

            st.session_state.update(options)
            st.session_state.page = "result"
            st.experimental_rerun()

# 📌 RESULT PAGE
elif st.session_state.page == "result":

    # Apply White Background for Result Page
    st.markdown("""
        <style>
            body { background-color: black !important; color: white !important; }
        </style>
    """, unsafe_allow_html=True)

    st.markdown('<div style="text-align:center; color:#007BFF; font-size:30px; font-weight:bold;">📑 Submission Summary</div>', unsafe_allow_html=True)

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
