import streamlit as st
import requests
import prompts
import chat_with_gpt

# Function to fetch drug structure image from PubChem API
def get_drug_structure(drug_name):
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{drug_name}/PNG"
        response = requests.get(url)
        if response.status_code == 200:
            return response.url
        else:
            return None
    except Exception as e:
        return None

# Set Page Configuration
st.set_page_config(page_title="QAI Model", layout="wide", page_icon="🧪")

# UI Styling
st.markdown("""
    <style>
        body { background-color: #0e1117; color: white; font-family: 'Arial', sans-serif; }
        .stTextInput>div>div>input { 
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
    </style>
""", unsafe_allow_html=True)

# Page Navigation
if "page" not in st.session_state:
    st.session_state.page = "form"
if "api_response" not in st.session_state:
    st.session_state.api_response = None

options = dict()

# 📌 FORM PAGE
if st.session_state.page == "form":
    st.markdown('<h2 style="text-align:center; color:#00D4FF;">🧪 QAI Model - AI-Powered Quality Assurance</h2>', unsafe_allow_html=True)

    col1, col2 = st.columns([3, 1])  # Creates two columns: Product Name (wide) & Button (small)

    with col1:
        options["product_name"] = st.text_input("📌 Product Name", placeholder="e.g., Paracetamol")

    with col2:
        if st.button("🔍 Get Structure"):
            structure_url = get_drug_structure(options["product_name"])
            if structure_url:
                st.image(structure_url, caption=f"Structure of {options['product_name']}", use_column_width=True)
            else:
                st.error("⚠️ Structure not found!")

    options["typeOfInfo"] = st.radio("📊 Select Information Required:", 
            ["METHOD OF PREPARATION", "CHARACTARIZATION/EVALUATION", "Both of above", "CHECK RESULTS", "FTIR ANALYSIS"])

    if options["typeOfInfo"] == "FTIR ANALYSIS":
        st.markdown("### 🔬 Fetching FTIR Data...")
        ftir_prompt = prompts.getPromptForOptions(options)
        with st.spinner("📡 Generating FTIR Report..."):
            ftir_data = chat_with_gpt.chatWithGpt(ftir_prompt)
            st.markdown(ftir_data, unsafe_allow_html=True)

    if st.button("🚀 Submit & Generate Report"):
        st.session_state.update(options)
        st.session_state.page = "result"
        st.experimental_rerun()
