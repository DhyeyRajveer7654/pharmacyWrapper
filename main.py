import streamlit as st
import streamlit.components.v1 as components
import prompts
import chat_with_gpt
from string import Template
from rdkit import Chem
from rdkit.Chem import Draw
import requests

size = (250, 250)

# Set Page Configuration
st.set_page_config(page_title="QAI Model - AI-Powered Quality Assistance", layout="wide", page_icon="🧪")

# Apply Custom Styles
st.markdown("""
    <style>
        @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap');

        body { 
            font-family: 'Inter', sans-serif;
            background-color: #E3F2FD; /* Light sky blue */
            color: #0B3D91;
        }
        
        /* Header */
        .header {
            background: linear-gradient(90deg, #0B3D91 0%, #1E4D9E 100%);
            padding: 2rem;
            border-radius: 12px;
            box-shadow: 0 4px 20px rgba(11,61,145,0.15);
            color: white;
            text-align: center;
            font-size: 24px;
            font-weight: 700;
        }

        /* Input Fields */
        .stTextInput > div > div > input,
        .stSelectbox > div > div > select,
        .stTextArea > div > textarea {
            background-color: #FFFFFF;
            border: 2px solid #007BFF;
            border-radius: 12px;
            padding: 12px;
            font-size: 16px;
            color: #0B3D91;
            transition: all 0.2s ease;
        }

        .stTextInput > div > div > input:focus,
        .stSelectbox > div > div > select:focus,
        .stTextArea > div > textarea:focus {
            border-color: #004085;
            box-shadow: 0 0 6px rgba(11,61,145,0.2);
        }

        /* Button Styling */
        .stButton > button {
            background: linear-gradient(90deg, #00B4DB, #0083B0);
            color: white;
            padding: 14px;
            border-radius: 8px;
            font-weight: 600;
            transition: all 0.3s ease;
            border: none;
            width: 100%;
            text-transform: uppercase;
            font-size: 16px;
        }

        .stButton > button:hover {
            background: linear-gradient(90deg, #007BFF, #004085);
            transform: scale(1.05);
            box-shadow: 0 8px 20px rgba(11,61,145,0.2);
        }

        /* Form Section */
        .form-container {
            background: white;
            padding: 2rem;
            border-radius: 12px;
            box-shadow: 0 6px 20px rgba(0,0,0,0.08);
            border: 1px solid #E2E8F0;
            margin-bottom: 1.5rem;
        }

        /* Results Table */
        .results-table {
            border-collapse: collapse;
            width: 100%;
            border-radius: 8px;
            box-shadow: 0 4px 15px rgba(0,0,0,0.08);
        }

        .results-table th {
            background: #0B3D91;
            color: white;
            padding: 1rem;
            text-align: left;
            font-weight: 600;
        }

        .results-table td {
            padding: 1rem;
            border-bottom: 1px solid #E2E8F0;
            background: white;
        }

        .results-table tr:hover td {
            background-color: #F8FAFC;
        }

        /* Spinner Animation */
        .stSpinner {
            border: 4px solid rgba(11,61,145,0.1);
            border-top: 4px solid #0B3D91;
            border-radius: 50%;
            width: 40px;
            height: 40px;
            animation: spin 1s linear infinite;
        }

        @keyframes spin {
            0% { transform: rotate(0deg); }
            100% { transform: rotate(360deg); }
        }

    </style>
""", unsafe_allow_html=True)

# Store session state
if "page" not in st.session_state:
    st.session_state.page = "form"
if "api_response" not in st.session_state:
    st.session_state.api_response = None

options = dict()

# 📌 FORM PAGE
if st.session_state.page == "form":
    # Header
    st.markdown('<div class="header">🧪 QAI Model - AI-Powered Quality Assurance</div>', unsafe_allow_html=True)

    # Form Section
    with st.container():
        col1, col2 = st.columns(2)

        with col1:
            st.markdown('<div class="form-container">', unsafe_allow_html=True)
            options["product_name"] = st.text_input("💊 Product Name", placeholder="e.g., Paracetamol")
            st.markdown('</div>', unsafe_allow_html=True)

        with col2:
            st.markdown('<div class="form-container">', unsafe_allow_html=True)
            options["quanOfMed"] = st.text_input("📦 Quantity", placeholder="e.g., 1000 tablets")
            options["powerOfDrug"] = st.text_input("⚡ Strength", placeholder="e.g., 500 mg")
            options["jurisdiction"] = st.selectbox("🌎 Pharmacopoeia Reference",
                ["INDIAN PHARMACOPIEA", "BRITISH PHARMACOPIEA", "UNITED STATES PHARMACOPOEIA",
                 "MARTINDALE-EXTRA PHARMACOPIEA", "COMPARE WITH ALL"])
            st.markdown('</div>', unsafe_allow_html=True)

    # Analysis Options
    st.markdown('<div class="form-container">', unsafe_allow_html=True)
    options["typeOfInfo"] = st.selectbox("📊 Analysis Type",
        ["METHOD OF PREPARATION", "CHARACTARIZATION/EVALUATION", "Both of above", "CHECK RESULTS"])

    if options["typeOfInfo"] == "CHECK RESULTS":
        options["resultsToCheck"] = st.text_area("🔍 Laboratory Results", height=200, placeholder="Enter lab results here...")

    options["ftir_required"] = st.checkbox("📡 Include FTIR Analysis")
    
    if st.button("🚀 Generate Report"):
        if all([options["product_name"], options["quanOfMed"], options["powerOfDrug"]]):
            with st.spinner("Analyzing data and generating report..."):
                prompt = prompts.getPromptForOptions(options)
                api_response = chat_with_gpt.chatWithGpt(prompt)
                st.session_state.api_response = api_response
                st.session_state.page = "result"
                st.experimental_rerun()
        else:
            st.error("⚠️ Please fill in all required fields.")
    st.markdown('</div>', unsafe_allow_html=True)

# 📌 RESULT PAGE
elif st.session_state.page == "result":
    if st.button("🔙 Back to Form"):
        st.session_state.page = "form"
        st.experimental_rerun()

    st.markdown('<div class="header">📑 Analysis Report</div>', unsafe_allow_html=True)

    if st.session_state.api_response:
        st.markdown('<div class="form-container">', unsafe_allow_html=True)
        st.markdown(st.session_state.api_response, unsafe_allow_html=True)
        st.markdown('</div>', unsafe_allow_html=True)
    else:
        st.error("⚠️ No analysis results available.")
