import streamlit as st
import streamlit.components.v1 as components
import prompts
import chat_with_gpt
from string import Template
import requests
from datetime import datetime

# Initialize session state for API key status
if 'api_key_status' not in st.session_state:
    st.session_state.api_key_status = False

# Add this function after the imports
def check_api_key():
    """Check if OpenAI API key is configured and valid."""
    try:
        if ("api" not in st.secrets or 
            "key" not in st.secrets["api"] or 
            not st.secrets["api"]["key"] or 
            not st.secrets["api"]["key"].strip()):
            return False
        return True
    except Exception:
        return False

# Set Page Configuration
st.set_page_config(
    page_title="QAI Model - Pharmaceutical Quality Assurance",
    layout="wide",
    page_icon="🧪",
    initial_sidebar_state="expanded"
)

# Apply Professional Styling
st.markdown("""
<style>
    /* Global Styles */
    body {
        background-color: #F8F9FA;
        color: #0B3D91;
        font-family: 'Helvetica Neue', Arial, sans-serif;
    }

    /* Header Styling */
    .main-header {
        background: linear-gradient(135deg, #0B3D91, #1E88E5);
        color: white;
        padding: 2rem;
        border-radius: 15px;
        margin-bottom: 2rem;
        text-align: center;
        box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
    }

    .main-title {
        font-size: 2.5rem;
        font-weight: 700;
        margin-bottom: 1rem;
        text-transform: uppercase;
        letter-spacing: 2px;
    }

    .subtitle {
        font-size: 1.2rem;
        opacity: 0.9;
        font-weight: 300;
    }

    /* Form Elements */
    .stTextInput>div>div>input,
    .stSelectbox>div>div>select,
    .stTextArea>div>textarea {
        border-radius: 10px !important;
        border: 2px solid #E3F2FD !important;
        padding: 1rem !important;
        background-color: white !important;
        transition: all 0.3s ease !important;
        font-size: 1rem !important;
    }

    .stTextInput>div>div>input:focus,
    .stSelectbox>div>div>select:focus,
    .stTextArea>div>textarea:focus {
        border-color: #1E88E5 !important;
        box-shadow: 0 0 0 2px rgba(30, 136, 229, 0.2) !important;
    }

    /* Button Styling */
    .stButton>button {
        background: linear-gradient(135deg, #0B3D91, #1E88E5);
        color: white;
        border: none;
        padding: 0.75rem 2rem;
        border-radius: 10px;
        font-weight: 600;
        transition: all 0.3s ease;
        width: 100%;
        text-transform: uppercase;
        letter-spacing: 1px;
        margin: 1rem 0;
    }

    .stButton>button:hover {
        transform: translateY(-2px);
        box-shadow: 0 4px 12px rgba(30, 136, 229, 0.3);
    }

    /* Card Styling */
    .card {
        background: white;
        border-radius: 15px;
        padding: 2rem;
        margin: 1rem 0;
        box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
        border-left: 5px solid #0B3D91;
    }

    /* Status Messages */
    .success-message {
        background-color: #E8F5E9;
        color: #2E7D32;
        padding: 1rem;
        border-radius: 10px;
        margin: 1rem 0;
        border-left: 5px solid #2E7D32;
    }

    .info-message {
        background-color: #E3F2FD;
        color: #1565C0;
        padding: 1rem;
        border-radius: 10px;
        margin: 1rem 0;
        border-left: 5px solid #1565C0;
    }

    .warning-message {
        background-color: #FFF3E0;
        color: #EF6C00;
        padding: 1rem;
        border-radius: 10px;
        margin: 1rem 0;
        border-left: 5px solid #EF6C00;
    }

    .error-message {
        background-color: #FFEBEE;
        color: #C62828;
        padding: 1rem;
        border-radius: 10px;
        margin: 1rem 0;
        border-left: 5px solid #C62828;
    }

    /* Loading Spinner */
    .stSpinner {
        text-align: center;
        padding: 2rem;
    }

    /* Tooltip */
    .tooltip {
        position: relative;
        display: inline-block;
        cursor: help;
    }

    .tooltip .tooltiptext {
        visibility: hidden;
        width: 200px;
        background-color: #0B3D91;
        color: white;
        text-align: center;
        border-radius: 6px;
        padding: 5px;
        position: absolute;
        z-index: 1;
        bottom: 125%;
        left: 50%;
        margin-left: -100px;
        opacity: 0;
        transition: opacity 0.3s;
    }

    .tooltip:hover .tooltiptext {
        visibility: visible;
        opacity: 1;
    }
</style>
""", unsafe_allow_html=True)

# Initialize Session State
if "page" not in st.session_state:
    st.session_state.page = "form"
if "api_response" not in st.session_state:
    st.session_state.api_response = None

# Helper Functions
def get_pubchem_structure_url(drug_name):
    """Get the 2D structure image URL from PubChem."""
    if not drug_name:
        st.warning("⚠️ Please enter a drug name first.")
        return None

    try:
        search_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{drug_name}/cids/JSON"
        response = requests.get(search_url, timeout=10)
        if response.status_code == 200:
            try:
                cid = response.json()["IdentifierList"]["CID"][0]
                return f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/{cid}/PNG"
            except (KeyError, IndexError):
                st.error("❌ Could not find structure for this drug.")
                return None
        else:
            st.error("❌ Error connecting to PubChem.")
            return None
    except Exception as e:
        st.error(f"❌ Error fetching structure: {str(e)}")
        return None

# Main Application Logic
if st.session_state.page == "form":
    # Header
    st.markdown("""
        <div class="main-header">
            <h1 class="main-title">🧪 QAI Model</h1>
            <p class="subtitle">AI-Powered Pharmaceutical Quality Assurance</p>
        </div>
    """, unsafe_allow_html=True)

    # Information Message
    st.markdown("""
        <div class="info-message">
            <strong>👋 Welcome!</strong> This tool helps you generate comprehensive pharmaceutical quality reports.
            Simply fill in the details below and let our AI assistant analyze your requirements.
        </div>
    """, unsafe_allow_html=True)

    # Add this at the start of the form page (after the header)
    st.session_state.api_key_status = check_api_key()

    # Form Layout
    col1, col2 = st.columns([2, 1])

    with col1:
        with st.container():
            st.markdown('<div class="card">', unsafe_allow_html=True)
            options = {}
            options["product_name"] = st.text_input("💊 Product Name",
                placeholder="e.g., Paracetamol",
                help="Enter the name of the pharmaceutical product")

            if st.button("View Molecular Structure", key="structure_button"):
                if not options["product_name"]:
                    st.warning("⚠️ Please enter a product name first.")
                else:
                    with st.spinner("🔍 Fetching molecular structure..."):
                        structure_url = get_pubchem_structure_url(options["product_name"])
                        if structure_url:
                            st.success("✅ Structure found!")
                            st.image(structure_url, caption=f"{options['product_name']} Molecular Structure")

            options["quanOfMed"] = st.text_input("📦 Quantity of Medicine",
                placeholder="e.g., 1000 tablets",
                help="Specify the total quantity to be manufactured")

            options["powerOfDrug"] = st.text_input("⚡ Power of Drug",
                placeholder="e.g., 500 mg",
                help="Enter the strength of each unit")
            st.markdown('</div>', unsafe_allow_html=True)

    with col2:
        with st.container():
            st.markdown('<div class="card">', unsafe_allow_html=True)
            options["jurisdiction"] = st.selectbox("🌎 Select Jurisdiction",
                ["INDIAN PHARMACOPIEA",
                 "BRITISH PHARMACOPIEA",
                 "UNITED STATES PHARMACOPOEIA",
                 "MARTINDALE-EXTRA PHARMACOPIEA",
                 "COMPARE WITH ALL"],
                help="Choose the regulatory standard to follow")

            options["typeOfInfo"] = st.selectbox("📊 Select Information Required:",
                ["METHOD OF PREPARATION",
                 "CHARACTARIZATION/EVALUATION",
                 "Both of above",
                 "CHECK RESULTS"],
                help="Choose the type of report you need")

            if options["typeOfInfo"] == "CHECK RESULTS":
                options["resultsToCheck"] = st.text_area("🔍 Enter Your Results:",
                    height=200,
                    placeholder="Paste lab results here...",
                    key="checkResults",
                    help="Enter your test results for analysis")

            options["ftir_required"] = st.checkbox("📡 Include FTIR Data Analysis",
                help="Include detailed FTIR spectral analysis in the report")
            st.markdown('</div>', unsafe_allow_html=True)

    # Submit Button
    st.markdown('<div class="card">', unsafe_allow_html=True)
    if st.button("🚀 Generate Quality Report", type="primary"):
        api_key_configured = check_api_key()
        if not api_key_configured:
            st.error("❌ OpenAI API key is not configured. Please check your configuration.")
        elif not all([options["product_name"], options["quanOfMed"], options["powerOfDrug"]]):
            st.error("⚠️ Please fill in all required fields!")
        else:
            with st.spinner("🔄 Generating comprehensive report..."):
                prompt = prompts.getPromptForOptions(options)
                api_response = chat_with_gpt.chatWithGpt(prompt)

                if "Error:" in api_response:
                    if "API key" in api_response:
                        st.error("🔑 OpenAI API key error. Please check your configuration.")
                    else:
                        st.error(api_response)
                else:
                    st.success("✅ Report generated successfully!")
                    st.session_state.api_response = api_response
                    st.session_state.update(options)
                    st.session_state.page = "result"
                    st.experimental_rerun()
    st.markdown('</div>', unsafe_allow_html=True)

elif st.session_state.page == "result":
    # Results Page Header
    st.markdown("""
        <div class="main-header">
            <h1 class="main-title">📊 Quality Analysis Report</h1>
            <p class="subtitle">Comprehensive Evaluation Results</p>
        </div>
    """, unsafe_allow_html=True)

    # Navigation Buttons
    col1, col2 = st.columns([1, 4])
    with col1:
        if st.button("← Back to Form"):
            st.session_state.page = "form"
            st.experimental_rerun()

    # Display Results
    with st.container():
        st.markdown('<div class="card">', unsafe_allow_html=True)
        st.markdown(f"""
            <div style='margin-bottom: 1rem;'>
                <strong>💊 Product:</strong> {st.session_state.product_name}<br>
                <strong>📦 Quantity:</strong> {st.session_state.quanOfMed}<br>
                <strong>⚡ Strength:</strong> {st.session_state.powerOfDrug}
            </div>
        """, unsafe_allow_html=True)

        # Download Button
        if st.session_state.api_response:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            filename = f"QA_Report_{st.session_state.product_name}_{timestamp}.html"

            # Create downloadable HTML content
            download_content = f"""
            <html>
            <head>
                <style>
                    body {{
                        font-family: Arial, sans-serif;
                        margin: 2rem;
                        color: #333;
                    }}
                    .header {{
                        text-align: center;
                        margin-bottom: 2rem;
                    }}
                    .report-info {{
                        margin-bottom: 2rem;
                        padding: 1rem;
                        background-color: #f8f9fa;
                        border-radius: 8px;
                    }}
                </style>
            </head>
            <body>
                <div class="header">
                    <h1>Quality Analysis Report</h1>
                    <p>Generated on: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</p>
                </div>
                <div class="report-info">
                    <p><strong>Product:</strong> {st.session_state.product_name}</p>
                    <p><strong>Quantity:</strong> {st.session_state.quanOfMed}</p>
                    <p><strong>Strength:</strong> {st.session_state.powerOfDrug}</p>
                </div>
                {st.session_state.api_response}
            </body>
            </html>
            """

            st.download_button(
                label="📥 Download Report",
                data=download_content,
                file_name=filename,
                mime="text/html",
                help="Download the report as an HTML file"
            )
        st.markdown('</div>', unsafe_allow_html=True)

        if st.session_state.api_response:
            st.markdown(st.session_state.api_response, unsafe_allow_html=True)

            if st.session_state.ftir_required:
                with st.spinner("📡 Analyzing FTIR Data..."):
                    ftir_data = chat_with_gpt.get_ftir_from_gpt(st.session_state.product_name)
                    if "Error:" not in ftir_data:
                        st.markdown("### 🔬 FTIR Analysis")
                        st.markdown(ftir_data, unsafe_allow_html=True)
                    else:
                        st.error(ftir_data)
        else:
            st.error("⚠️ No response received. Please try again.")