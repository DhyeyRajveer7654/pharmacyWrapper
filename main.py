import prompts
import chat_with_gpt
from string import Template
import os

# Configure the page
st.set_page_config(
    page_title="QRx - Pharmaceutical Quality & Regulatory Experts",
    page_icon="💊",
    layout="wide",
    initial_sidebar_state="collapsed"
)

# Initialize session state for page tracking
if 'current_page' not in st.session_state:
    st.session_state.current_page = 'home'
# Not showing popup by default
if 'show_popup' not in st.session_state:
    st.session_state.show_popup = False

# Custom CSS styling
st.markdown("""
<style>
    /* Global styles */
    body {
        font-family: 'Segoe UI', sans-serif;
        color: #333;
        background-color: #f9f9f9;
    }
    
    /* Header styling */
    .header {
        background-color: white;
        color: #1e40af;
        padding: 1rem 2rem;
        display: flex;
        justify-content: space-between;
        align-items: center;
        border-bottom: 1px solid #e0e0e0;
        margin-bottom: 2rem;
    }
    
    .header h1 {
        font-size: 2.5rem;
        font-weight: 700;
        margin: 0;
    }
    
    .nav-links {
        display: flex;
        gap: 1.5rem;
    }
    
    .nav-link {
        color: #334155;
        text-decoration: none;
        font-weight: 600;
        font-size: 1rem;
        letter-spacing: 0.5px;
    }
    
    .nav-link:hover {
        color: #1e40af;
    }
    
    /* Hero section */
    .hero {
        background: linear-gradient(135deg, #1e40af, #3b82f6);
        color: white;
        padding: 3rem 2rem;
        border-radius: 10px;
        text-align: center;
        margin-bottom: 2rem;
    }
    
    .hero h2 {
        font-size: 2rem;
        margin-bottom: 1rem;
    }
    
    .hero p {
        font-size: 1.1rem;
        margin: 0 auto;
        max-width: 800px;
    }
    
    /* Card styling */
    .card {
        background-color: white;
        border-radius: 10px;
        padding: 1.5rem;
        box-shadow: 0 4px 6px rgba(0,0,0,0.1);
        margin-bottom: 1.5rem;
    }
    
    .card h3 {
        color: #1e40af;
        margin-bottom: 1rem;
    }
    
    /* Section titles */
    .section-title {
        color: #1e40af;
        margin: 2rem 0 1rem 0;
        text-align: center;
        font-size: 1.8rem;
    }
    
    /* Footer */
    .footer {
        background-color: #1e3a8a;
        color: white;
        padding: 2rem;
        border-radius: 10px 10px 0 0;
        margin-top: 3rem;
    }
    
    .footer h3 {
        color: #60a5fa;
        margin-bottom: 1rem;
    }
    
    .footer-bottom {
        text-align: center;
        margin-top: 2rem;
        padding-top: 1rem;
        border-top: 1px solid #334155;
        color: #94a3b8;
    }
    
    /* Popup styling */
    .popup-container {
        background-color: white;
        border-radius: 10px;
        box-shadow: 0 0 20px rgba(0, 0, 0, 0.2);
        padding: 2rem;
        max-width: 500px;
        margin: 0 auto 2rem auto;
        position: relative;
    }
    
    .popup-container h3 {
        color: #1e40af;
        text-align: center;
        margin-bottom: 1.5rem;
    }
    
    .popup-close {
        position: absolute;
        top: 10px;
        right: 15px;
        cursor: pointer;
        font-size: 1.5rem;
        color: #94a3b8;
    }
    
    /* Override Streamlit elements */
    .stButton > button {
        background-color: #1e40af;
        color: white;
        font-weight: 600;
        border: none;
        width: 100%;
    }
    
    .stButton > button:hover {
        background-color: #2563eb;
    }
    
    div[data-testid="stToolbar"] {
        display: none;
    }
    
    div[data-testid="stDecoration"] {
        display: none;
    }
    
    section[data-testid="stSidebar"] {
        display: none;
    }
    
    #MainMenu {
        display: none;
    }
    
    footer {
        display: none;
    }
</style>
""", unsafe_allow_html=True)

# Custom header with navigation
header_html = """
<div class="header">
    <h1>QRx AI</h1>
    <div class="nav-links">
        <a href="#" class="nav-link" id="home-link">HOME</a>
        <a href="#" class="nav-link" id="services-link">SERVICES</a>
        <a href="#" class="nav-link" id="contact-link">CONTACT US</a>
        <a href="#" class="nav-link" id="about-link">ABOUT US</a>
    </div>
</div>
<script>
    // Navigation handling using JavaScript
    document.addEventListener('DOMContentLoaded', function() {
        document.getElementById('home-link').addEventListener('click', function(e) {
            e.preventDefault();
            window.parent.postMessage({
                type: "streamlit:setComponentValue",
                value: {"page": "home"}
            }, "*");
        });
        
        document.getElementById('services-link').addEventListener('click', function(e) {
            e.preventDefault();
            window.parent.postMessage({
                type: "streamlit:setComponentValue",
                value: {"page": "services"}
            }, "*");
        });
        
        document.getElementById('contact-link').addEventListener('click', function(e) {
            e.preventDefault();
            window.parent.postMessage({
                type: "streamlit:setComponentValue",
                value: {"page": "contact"}
            }, "*");
        });
        
        document.getElementById('about-link').addEventListener('click', function(e) {
            e.preventDefault();
            window.parent.postMessage({
                type: "streamlit:setComponentValue",
                value: {"page": "about"}
            }, "*");
        });
    });
</script>
"""
st.markdown(header_html, unsafe_allow_html=True)

# Navigation state handler
if st.checkbox("Navigation State", key="nav_state_handler", value=False, label_visibility="collapsed"):
    clicked_page = st.session_state.get("nav_state_handler", None)
    st.session_state.nav_state_handler = False

# Function to close popup
def close_popup():
    st.session_state.show_popup = False
    st.rerun()

# Function to change page
def change_page(page):
    st.session_state.current_page = page
    st.session_state.show_popup = False
    st.rerun()

# Function to handle JavaScript close click
def handle_js_close():
    # This will be called after the JS close action to update session state
    st.session_state.show_popup = False
    st.rerun()

# Popup implementation
if st.session_state.show_popup:
    # Insert custom CSS for the popup
    st.markdown("""
    <style>
    .popup-overlay {
        position: fixed;
        top: 0;
        left: 0;
        width: 100%;
        height: 100%;
        background-color: rgba(0, 0, 0, 0.5);
        z-index: 1000;
        display: flex;
        justify-content: center;
        align-items: center;
    }
    
    .popup-container {
        background-color: white;
        width: 90%;
        max-width: 600px;
        border-radius: 10px;
        box-shadow: 0 0 30px rgba(0, 0, 0, 0.2);
        position: relative;
        padding: 2rem;
    }
    
    .popup-title {
        color: #1e40af;
        text-align: center;
        margin-bottom: 2rem;
        font-size: 1.8rem;
        font-weight: 600;
    }
    
    .popup-close {
        position: absolute;
        top: 15px;
        right: 20px;
        font-size: 28px;
        font-weight: bold;
        color: #94a3b8;
        cursor: pointer;
        line-height: 1;
    }
    
    .popup-close:hover {
        color: #1e40af;
    }
    
    .popup-buttons {
        display: flex;
        gap: 20px;
        margin-top: 1.5rem;
    }
    
    .popup-button {
        flex: 1;
        background: linear-gradient(135deg, #1e40af, #3b82f6);
        color: white;
        border: none;
        padding: 12px 20px;
        border-radius: 8px;
        font-weight: 600;
        cursor: pointer;
        text-align: center;
        transition: transform 0.2s, box-shadow 0.2s;
    }
    
    .popup-button:hover {
        transform: translateY(-2px);
        box-shadow: 0 4px 12px rgba(59, 130, 246, 0.3);
    }
    
    /* Style the popup buttons to match specification */
    .stButton > button {
        background-color: #1e40af !important;
        color: white !important;
        font-weight: 600 !important;
        border: none !important;
        border-radius: 6px !important;
        transition: background-color 0.3s, transform 0.2s !important;
        padding: 0.5rem 1rem !important;
        text-transform: uppercase !important;
        letter-spacing: 0.5px !important;
    }
    
    .stButton > button:hover {
        background-color: #2563eb !important;
        transform: translateY(-2px) !important;
        box-shadow: 0 4px 8px rgba(37, 99, 235, 0.3) !important;
    }
    </style>
    """, unsafe_allow_html=True)
    
    # Create container for popup content
    popup_container = st.container()
    
    # Create popup overlay and content with close button
    st.markdown("""
    <div class="popup-overlay">
        <div class="popup-container">
            <div class="popup-close" id="closeButton">×</div>
            <h2 class="popup-title">Choose an Option</h2>
            <div id="popup-content">
                <!-- The buttons will be rendered by Streamlit below -->
            </div>
        </div>
    </div>
    <script>
        // Add event listener to close button
        document.addEventListener('DOMContentLoaded', function() {
            document.getElementById('closeButton').addEventListener('click', function() {
                document.querySelector('.popup-overlay').style.display = 'none';
                // After hiding, submit a form to trigger the Streamlit rerun
                setTimeout(function() {
                    window.parent.postMessage({
                        type: "streamlit:setComponentValue",
                        value: true
                    }, "*");
                }, 100);
            });
        });
    </script>
    """, unsafe_allow_html=True)
    
    # Add a callback for JavaScript close button
    if st.checkbox("JS Close Triggered", key="js_close_trigger", value=False, label_visibility="collapsed"):
        close_popup()
    
    # Add the buttons with more prominent styling
    st.markdown("<div style='padding: 20px;'></div>", unsafe_allow_html=True)
    col1, col2 = st.columns(2)
    with col1:
        if st.button("OPTION 1: REGULATORY COMPLIANCE", key="popup_option1", use_container_width=True):
            change_page('regulatory')
    with col2:
        if st.button("OPTION 2: QUALITY ASSURANCE", key="popup_option2", use_container_width=True):
            change_page('quality')
            
    # Add close button under the options
    # This is a backup in case the × close button doesn't work
    if st.button("Close", key="close_popup"):
        close_popup()

# Home page content
if st.session_state.current_page == 'home':
    # Hero section
    st.markdown("""
    <div class="hero">
        <h2>Welcome to QRx AI</h2>
        <p>Your trusted partner for pharmaceutical quality analysis and regulatory compliance. 
        We combine cutting-edge technology with expert knowledge to ensure your products meet the highest standards.</p>
    </div>
    """, unsafe_allow_html=True)
    
    # Introduction content
    st.markdown('<h3 class="section-title">Introduction to QRx AI</h3>', unsafe_allow_html=True)
    
    st.markdown("""
    <div class="card">
        <p>QRx AI is a leading provider of pharmaceutical quality and regulatory solutions. We help pharmaceutical companies navigate complex regulatory requirements while ensuring their products meet the highest quality standards.</p>
        <p>Our team of experts combines decades of industry experience with cutting-edge technology to deliver comprehensive solutions tailored to your specific needs.</p>
        <p>Whether you need help with regulatory compliance or quality assurance, we have the expertise and resources to support you at every stage of your product lifecycle.</p>
    </div>
    """, unsafe_allow_html=True)
    
    # Services section
    st.markdown('<h3 class="section-title">Our Services</h3>', unsafe_allow_html=True)
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("""
        <div class="card">
            <h3>Regulatory Compliance</h3>
            <p>Navigate complex regulatory landscapes with our expert guidance. We ensure your products meet all requirements for FDA, EMA, and other global regulatory bodies.</p>
        </div>
        """, unsafe_allow_html=True)
        if st.button("Learn More - Regulatory", key="reg_button"):
            change_page('regulatory')
    
    with col2:
        st.markdown("""
        <div class="card">
            <h3>Quality Assurance</h3>
            <p>Comprehensive quality analysis using AI-powered tools to evaluate pharmaceutical products against pharmacopeial standards.</p>
        </div>
        """, unsafe_allow_html=True)
        if st.button("Learn More - Quality", key="qa_button"):
            change_page('quality')
    
    
# Regulatory Compliance page
elif st.session_state.current_page == 'regulatory':
    st.markdown("""
    <div class="hero">
        <h2>Regulatory Compliance Services</h2>
        <p>Navigate complex regulatory landscapes with confidence. Our experts ensure your products meet all global standards and requirements.</p>
    </div>
    """, unsafe_allow_html=True)
    
    # Back button
    if st.button("← Back to Home", key="reg_back"):
        change_page('home')
    
    # Content
    st.markdown('<h3 class="section-title">FDA Compliance Services</h3>', unsafe_allow_html=True)
    
    st.markdown("""
    <div class="card">
        <p>Our FDA compliance services help pharmaceutical companies navigate the complex regulatory requirements of the U.S. Food and Drug Administration.</p>
        <p>We provide comprehensive guidance on:</p>
        <ul>
            <li>Drug application submissions (NDA, ANDA, BLA)</li>
            <li>GMP compliance</li>
            <li>Manufacturing practices</li>
            <li>Labeling requirements</li>
            <li>Post-market surveillance</li>
        </ul>
    </div>
    """, unsafe_allow_html=True)
    
    # Additional services
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("""
        <div class="card">
            <h3>Global Regulatory Strategy</h3>
            <p>Develop comprehensive regulatory strategies tailored to your products and target markets. We help navigate requirements for FDA, EMA, MHRA, TGA, and other global regulatory bodies.</p>
        </div>
        """, unsafe_allow_html=True)
    
    with col2:
        st.markdown("""
        <div class="card">
            <h3>Compliance Audits</h3>
            <p>Thorough audits of your facilities and processes to ensure compliance with GMP, GLP, GCP, and other regulatory requirements.</p>
        </div>
        """, unsafe_allow_html=True)

# Quality Assurance page
elif st.session_state.current_page == 'quality':
    # Import required libraries for QA tool functionality
    
    import streamlit as st # type: ignore
import streamlit.components.v1 as components # type: ignore
import prompts
import chat_with_gpt
from string import Template
from rdkit import Chem # type: ignore # type: ignore
from rdkit.Chem import Draw # type: ignore
import requests # type: ignore
import os
import streamlit as st # type: ignore

size = (250, 250)

# Set Page Configuration
st.set_page_config(page_title="QAI Model AI-Powered Quality Assistance", layout="wide", page_icon="🧪")

# Enhanced Professional UI Styling
st.markdown("""
    <style>
        /* Global Styles */
        body {
            background-color: #f0f2f6;
            color: #1e293b;
            font-family: 'Inter', 'sans serif';
        }

        /* Header Styling */
        .main-header {
            background: linear-gradient(135deg, #0052cc, #00a3bf);
            color: white;
            padding: 2rem;
            border-radius: 10px;
            margin-bottom: 2rem;
            text-align: center;
        }

        /* Form Elements */
        div[data-testid="stTextInput"] input,
        div[data-testid="stTextArea"] textarea,
        div[data-testid="stSelectbox"] > div[data-baseweb="select"] {
            background-color: white;
            border: 1px solid #e2e8f0;
            border-radius: 8px;
            padding: 0.75rem;
            font-size: 1rem;
            transition: all 0.3s ease;
            box-shadow: 0 2px 4px rgba(0,0,0,0.05);
        }

        div[data-testid="stTextInput"] input:focus,
        div[data-testid="stTextArea"] textarea:focus {
            border-color: #0052cc;
            box-shadow: 0 0 0 2px rgba(0,82,204,0.2);
        }

        /* Button Styling */
        .stButton > button {
            width: 100%;
            background: linear-gradient(135deg, #0052cc, #00a3bf);
            color: white;
            border: none;
            padding: 0.75rem 1.5rem;
            border-radius: 8px;
            font-weight: 600;
            transition: all 0.3s ease;
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }

        .stButton > button:hover {
            transform: translateY(-2px);
            box-shadow: 0 4px 12px rgba(0,82,204,0.2);
        }

        /* Card Styling */
        .card div[data-testid="stSelectbox"]{
            background: white;
            border-radius: 10px;
            padding: 1.5rem;
            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
            margin-bottom: 1.5rem;
        }

        /* Results Table */
        .table-container {
            background: white;
            border-radius: 10px;
            padding: 1rem;
            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
            border-collapse: collapse;
        }

        table {
            width: 100%;
            border-collapse: collapse;
            border-spacing: 0;
            margin: 1rem 0;
        }

        th {
            background: #0052cc;
            color: white;
            padding: 1rem;
            text-align: left;
            font-weight: 600;
        }

        td {
            padding: 1rem;
            border-bottom: 1px solid #e2e8f0;
        }

        tr:hover {
            background: #f8fafc;
        }

        /* Loading Spinner */
        .stSpinner > div {
            border-color: #0052cc !important;
        }

        /* Success Message */
        .success-message {
            background: #dcfce7;
            color: #166534;
            padding: 1rem;
            border-radius: 8px;
            margin: 1rem 0;
        }

        /* Error Message */
        .error-message {
            background: #fee2e2;
            color: #991b1b;
            padding: 1rem;
            border-radius: 8px;
            margin: 1rem 0;
        }
        
    </style>
""", unsafe_allow_html=True)

# Page Navigation
if "page" not in st.session_state:
    st.session_state.page = "form"
if "api_response" not in st.session_state:
    st.session_state.api_response = None

options = dict()

def get_cid_from_name(drug_name):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{drug_name}/cids/JSON"
    response = requests.get(url)

    if response.status_code == 200:
        try:
            cids = response.json()["IdentifierList"]["CID"]
            return cids[0]  # Return the first matching CID
        except (KeyError, IndexError):
            return None
    else:
        return None

def get_pubchem_product_code(product_name):
    product_code_from_pubchem = ""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{product_name}/property/CanonicalSMILES/JSON"
    response = requests.get(url)
    if response.status_code == 200:
        try:
            smiles = response.json()["PropertyTable"]["Properties"][0]["CanonicalSMILES"]
            product_code_from_pubchem=smiles
        except (KeyError, IndexError):
            product_code_from_pubchem = "NO DRUG FOUND"
    else:
        product_code_from_pubchem="NO DRUG FOUND"
    if product_code_from_pubchem=="NO DRUG FOUND":
        return ""
    else:
        return product_code_from_pubchem

def showStructure(product_name):
    product_code = ""
    product_code_from_pubchem = get_pubchem_product_code(product_name)
    if product_code_from_pubchem=="":
        product_code_prompt = prompts.STRUCTURE_PROMPT.substitute(product_name=product_name)
        print("Prompt is: "+product_code_prompt)
        product_code = chat_with_gpt.chatWithGpt(product_code_prompt)
        if product_code == "NO DRUG FOUND":
            return ""
    else:
        product_code = product_code_from_pubchem

    print("product code is: "+product_code)
    print("product code from pubchem: "+product_code_from_pubchem)
    m = Chem.MolFromSmiles(product_code)
    if m:
        return Draw.MolToImage(m, size=(400, 400))
    return None

# Directory where FTIR images are stored
FTIR_IMAGE_DIR = "./"

def get_ftir_image(product_name):
    """Fetches the corresponding FTIR image for the given product name."""
    image_filename = f"{product_name.lower()}.png"
    image_path = os.path.join(FTIR_IMAGE_DIR, image_filename)
    if os.path.exists(image_path):
        return image_path
    return None

# 📌 FORM PAGE
if st.session_state.page == "form":
    st.markdown('<div class="main-header"><h1>🧪 QAI Model AI-Powered Quality Assistance</h1><p> CREATED BY :- MEERA ACHARYA & RAJ PATEL</P><p>Enter details below to generate a comprehensive quality report</p></div>', unsafe_allow_html=True)

    # User Input Form in a card layout
    # st.markdown('<div class="card">', unsafe_allow_html=True)
    col1, col2 = st.columns(2)

    with col1:
        options["product_name"] = st.text_input("💊 Product Name", placeholder="e.g., Paracetamol")
        if st.button("🔬 Get Structure"):
            if not options["product_name"]:
                st.error("⚠️ Please write product name!")
            else:
                with st.spinner("🛠️ Processing... Please wait"):
                    fig = showStructure(options["product_name"])
                if fig == "":
                    st.error("⚠️ Drug not found, please input a valid drug name")
                else:
                    st.image(fig, caption=f"{options['product_name']} Molecule")

        if st.button("📊 Show FTIR Graph"):
            if options.get("product_name"):  # Ensure product name exists
                ftir_image = get_ftir_image(options["product_name"])
                if ftir_image:
                    st.image(ftir_image, caption=f"FTIR Graph for {options['product_name']}", use_column_width=True)
                else:
                    st.error(f"⚠️ No FTIR data available for {options['product_name']}.")
            else:
                st.error("⚠️ Please enter a product name.")            

    with col2:
        options["quanOfMed"] = st.text_input("📦 Quantity of Medicine", placeholder="e.g., 1000 tablets")
        options["jurisdiction"] = st.selectbox("🌎 Select Jurisdiction", 
            ["INDIAN PHARMACOPIEA", "BRITISH PHARMACOPIEA", "UNITED STATES PHARMACOPOEIA", "MARTINDALE-EXTRA PHARMACOPIEA", "COMPARE WITH ALL"])
        options["powerOfDrug"] = st.text_input("⚡ Power of Drug", placeholder="e.g., 500 mg")

    # st.markdown('</div>', unsafe_allow_html=True)

    # Analysis Options in a separate card
    # st.markdown('<div class="card">', unsafe_allow_html=True)
    options["typeOfInfo"] = st.selectbox("📊 Select Analysis Type:", 
            ["METHOD OF PREPARATION", "CHARACTARIZATION/EVALUATION", "Both of above", "CHECK RESULTS"])

    if options["typeOfInfo"] == "CHECK RESULTS":
        options["resultsToCheck"] = st.text_area("🔍 Enter Your Results:", height=200, placeholder="Paste lab results here...", key="checkResults")

    options["ftir_required"] = st.checkbox("📡 Include FTIR Analysis")
    # st.markdown('</div>', unsafe_allow_html=True)

    # Submit button with enhanced styling
    submit_button = st.button("🚀 Generate Report")
    if submit_button:
        if not all([options.get("product_name"), options.get("quanOfMed"), options.get("powerOfDrug")]):
            st.error("⚠️ Please fill in all required fields!")
        else:
            prompt = prompts.getPromptForOptions(options)
            with st.spinner("🛠️ Generating comprehensive report... Please wait"):
                api_response = chat_with_gpt.chatWithGpt(prompt)
                st.session_state.api_response = api_response

            st.session_state.update(options)
            st.session_state.page = "result"
            st.experimental_rerun()

# 📌 RESULT PAGE
elif st.session_state.page == "result":
    st.markdown('<div class="main-header"><h1>📑 Quality Analysis Report</h1></div>', unsafe_allow_html=True)
    
    if st.button("🔙 Return to Form", key="back_button"):
        st.session_state.page = "form"
        st.experimental_rerun()

    # st.markdown('<div class="card">', unsafe_allow_html=True)
    st.markdown("### 📋 Analysis Details")
    st.markdown(f"**💊 Product:** {st.session_state.product_name}",)
    st.markdown(f"**📦 Quantity:** {st.session_state.quanOfMed}")
    st.markdown(f"**⚡ Strength:** {st.session_state.powerOfDrug}")
    # st.markdown('</div>', unsafe_allow_ht ml=True)

    if st.session_state.get("ftir_required"):
        with st.spinner("📡 Analyzing FTIR Data..."):
            ftir_data = chat_with_gpt.get_ftir_from_gpt(st.session_state.product_name)
            components.html("### 🔬 FTIR Analysis")
            st.markdown(ftir_data, unsafe_allow_html=True)
            # components.html(ftir_data)

    if st.session_state.api_response:
        components.html("<div class='table-container'>"+st.session_state.api_response+"</div>",height=800,width=1000,scrolling=True)
    else:
        st.warning("⚠️ No response received. Please try again.")

    # Back button
    if st.button("← Back to Home", key="method_back"):
        change_page('home')
    
    
# Contact Us page
elif st.session_state.current_page == 'contact':
    st.markdown("""
    <div class="hero">
        <h2>Contact Us</h2>
        <p>Reach out to our team of experts with your questions or consultation requests.</p>
    </div>
    """, unsafe_allow_html=True)
    
    # Back button
    if st.button("← Back to Home", key="contact_back"):
        change_page('home')
    
    # Content
    st.markdown('<h3 class="section-title">Get in Touch</h3>', unsafe_allow_html=True)
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("""
        <div class="card">
            <h3>Contact Information</h3>
            <p><strong>Email:</strong> info@qrx-pharma.com</p>
            <p><strong>Phone:</strong> +1 (555) 123-4567</p>
            <p><strong>Address:</strong><br>
            QRx Headquarters<br>
            123 Pharmaceutical Lane<br>
            Boston, MA 02115<br>
            United States</p>
        </div>
        """, unsafe_allow_html=True)
    
    with col2:
        st.markdown("""
        <div class="card">
            <h3>Contact Form</h3>
        </div>
        """, unsafe_allow_html=True)
        
        with st.form("contact_form"):
            st.text_input("Full Name", placeholder="Enter your full name")
            st.text_input("Email", placeholder="Enter your email address")
            st.text_input("Company", placeholder="Enter your company name")
            st.selectbox("Subject", ["General Inquiry", "Regulatory Services", "Quality Assurance", "Other"])
            st.text_area("Message", placeholder="How can we help you?")
            
            if st.form_submit_button("Submit"):
                st.success("Thank you for your message! Our team will contact you shortly.")

# About Us page
elif st.session_state.current_page == 'about':
    st.markdown("""
    <div class="hero">
        <h2>About QRx</h2>
        <p>Learn about our mission, vision, and the expert team behind our services.</p>
    </div>
    """, unsafe_allow_html=True)
    
    # Back button
    if st.button("← Back to Home", key="about_back"):
        change_page('home')
    
    # Content
    st.markdown('<h3 class="section-title">Our Story</h3>', unsafe_allow_html=True)
    
    st.markdown("""
    <div class="card">
        <p>QRx was founded in 2015 by a team of pharmaceutical industry veterans who recognized the need for more efficient and technology-driven approaches to quality assurance and regulatory compliance.</p>
        <p>Today, we are a leading provider of pharmaceutical quality and regulatory solutions, serving clients across the globe. Our team combines decades of industry experience with cutting-edge technology to deliver comprehensive solutions tailored to our clients' specific needs.</p>
    </div>
    """, unsafe_allow_html=True)
    
    # Mission and vision
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("""
        <div class="card">
            <h3>Our Mission</h3>
            <p>To enhance pharmaceutical quality and safety through innovative solutions that combine expert knowledge with cutting-edge technology.</p>
        </div>
        """, unsafe_allow_html=True)
    
    with col2:
        st.markdown("""
        <div class="card">
            <h3>Our Vision</h3>
            <p>To be the global leader in pharmaceutical quality and regulatory solutions, recognized for our expertise, innovation, and commitment to excellence.</p>
        </div>
        """, unsafe_allow_html=True)
    
    # Team section
    st.markdown('<h3 class="section-title">Our Leadership Team</h3>', unsafe_allow_html=True)
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.markdown("""
        <div class="card">
            <h3>Dr. Sarah Johnson</h3>
            <p><em>CEO & Founder</em></p>
            <p>Dr. Johnson has over 20 years of experience in pharmaceutical quality and regulatory affairs, having held leadership positions at major pharmaceutical companies.</p>
        </div>
        """, unsafe_allow_html=True)
    
    with col2:
        st.markdown("""
        <div class="card">
            <h3>Dr. Michael Chen</h3>
            <p><em>CTO</em></p>
            <p>Dr. Chen leads our technology team, combining his expertise in pharmaceutical science with advanced data analytics and AI to develop our cutting-edge solutions.</p>
        </div>
        """, unsafe_allow_html=True)
    
    with col3:
        st.markdown("""
        <div class="card">
            <h3>Dr. Emily Rodriguez</h3>
            <p><em>Head of Regulatory Affairs</em></p>
            <p>Dr. Rodriguez brings extensive experience in global regulatory strategy, having successfully guided numerous products through approval processes worldwide.</p>
        </div>
        """, unsafe_allow_html=True)

# Footer
st.markdown("""
<div class="footer">
    <div style="display: flex; flex-wrap: wrap; gap: 2rem;">
        <div style="flex: 1 1 200px;">
            <h3>QRx</h3>
            <p>Your trusted partner for pharmaceutical quality analysis and regulatory compliance.</p>
        </div>
        
        <div style="flex: 1 1 200px;">
            <h3>Quick Links</h3>
            <ul>
                <li>Home</li>
                <li>Services</li>
                <li>About Us</li>
                <li>Contact Us</li>
            </ul>
        </div>
        
        <div style="flex: 1 1 200px;">
            <h3>Our Services</h3>
            <ul>
                <li>Regulatory Compliance</li>
                <li>Quality Assurance</li>
                <li>Analytical Testing</li>
            </ul>
        </div>
        
        <div style="flex: 1 1 200px;">
            <h3>Contact</h3>
            <p>123 Pharma Drive, Research Park<br>
            CA 94103, United States<br>
            Email: info@qrx-pharma.com<br>
            Phone: +1 (555) 123-4567</p>
        </div>
    </div>
    
    <div class="footer-bottom">
        <p>&copy; 2023 QRx Pharmaceuticals. All rights reserved.</p>
    </div>
</div>
""", unsafe_allow_html=True)