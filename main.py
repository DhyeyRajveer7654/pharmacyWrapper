import streamlit as st
import streamlit.components.v1 as components
from string import Template
import os
import requests
from rdkit import Chem 
from rdkit.Chem import Draw

# Set flag for RDKit availability
RDKIT_AVAILABLE = True

###############################################################################
# UTILITY FUNCTIONS
###############################################################################

def chatWithGpt(prompt):
    """Simulate a GPT API call with a deterministic response"""
    # In a real implementation, this would call an API
    # For demo purposes, return a formatted HTML response
    html_response = f"""
    <h3>Quality Analysis Report</h3>
    <table style="width:100%; border-collapse: collapse;">
        <tr style="background-color: #1e40af; color: white;">
            <th style="padding: 8px; text-align: left;">Parameter</th>
            <th style="padding: 8px; text-align: left;">Specification</th>
            <th style="padding: 8px; text-align: left;">Result</th>
        </tr>
        <tr>
            <td style="padding: 8px; border-bottom: 1px solid #ddd;">Description</td>
            <td style="padding: 8px; border-bottom: 1px solid #ddd;">White to off-white crystalline powder</td>
            <td style="padding: 8px; border-bottom: 1px solid #ddd;">Complies</td>
        </tr>
        <tr>
            <td style="padding: 8px; border-bottom: 1px solid #ddd;">Identification</td>
            <td style="padding: 8px; border-bottom: 1px solid #ddd;">IR spectrum matches reference standard</td>
            <td style="padding: 8px; border-bottom: 1px solid #ddd;">Complies</td>
        </tr>
        <tr>
            <td style="padding: 8px; border-bottom: 1px solid #ddd;">Assay</td>
            <td style="padding: 8px; border-bottom: 1px solid #ddd;">98.0% - 102.0%</td>
            <td style="padding: 8px; border-bottom: 1px solid #ddd;">99.7%</td>
        </tr>
        <tr>
            <td style="padding: 8px; border-bottom: 1px solid #ddd;">Dissolution</td>
            <td style="padding: 8px; border-bottom: 1px solid #ddd;">NLT 80% in 30 minutes</td>
            <td style="padding: 8px; border-bottom: 1px solid #ddd;">92% in 30 minutes</td>
        </tr>
        <tr>
            <td style="padding: 8px; border-bottom: 1px solid #ddd;">Related Substances</td>
            <td style="padding: 8px; border-bottom: 1px solid #ddd;">Any individual impurity: NMT 0.5%</td>
            <td style="padding: 8px; border-bottom: 1px solid #ddd;">Maximum individual impurity: 0.3%</td>
        </tr>
        <tr>
            <td style="padding: 8px; border-bottom: 1px solid #ddd;">Water Content</td>
            <td style="padding: 8px; border-bottom: 1px solid #ddd;">NMT 0.5%</td>
            <td style="padding: 8px; border-bottom: 1px solid #ddd;">0.3%</td>
        </tr>
    </table>
    """
    return html_response

def get_ftir_from_gpt(product_name):
    """Generate FTIR analysis data for a product"""
    ftir_html = f"""
    <h3>FTIR Analysis for {product_name}</h3>
    <p>Key peaks identified:</p>
    <ul>
        <li>3400-3200 cm<sup>-1</sup>: O-H stretching</li>
        <li>2960-2850 cm<sup>-1</sup>: C-H stretching</li>
        <li>1700-1680 cm<sup>-1</sup>: C=O stretching</li>
        <li>1600-1450 cm<sup>-1</sup>: Aromatic ring vibrations</li>
        <li>1300-1000 cm<sup>-1</sup>: C-O stretching</li>
    </ul>
    <p>All characteristic peaks match the reference standard for {product_name}.</p>
    """
    return ftir_html

# Template for drug structure prompt
STRUCTURE_PROMPT = Template("What is the SMILES notation for $product_name?")

def getPromptForOptions(options):
    """Generate a prompt based on the provided options"""
    prompt = f"Generate a quality analysis report for {options.get('product_name', 'Unknown')} "
    prompt += f"with strength {options.get('powerOfDrug', 'Unknown')} "
    prompt += f"in quantity {options.get('quanOfMed', 'Unknown')} "
    prompt += f"according to {options.get('jurisdiction', 'Unknown')} standards."
    
    if options.get('typeOfInfo') == "METHOD OF PREPARATION":
        prompt += " Focus on method of preparation."
    elif options.get('typeOfInfo') == "CHARACTARIZATION/EVALUATION":
        prompt += " Focus on characterization and evaluation."
    elif options.get('typeOfInfo') == "Both of above":
        prompt += " Include both method of preparation and characterization/evaluation."
    elif options.get('typeOfInfo') == "CHECK RESULTS":
        prompt += f" Evaluate the following results: {options.get('resultsToCheck', '')}."
    
    return prompt

def get_cid_from_name(drug_name):
    """Get compound ID from PubChem by compound name"""
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
    """Get SMILES notation for a compound from PubChem"""
    product_code_from_pubchem = ""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{product_name}/property/CanonicalSMILES/JSON"
    response = requests.get(url)
    if response.status_code == 200:
        try:
            smiles = response.json()["PropertyTable"]["Properties"][0]["CanonicalSMILES"]
            product_code_from_pubchem = smiles
        except (KeyError, IndexError):
            product_code_from_pubchem = "NO DRUG FOUND"
    else:
        product_code_from_pubchem = "NO DRUG FOUND"
    if product_code_from_pubchem == "NO DRUG FOUND":
        return ""
    else:
        return product_code_from_pubchem
        
def showStructure(product_name):
    """Generate molecular structure image from compound name"""
    product_code = ""
    product_code_from_pubchem = get_pubchem_product_code(product_name)
    if product_code_from_pubchem == "":
        product_code_prompt = STRUCTURE_PROMPT.substitute(product_name=product_name)
        print("Prompt is: " + product_code_prompt)
        product_code = chatWithGpt(product_code_prompt)
        if product_code == "NO DRUG FOUND":
            return ""
    else:
        product_code = product_code_from_pubchem

    print("product code is: " + product_code)
    print("product code from pubchem: " + product_code_from_pubchem)
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

###############################################################################
# MAIN APPLICATION
###############################################################################

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

# Custom CSS for the navigation buttons
st.markdown("""
<style>
.nav-button {
    background-color: transparent !important;
    color: #1e40af !important;
    border: none !important;
    font-weight: 600 !important;
    padding: 8px 16px !important;
    text-transform: uppercase !important;
    letter-spacing: 0.5px !important;
    border-radius: 4px !important;
    transition: all 0.3s ease !important;
}

.nav-button:hover {
    background-color: #e0f2fe !important;
    color: #1e3a8a !important;
    transform: translateY(-2px) !important;
}

/* Make sure the active page button looks different */
.stButton button[data-testid="BaseButton"] {
    width: 100%;
}
</style>
""", unsafe_allow_html=True)

st.markdown('<div style="background-color: white; padding: 1rem 0; border-bottom: 1px solid #e0e0e0; margin-bottom: 1rem;">', unsafe_allow_html=True)
col1, col2, col3, col4, col5, col6, col7 = st.columns(7)
with col1:
    st.markdown('<h1 style="color:#1e40af; margin:0; padding:0; font-size:2.5rem; font-weight:700;">QRx</h1>', unsafe_allow_html=True)
with col2:
    if st.button("HOME", key="nav_home", use_container_width=True, type="secondary", help="Go to home page"):
        st.session_state.current_page = 'home'
        st.rerun()
with col3:
    if st.button("SERVICES", key="nav_services", use_container_width=True, type="secondary", help="View our services"):
        st.session_state.current_page = 'services'
        st.rerun()
with col4:
    if st.button("CONTACT", key="nav_contact", use_container_width=True, type="secondary", help="Contact us"):
        st.session_state.current_page = 'contact'
        st.rerun()
with col5:
    if st.button("ABOUT", key="nav_about", use_container_width=True, type="secondary", help="About QRx"):
        st.session_state.current_page = 'about'
        st.rerun()
with col6:
    if st.button("REGULATORY", key="nav_regulatory", use_container_width=True, type="secondary", help="Regulatory Info"):
        st.session_state.current_page = 'regulatory'
        st.rerun()
with col7:
    if st.button("QUALITY", key="nav_quality", use_container_width=True, type="secondary", help="Quality Info"):
        st.session_state.current_page = 'quality'
        st.rerun()
st.markdown('</div>', unsafe_allow_html=True)

# Define navigation functions
def close_popup():
    st.session_state.show_popup = False
    st.rerun()

def change_page(page):
    st.session_state.current_page = page
    st.session_state.show_popup = False
    st.rerun()

def handle_js_close():
    # This will be called after the JS close action to update session state
    st.session_state.show_popup = False
    st.rerun()

# Custom CSS styling
st.markdown("""
<style>
    /* Global styles */
    body {
        font-family: 'Segoe UI', sans-serif;
        color: #333;
        background-color: #ffffff;
    }
    
    /* Header styling */
    .header-container {
        position: fixed;
        top: 0;
        left: 0;
        right: 0;
        z-index: 9999;
        background-color: white;
        box-shadow: 0 2px 10px rgba(0, 0, 0, 0.1);
        padding: 0.8rem 2rem;
        height: 70px;
    }
    
    .header {
        display: flex;
        justify-content: space-between;
        align-items: center;
        max-width: 1400px;
        margin: 0 auto;
        height: 100%;
    }
    
    .logo {
        color: #1e40af;
        font-size: 2.5rem;
        font-weight: 700;
        margin: 0;
        letter-spacing: -1px;
    }
    
    .nav-links {
        display: flex;
        gap: 1.5rem;
        align-items: center;
        height: 100%;
    }
    
    .nav-link {
        color: #1e40af;
        text-decoration: none;
        font-weight: 600;
        font-size: 0.9rem;
        letter-spacing: 0.5px;
        padding: 0.5rem 1rem;
        border-radius: 4px;
        transition: all 0.3s ease;
        display: inline-block;
    }
    
    .nav-link:hover {
        background-color: #e0f2fe;
        color: #1e3a8a;
    }
    
    /* Add margin to content to prevent overlap with fixed header */
    .main-content {
        margin-top: 90px; /* Increased to ensure content doesn't hide under header */
        padding-top: 1rem;
    }
    
    /* Hero section */
    .hero {
        background: linear-gradient(135deg, #e0f2fe, #bfdbfe);
        color: #1e40af;
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
        background-color: #f0f9ff;
        color: #1e3a8a;
        padding: 2rem;
        border-radius: 10px 10px 0 0;
        margin-top: 3rem;
    }
    
    .footer h3 {
        color: #1e40af;
        margin-bottom: 1rem;
    }
    
    .footer-bottom {
        text-align: center;
        margin-top: 2rem;
        padding-top: 1rem;
        border-top: 1px solid #bfdbfe;
        color: #1e40af;
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
    
    /* Enhanced Form Elements from OLD AI.py */
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
    
    /* Results Table */
    .table-container {
        background: white;
        border-radius: 10px;
        padding: 1rem;
        box-shadow: 0 4px 6px rgba(0,0,0,0.1);
        border-collapse: collapse;
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
    
    /* Main header styling */
    .main-header {
        text-align: center;
        margin-bottom: 2rem;
    }

    .main-header h1 {
        color: #0ea5e9;
        font-size: 2rem;
        margin-bottom: 0.5rem;
    }

    .main-header p {
        color: #64748b;
        font-size: 1.1rem;
    }

    /* Results styling */
    .result-section {
        margin-top: 2rem;
    }

    .result-header {
        color: #0ea5e9;
        font-size: 1.5rem;
        margin-bottom: 1rem;
        border-bottom: 2px solid #e2e8f0;
        padding-bottom: 0.5rem;
    }
</style>
""", unsafe_allow_html=True)

# Also handle direct page changes via URLs or buttons
if 'direct_nav' in st.session_state and st.session_state.direct_nav:
    page = st.session_state.direct_nav
    st.session_state.direct_nav = None
    change_page(page)

# Main content wrapper
st.markdown('<div class="main-content">', unsafe_allow_html=True)

###############################################################################
# PAGE CONTENT BASED ON CURRENT PAGE
###############################################################################

# Home Page
if st.session_state.current_page == 'home':
    # Hero section
    st.markdown("""
    <div class="hero">
        <h2>Pharmaceutical Quality & Regulatory Excellence</h2>
        <p>QRx provides comprehensive quality assurance and regulatory compliance solutions for the pharmaceutical industry. 
        With our expert team and advanced AI-powered tools, we ensure your products meet the highest standards at every stage.</p>
    </div>
    """, unsafe_allow_html=True)
    
    # Services section
    st.markdown('<h2 class="section-title">Our Services</h2>', unsafe_allow_html=True)
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("""
        <div class="card">
            <h3>Regulatory Compliance</h3>
            <p>Navigate complex regulatory requirements with our comprehensive compliance services:</p>
            <ul>
                <li>Regulatory strategy development</li>
                <li>Product registration support</li>
                <li>Regulatory submissions management</li>
                <li>Post-approval regulatory maintenance</li>
            </ul>
        </div>
        """, unsafe_allow_html=True)
        
    with col2:
        st.markdown("""
        <div class="card">
            <h3>Quality Assurance</h3>
            <p>Ensure product safety, efficacy, and compliance with our quality assurance services:</p>
            <ul>
                <li>Quality management system development</li>
                <li>GMP compliance audits</li>
                <li>Quality control testing</li>
                <li>CAPA system implementation</li>
            </ul>
        </div>
        """, unsafe_allow_html=True)
    
    # CTA section
    st.markdown("""
    <div style="text-align: center; margin: 3rem 0;">
        <h2>Ready to elevate your pharmaceutical quality and compliance?</h2>
    </div>
    """, unsafe_allow_html=True)
    
    col1, col2 = st.columns(2)
    with col1:
        if st.button("Contact Us", key="home_contact"):
            st.session_state.current_page = 'contact'
            st.rerun()
    with col2:
        if st.button("Learn More About Our Services", key="home_services"):
            st.session_state.current_page = 'services'
            st.rerun()
    
    # Footer
    st.markdown("""
    <div class="footer">
        <p>© 2025 QRx Pharmaceutical Consultants. All rights reserved.</p>
    </div>
    """, unsafe_allow_html=True)

# Services Page
elif st.session_state.current_page == 'services':
    st.markdown('<h2 class="section-title">Our Services</h2>', unsafe_allow_html=True)
    
    st.markdown("""
    <div class="hero">
        <h2>Comprehensive Pharmaceutical Services</h2>
        <p>QRx offers a wide range of services designed to support pharmaceutical companies throughout the product lifecycle, from development to post-marketing activities.</p>
    </div>
    """, unsafe_allow_html=True)
    
    # Regulatory Services
    st.markdown('<h3 class="section-title">Regulatory Services</h3>', unsafe_allow_html=True)
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("""
        <div class="card">
            <h3>Regulatory Strategy</h3>
            <p>Develop comprehensive strategies for product registration and lifecycle management:</p>
            <ul>
                <li>Global regulatory pathway assessment</li>
                <li>Strategic planning for market access</li>
                <li>Regulatory gap analysis</li>
                <li>Product development roadmaps</li>
            </ul>
        </div>
        """, unsafe_allow_html=True)
        
    with col2:
        st.markdown("""
        <div class="card">
            <h3>Regulatory Submissions</h3>
            <p>Complete preparation and management of regulatory dossiers:</p>
            <ul>
                <li>CTD/eCTD compilation</li>
                <li>CMC documentation</li>
                <li>Clinical and non-clinical summaries</li>
                <li>Agency meeting preparation</li>
            </ul>
        </div>
        """, unsafe_allow_html=True)
    
    # Quality Services
    st.markdown('<h3 class="section-title">Quality Services</h3>', unsafe_allow_html=True)
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("""
        <div class="card">
            <h3>Quality Management Systems</h3>
            <p>Design and implement robust quality systems for GMP compliance:</p>
            <ul>
                <li>QMS development and implementation</li>
                <li>SOP creation and management</li>
                <li>Quality metrics implementation</li>
                <li>Training program development</li>
            </ul>
        </div>
        """, unsafe_allow_html=True)
        
    with col2:
        st.markdown("""
        <div class="card">
            <h3>Analytical Services</h3>
            <p>Comprehensive quality control and product characterization:</p>
            <ul>
                <li>Method development and validation</li>
                <li>Product characterization</li>
                <li>Stability testing</li>
                <li>Reference standard qualification</li>
            </ul>
        </div>
        """, unsafe_allow_html=True)
    
    # CTA section
    st.markdown("""
    <div style="text-align: center; margin: 3rem 0;">
        <h2>Interested in our services?</h2>
    </div>
    """, unsafe_allow_html=True)
    
    col1, col2 = st.columns(2)
    with col1:
        if st.button("Contact Our Regulatory Team", key="services_reg_contact"):
            st.session_state.current_page = 'regulatory'
            st.rerun()
    with col2:
        if st.button("Contact Our Quality Team", key="services_quality_contact"):
            st.session_state.current_page = 'quality'
            st.rerun()
            
    # Footer
    st.markdown("""
    <div class="footer">
        <p>© 2025 QRx Pharmaceutical Consultants. All rights reserved.</p>
    </div>
    """, unsafe_allow_html=True)

# Contact Page
elif st.session_state.current_page == 'contact':
    st.markdown('<h2 class="section-title">Contact Us</h2>', unsafe_allow_html=True)
    
    st.markdown("""
    <div class="hero">
        <h2>Get in Touch</h2>
        <p>Have questions about our services? Need expert assistance with your pharmaceutical quality and regulatory challenges? Contact our team and get the help you need.</p>
    </div>
    """, unsafe_allow_html=True)
    
    # Contact Form
    st.markdown("""
    <div class="card">
        <h3>Contact Form</h3>
        <p>Fill out the form below and one of our experts will get back to you within 24 hours.</p>
    </div>
    """, unsafe_allow_html=True)
    
    col1, col2 = st.columns(2)
    
    with col1:
        name = st.text_input("Name")
        email = st.text_input("Email")
        phone = st.text_input("Phone Number")
        
    with col2:
        company = st.text_input("Company")
        topic = st.selectbox("Topic", ["Regulatory Services", "Quality Services", "General Inquiry", "Other"])
        
    message = st.text_area("Message", height=150)
    
    if st.button("Submit", key="contact_submit"):
        # In a real implementation, this would send the form data
        st.success("Thank you for contacting us! We'll get back to you shortly.")
    
    # Office Locations 
    st.markdown("""
    <div class="card">
        <h3>Our Offices</h3>
    </div>
    """, unsafe_allow_html=True)
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("""
        <div class="card">
            <h3>Main Office</h3>
            <p>123 Pharma Boulevard<br>
            Boston, MA 02110<br>
            United States</p>
            <p>Phone: +1 (555) 123-4567<br>
            Email: info@qrxpharma.com</p>
        </div>
        """, unsafe_allow_html=True)
        
    with col2:
        st.markdown("""
        <div class="card">
            <h3>European Office</h3>
            <p>45 Regent Square<br>
            London, EC1V 2NP<br>
            United Kingdom</p>
            <p>Phone: +44 (0) 20 7123 4567<br>
            Email: europe@qrxpharma.com</p>
        </div>
        """, unsafe_allow_html=True)
    
    # Footer
    st.markdown("""
    <div class="footer">
        <p>© 2025 QRx Pharmaceutical Consultants. All rights reserved.</p>
    </div>
    """, unsafe_allow_html=True)

# About Page
elif st.session_state.current_page == 'about':
    st.markdown('<h2 class="section-title">About QRx</h2>', unsafe_allow_html=True)
    
    st.markdown("""
    <div class="hero">
        <h2>Our Story</h2>
        <p>QRx was founded by pharmaceutical industry veterans with a mission to simplify quality and regulatory compliance for pharmaceutical companies of all sizes.</p>
    </div>
    """, unsafe_allow_html=True)
    
    # Company Profile
    st.markdown("""
    <div class="card">
        <h3>Who We Are</h3>
        <p>QRx is a specialized pharmaceutical consulting firm focused on quality and regulatory excellence. Founded in 2015, we've helped over 200 companies across 25 countries navigate complex regulatory pathways and implement robust quality systems.</p>
        <p>Our team consists of former regulatory agency officials, quality directors, and pharmaceutical development experts who bring decades of real-world experience to every client project.</p>
    </div>
    """, unsafe_allow_html=True)
    
    # Mission and Values
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("""
        <div class="card">
            <h3>Our Mission</h3>
            <p>To advance global health by helping pharmaceutical companies navigate complex quality and regulatory challenges efficiently and effectively.</p>
            <p>We believe in making compliance accessible to companies of all sizes, from startups to multinational corporations, through practical solutions and technology-enabled services.</p>
        </div>
        """, unsafe_allow_html=True)
        
    with col2:
        st.markdown("""
        <div class="card">
            <h3>Our Values</h3>
            <ul>
                <li><strong>Excellence:</strong> We maintain the highest standards in all our work</li>
                <li><strong>Integrity:</strong> We operate with honesty and transparency</li>
                <li><strong>Innovation:</strong> We embrace new technologies and approaches</li>
                <li><strong>Partnership:</strong> We work as an extension of our clients' teams</li>
                <li><strong>Impact:</strong> We measure success by our clients' outcomes</li>
            </ul>
        </div>
        """, unsafe_allow_html=True)
    
    # Team Section
    st.markdown('<h3 class="section-title">Our Leadership</h3>', unsafe_allow_html=True)
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.markdown("""
        <div class="card">
            <h3>Dr. Sarah Johnson</h3>
            <p><em>CEO & Founder</em></p>
            <p>Former FDA reviewer with 20+ years of experience in pharmaceutical regulatory affairs. Ph.D. in Pharmaceutical Sciences.</p>
        </div>
        """, unsafe_allow_html=True)
        
    with col2:
        st.markdown("""
        <div class="card">
            <h3>Dr. Michael Chen</h3>
            <p><em>Chief Quality Officer</em></p>
            <p>Previously Quality Director at a global pharmaceutical company. Expert in QMS implementation and remediation.</p>
        </div>
        """, unsafe_allow_html=True)
        
    with col3:
        st.markdown("""
        <div class="card">
            <h3>Dr. Elena Rodriguez</h3>
            <p><em>Head of Regulatory Strategy</em></p>
            <p>Former EMA scientific advisor with expertise in global regulatory strategy and submissions.</p>
        </div>
        """, unsafe_allow_html=True)
    
    # Footer
    st.markdown("""
    <div class="footer">
        <p>© 2025 QRx Pharmaceutical Consultants. All rights reserved.</p>
    </div>
    """, unsafe_allow_html=True)

# Regulatory Page - ONLY SHOWS REGULATORY CONTENT
elif st.session_state.current_page == 'regulatory':
    st.markdown('<h2 class="section-title">Regulatory Compliance Services</h2>', unsafe_allow_html=True)
    
    # Introduction
    st.markdown("""
    <div class="hero">
        <h2>Navigate Complex Regulatory Frameworks</h2>
        <p>QRx offers comprehensive regulatory compliance services to help pharmaceutical companies navigate complex regulatory landscapes across global markets. Our team of regulatory experts provides strategic guidance and practical support throughout the product lifecycle.</p>
    </div>
    """, unsafe_allow_html=True)
    
    # Regulatory services
    st.markdown('<h3 class="section-title">Our Regulatory Services</h3>', unsafe_allow_html=True)
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("""
        <div class="card">
            <h3>Regulatory Strategy</h3>
            <p>We develop comprehensive regulatory strategies tailored to your specific products and target markets. Our approach helps you navigate complex regulatory pathways efficiently, saving time and resources while maximizing chances for approval.</p>
            <ul>
                <li>Global regulatory pathway assessment</li>
                <li>Regulatory risk evaluation</li>
                <li>Strategic planning for submissions</li>
                <li>Regulatory intelligence and monitoring</li>
            </ul>
        </div>
        """, unsafe_allow_html=True)
        
    with col2:
        st.markdown("""
        <div class="card">
            <h3>Submission Support</h3>
            <p>Our experts assist with preparation, review, and management of regulatory submissions across multiple jurisdictions, ensuring compliance with all requirements and standards.</p>
            <ul>
                <li>IND/CTA preparation and submissions</li>
                <li>NDA/MAA preparation and submissions</li>
                <li>DMF preparation and maintenance</li>
                <li>Responses to regulatory authority queries</li>
            </ul>
        </div>
        """, unsafe_allow_html=True)
    
    # REGULATORY FORM - ONLY SHOWN ON REGULATORY PAGE
    st.markdown('<div class="main-header"><h1>📋 Regulatory Compliance Assistant</h1><p>Enter details below to generate comprehensive regulatory information</p></div>', unsafe_allow_html=True)

    # User Input Form
    options = dict()
    
    col1, col2 = st.columns(2)

    with col1:
        options["prodct_type"] = st.selectbox("🌎 Select Product Type", 
            ["API", "Tablets(com)", "Syrups", "Infusion", "Capsules", "Injectables","Other"])
            
    with col2:
        options["report_type"] = st.selectbox("📝 Report Type", 
            ["Regulatory Overview", "Detailed Information", "License Requirements", "Documentation Needed"])
            
    if options["report_type"] == "Detailed Information":
        options["resultsToCheck"] = st.text_area("🔍 Enter Your Results:", height=200, placeholder="Provide License you need info about here...", key="checkResults")
    options["regulatory"] = st.selectbox("🌎 Select Regulatory Authority", 
        ["CDSCO", "United States (FDA)", "European Union (EMA)","Brazil (ANVISA)", "Australia (TGA)"])
    
    # Submit button
    submit_button = st.button("🚀 Generate Regulatory Report")
    if submit_button:
        if not all([options.get("prodct_type"), options.get("report_type"), options.get("regulatory")]):
            st.error("⚠️ Please fill in all required fields!")
        else:
            # Display the results - simulating API response
            st.markdown('<div class="main-header"><h1>📑 Regulatory Compliance Report</h1></div>', unsafe_allow_html=True)
            
            st.markdown("### 📋 Analysis Details")
            st.markdown(f"**💊 Product Type:** {options.get('prodct_type')}",)
            st.markdown(f"**📝 Report Type:** {options.get('report_type')}")
            if options.get('resultsToCheck'):
                st.markdown(f"**🔍 Results Analyzed:** Yes")
            st.markdown(f"**🌎 Regulatory Authority:** {options.get('regulatory')}")
            
            st.markdown('<div class="result-header">🔬 Regulatory Analysis</div>', unsafe_allow_html=True)
            
            # Simulate response
            api_response = """
            <div style="background-color: white; padding: 1.5rem; border-radius: 0.5rem; box-shadow: 0 4px 6px rgba(0,0,0,0.1);">
                <h3>Regulatory Requirements Summary</h3>
                <p>For <strong>{}</strong> in <strong>{}</strong>:</p>
                <ul>
                    <li><strong>Marketing Authorization:</strong> Required - Standard submission pathway</li>
                    <li><strong>GMP Certification:</strong> Required - On-site inspection typically conducted</li>
                    <li><strong>Product Testing:</strong> Required - Local laboratory testing may be necessary</li>
                    <li><strong>Labeling Requirements:</strong> Package insert and patient information must be in local language</li>
                    <li><strong>Post-Market Surveillance:</strong> Active reporting of adverse events required</li>
                </ul>
                <h4>Key Documentation Required:</h4>
                <ol>
                    <li>Common Technical Document (CTD) or equivalent</li>
                    <li>Quality Control Methods and Validation Reports</li>
                    <li>Stability Data (Long-term and Accelerated)</li>
                    <li>Process Validation Reports</li>
                    <li>GMP Certificate</li>
                </ol>
                <p><strong>Estimated Timeline:</strong> 12-18 months for full approval process</p>
                <p><strong>Key Challenges:</strong> Local language requirements, potential for additional testing, inspection scheduling can cause delays</p>
            </div>
            """.format(options.get('prodct_type'), options.get('regulatory'))
            st.markdown(api_response, unsafe_allow_html=True)

    # CTA section
    st.markdown("""
    <div style="text-align: center; margin: 3rem 0;">
        <h2>Need personalized regulatory support?</h2>
    </div>
    """, unsafe_allow_html=True)
    
    if st.button("Contact Our Regulatory Team", key="contact_reg_btn"):
        change_page('contact')
        
    # Footer
    st.markdown("""
    <div class="footer">
        <p>© 2025 QRx Pharmaceutical Consultants. All rights reserved.</p>
    </div>
    """, unsafe_allow_html=True)

# Quality Page - ONLY SHOWS QUALITY CONTENT
elif st.session_state.current_page == 'quality':
    st.markdown('<h2 class="section-title">Pharmaceutical Quality Services</h2>', unsafe_allow_html=True)
    
    # Introduction
    st.markdown("""
    <div class="hero">
        <h2>Excellence in Pharmaceutical Quality</h2>
        <p>QRx delivers comprehensive quality analysis and assurance services to help pharmaceutical companies meet global quality standards. Our experienced team ensures your products are manufactured with the highest level of quality and consistency.</p>
    </div>
    """, unsafe_allow_html=True)
    
    # QUALITY FORM - ONLY SHOWN ON QUALITY PAGE
    st.markdown('<div class="main-header"><h1>🧪 Quality Analysis Assistant</h1><p>Enter details below to generate a comprehensive quality report</p></div>', unsafe_allow_html=True)

    # Quality form
    col1, col2 = st.columns(2)
    
    with col1:
        product_name = st.text_input("Enter product name (e.g., Acetaminophen)", "Acetaminophen")
        jurisdiction = st.selectbox("Select regulatory jurisdiction", ["US FDA", "EMA (Europe)", "Health Canada", "WHO"])
        
    with col2:
        strength = st.text_input("Strength (e.g., 500 mg)", "500 mg")
        quantity = st.text_input("Quantity (e.g., 100 tablets)", "100 tablets")
        
    analysis_type = st.radio("Select analysis type", ["CHARACTARIZATION/EVALUATION", "METHOD OF PREPARATION", "Both of above"])
    
    if st.button("Generate Quality Analysis", key="gen_quality_analysis"):
        options = {
            'product_name': product_name,
            'powerOfDrug': strength,
            'quanOfMed': quantity,
            'jurisdiction': jurisdiction,
            'typeOfInfo': analysis_type
        }
        
        # Generate and display analysis report
        st.markdown('### Quality Analysis Results', unsafe_allow_html=True)
        
        # Structure display (if RDKit is available)
        col1, col2 = st.columns(2)
        
        with col1:
            if RDKIT_AVAILABLE:
                structure_img = showStructure(product_name)
                if structure_img is not None:
                    st.image(structure_img, caption=f"Molecular structure of {product_name}")
                else:
                    st.warning(f"Could not generate structure for {product_name}")
            else:
                st.warning("Molecular structure display is currently unavailable")
        
        with col2:
            # FTIR analysis
            ftir_analysis = get_ftir_from_gpt(product_name)
            st.markdown(ftir_analysis, unsafe_allow_html=True)
        
        # Quality report
        quality_report = chatWithGpt(getPromptForOptions(options))
        st.markdown(quality_report, unsafe_allow_html=True)
    
    # Quality services overview - after the analysis form
    st.markdown('<h3 class="section-title">Our Quality Services</h3>', unsafe_allow_html=True)
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("""
        <div class="card">
            <h3>Quality Management Systems</h3>
            <p>We help design, implement, and improve quality management systems that ensure compliance with global regulations while optimizing operational efficiency.</p>
            <ul>
                <li>QMS development and implementation</li>
                <li>Gap analysis and remediation</li>
                <li>SOP development and review</li>
                <li>Quality metrics and performance monitoring</li>
            </ul>
        </div>
        """, unsafe_allow_html=True)
        
    with col2:
        st.markdown("""
        <div class="card">
            <h3>Quality Assurance</h3>
            <p>Our quality assurance services help maintain ongoing compliance and product quality throughout the manufacturing process.</p>
            <ul>
                <li>Batch record review</li>
                <li>Change control management</li>
                <li>Deviation investigation and CAPA implementation</li>
                <li>Release testing oversight</li>
            </ul>
        </div>
        """, unsafe_allow_html=True)
    
    # CTA section
    st.markdown("""
    <div style="text-align: center; margin: 3rem 0;">
        <h2>Ready to elevate your quality assurance program?</h2>
    </div>
    """, unsafe_allow_html=True)
    
    if st.button("Contact Our Quality Experts", key="quality_contact_button"):
        st.session_state.current_page = 'contact'
        st.rerun()
        
    # Footer
    st.markdown("""
    <div class="footer">
        <p>© 2025 QRx Pharmaceutical Consultants. All rights reserved.</p>
    </div>
    """, unsafe_allow_html=True)

# Close the main content div
st.markdown('</div>', unsafe_allow_html=True)