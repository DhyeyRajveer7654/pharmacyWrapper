import streamlit as st
import streamlit.components.v1 as components
from string import Template
from rdkit import Chem
from rdkit.Chem import Draw
import requests
import os

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
                <li>Documentation preparation for regulatory submissions</li>
                <li>Gap analysis against regional regulatory requirements</li>
                <li>Regulatory strategy development</li>
                <li>Compliance training for your team</li>
            </ul>
        </div>
        """, unsafe_allow_html=True)
        
        if st.button("Learn More About Regulatory Services", key="reg_services_btn"):
            change_page('regulatory')
    
    with col2:
        st.markdown("""
        <div class="card">
            <h3>Quality Assurance</h3>
            <p>Ensure pharmaceutical product quality with our comprehensive QA services:</p>
            <ul>
                <li>Method development and validation</li>
                <li>QA system design and implementation</li>
                <li>Advanced chemical structure analysis</li>
                <li>FTIR and analytical testing for raw materials</li>
            </ul>
        </div>
        """, unsafe_allow_html=True)
        
        if st.button("Explore Quality Services", key="qa_services_btn"):
            change_page('quality')
    
    # About us blurb
    st.markdown("""
    <div class="card" style="background: linear-gradient(135deg, #f0f9ff, #e0f2fe); text-align: center; padding: 2rem;">
        <h3>Why Choose QRx?</h3>
        <p style="margin-bottom: 1.5rem;">Our team of industry experts combines decades of experience with cutting-edge technology to deliver tailored solutions for your pharmaceutical quality and regulatory needs.</p>
        <p>We've helped over 200 pharmaceutical companies achieve and maintain compliance while optimizing their quality processes.</p>
    </div>
    """, unsafe_allow_html=True)
    
    if st.button("About Our Team", key="about_btn_home"):
        change_page('about')
    
    # Contact section
    st.markdown("""
    <div class="card" style="background: linear-gradient(135deg, #f0f9ff, #e0f2fe); text-align: center; padding: 2rem; margin-top: 2rem;">
        <h3>Ready to Elevate Your Pharmaceutical Quality & Compliance?</h3>
        <p style="margin-bottom: 1.5rem;">Contact our team of experts to discuss how we can support your specific needs.</p>
    </div>
    """, unsafe_allow_html=True)
    
    if st.button("Contact Us", key="contact_btn_home"):
        change_page('contact')

# Services Page
elif st.session_state.current_page == 'services':
    st.markdown('<h2 class="section-title">Our Services</h2>', unsafe_allow_html=True)
    
    # Introduction
    st.markdown("""
    <div class="card">
        <p>QRx offers a comprehensive suite of services designed to meet the complex needs of pharmaceutical companies. Our solutions combine industry expertise with innovative technology to ensure quality, compliance, and operational excellence.</p>
    </div>
    """, unsafe_allow_html=True)
    
    # Regulatory Services
    st.markdown('<h3 class="section-title">Regulatory Compliance</h3>', unsafe_allow_html=True)
    
    st.markdown("""
    <div class="card">
        <h3>Documentation and Submission Support</h3>
        <p>Our team provides expert assistance with preparing and reviewing regulatory documentation for submissions to health authorities worldwide, including:</p>
        <ul>
            <li>Common Technical Document (CTD) preparation</li>
            <li>Chemistry, Manufacturing, and Controls (CMC) documentation</li>
            <li>Response to regulatory queries</li>
            <li>Electronic submission preparation</li>
        </ul>
    </div>
    """, unsafe_allow_html=True)
    
    st.markdown("""
    <div class="card">
        <h3>Regulatory Strategy</h3>
        <p>Develop a clear regulatory pathway with our strategic services:</p>
        <ul>
            <li>Global regulatory roadmap development</li>
            <li>Regulatory gap analysis</li>
            <li>Product classification guidance</li>
            <li>Regulatory intelligence and monitoring</li>
        </ul>
    </div>
    """, unsafe_allow_html=True)
    
    st.markdown("""
    <div class="card">
        <h3>Compliance Management</h3>
        <p>Maintain ongoing compliance with regulatory requirements:</p>
        <ul>
            <li>GMP compliance assessments</li>
            <li>Regulatory inspection preparation</li>
            <li>Quality system development and review</li>
            <li>Compliance training programs</li>
        </ul>
    </div>
    """, unsafe_allow_html=True)
    
    if st.button("Explore Regulatory Solutions", key="reg_btn"):
        change_page('regulatory')
    
    # Quality Services
    st.markdown('<h3 class="section-title">Quality Assurance</h3>', unsafe_allow_html=True)
    
    st.markdown("""
    <div class="card">
        <h3>Product Quality Analysis</h3>
        <p>Comprehensive analytical services to ensure product quality:</p>
        <ul>
            <li>Chemical structure verification and visualization</li>
            <li>FTIR spectroscopy analysis</li>
            <li>Impurity profiling and identification</li>
            <li>Method development and validation</li>
        </ul>
    </div>
    """, unsafe_allow_html=True)
    
    st.markdown("""
    <div class="card">
        <h3>Quality System Development</h3>
        <p>Design and implement robust quality systems:</p>
        <ul>
            <li>Quality management system (QMS) development</li>
            <li>Standard operating procedures (SOPs) creation</li>
            <li>Quality risk management implementation</li>
            <li>Quality metrics and performance monitoring</li>
        </ul>
    </div>
    """, unsafe_allow_html=True)
    
    st.markdown("""
    <div class="card">
        <h3>Quality Audits</h3>
        <p>Comprehensive audit services to identify and address gaps:</p>
        <ul>
            <li>Internal quality audits</li>
            <li>Supplier and vendor qualification audits</li>
            <li>Mock regulatory inspections</li>
            <li>GMP/GDP compliance audits</li>
        </ul>
    </div>
    """, unsafe_allow_html=True)
    
    if st.button("Explore Quality Solutions", key="qa_btn"):
        change_page('quality')
    
    # Contact Call to Action
    st.markdown("""
    <div class="card" style="background: linear-gradient(135deg, #f0f9ff, #e0f2fe); text-align: center; padding: 2rem; margin-top: 2rem;">
        <h3>Need a Customized Solution?</h3>
        <p style="margin-bottom: 1.5rem;">Contact our specialists to discuss how we can tailor our services to your specific requirements.</p>
    </div>
    """, unsafe_allow_html=True)
    
    if st.button("Contact Our Team", key="contact_btn_services"):
        change_page('contact')

# Contact Page
elif st.session_state.current_page == 'contact':
    st.markdown('<h2 class="section-title">Contact Us</h2>', unsafe_allow_html=True)
    
    # Contact form
    st.markdown("""
    <div class="card">
        <p>We're here to help with your pharmaceutical quality and regulatory needs. Fill out the form below and one of our experts will get back to you within 24 hours.</p>
    </div>
    """, unsafe_allow_html=True)
    
    col1, col2 = st.columns(2)
    
    with col1:
        name = st.text_input("Full Name*")
        email = st.text_input("Email Address*")
        company = st.text_input("Company Name*")
        phone = st.text_input("Phone Number")
    
    with col2:
        service = st.selectbox(
            "Service of Interest*",
            ["Regulatory Compliance", "Quality Assurance", "Consulting Services", "Training Programs", "Other"]
        )
        
        if service == "Other":
            other_service = st.text_input("Please specify")
            
        urgency = st.selectbox(
            "Urgency",
            ["Standard (Response within 24 hours)", "Urgent (Response within 4 hours)", "Just exploring options"]
        )
        
        message = st.text_area("Message*", height=123, placeholder="Tell us about your specific needs...")
    
    terms_agree = st.checkbox("I agree to the processing of my personal data in accordance with the Privacy Policy")
    
    submit_col1, submit_col2, submit_col3 = st.columns([1, 2, 1])
    with submit_col2:
        if st.button("Submit", use_container_width=True):
            if not name or not email or not company or not message:
                st.error("Please fill in all required fields marked with *")
            elif not terms_agree:
                st.error("Please agree to the Privacy Policy to submit the form")
            else:
                st.success("Thank you for contacting us! We'll get back to you shortly.")
    
    # Contact information
    st.markdown("""
    <div class="card">
        <h3>Contact Information</h3>
        <div style="display: flex; flex-wrap: wrap; gap: 30px;">
            <div style="flex: 1; min-width: 200px;">
                <p><strong>Main Office:</strong><br>
                123 Pharma Boulevard, Suite 200<br>
                Boston, MA 02110</p>
            </div>
            <div style="flex: 1; min-width: 200px;">
                <p><strong>Phone:</strong><br>
                +1 (617) 555-0123</p>
                <p><strong>Email:</strong><br>
                info@qrxpharma.com</p>
            </div>
            <div style="flex: 1; min-width: 200px;">
                <p><strong>Hours of Operation:</strong><br>
                Monday - Friday: 8:00 AM - 6:00 PM EST<br>
                Saturday - Sunday: Closed</p>
            </div>
        </div>
    </div>
    """, unsafe_allow_html=True)
    
    # FAQs section
    st.markdown("""
    <div class="card">
        <h3>Frequently Asked Questions</h3>
        <p><strong>What is your typical response time for inquiries?</strong><br>
        We aim to respond to all inquiries within 24 hours during business days.</p>
        
        <p><strong>Do you offer services internationally?</strong><br>
        Yes, we provide services to pharmaceutical companies worldwide with expertise in global regulatory requirements.</p>
        
        <p><strong>Can you assist with urgent regulatory issues?</strong><br>
        Absolutely. We offer priority support for urgent regulatory matters with expedited response times.</p>
        
        <p><strong>Do you provide customized training for our team?</strong><br>
        Yes, we develop tailored training programs based on your specific needs and regulatory requirements.</p>
    </div>
    """, unsafe_allow_html=True)

# About Page
elif st.session_state.current_page == 'about':
    st.markdown('<h2 class="section-title">About QRx</h2>', unsafe_allow_html=True)
    
    # Company overview
    st.markdown("""
    <div class="card">
        <h3>Our Story</h3>
        <p>QRx was founded in 2010 by a team of pharmaceutical industry veterans who recognized the need for integrated quality and regulatory expertise in a rapidly evolving industry. Our mission is to help pharmaceutical companies navigate complex regulatory landscapes while maintaining the highest quality standards.</p>
        <p>With over a decade of experience, we've grown to become a trusted partner for pharmaceutical companies ranging from emerging startups to global enterprises. Our team combines deep industry knowledge with innovative technology solutions to deliver exceptional results for our clients.</p>
    </div>
    """, unsafe_allow_html=True)
    
    # Our approach
    st.markdown("""
    <div class="card">
        <h3>Our Approach</h3>
        <p>At QRx, we believe in a collaborative and tailored approach to addressing pharmaceutical quality and regulatory challenges. Our process includes:</p>
        <ol>
            <li><strong>Comprehensive Assessment:</strong> We begin by thoroughly understanding your specific needs, challenges, and objectives.</li>
            <li><strong>Strategic Planning:</strong> Our experts develop customized strategies that align with your business goals and regulatory requirements.</li>
            <li><strong>Implementation Support:</strong> We provide hands-on assistance to implement solutions effectively and efficiently.</li>
            <li><strong>Continuous Improvement:</strong> We monitor outcomes and make adjustments to ensure ongoing compliance and quality excellence.</li>
        </ol>
    </div>
    """, unsafe_allow_html=True)
    
    # Our team
    st.markdown("""
    <div class="card">
        <h3>Our Leadership Team</h3>
        <div style="display: flex; flex-wrap: wrap; gap: 20px;">
            <div style="flex: 1; min-width: 250px;">
                <h4>Dr. Sarah Johnson</h4>
                <p><em>Founder & CEO</em></p>
                <p>With over 25 years of experience in pharmaceutical quality and regulatory affairs, Dr. Johnson has led regulatory strategy for multiple successful FDA and EMA submissions. She holds a Ph.D. in Pharmaceutical Sciences and is a recognized industry thought leader.</p>
            </div>
            <div style="flex: 1; min-width: 250px;">
                <h4>Dr. Michael Chen</h4>
                <p><em>Chief Scientific Officer</em></p>
                <p>Dr. Chen leads our scientific and quality assurance services. With expertise in analytical chemistry and quality systems, he has helped clients resolve complex quality challenges for more than 20 years. He holds a Ph.D. in Analytical Chemistry.</p>
            </div>
            <div style="flex: 1; min-width: 250px;">
                <h4>Jennifer Williams</h4>
                <p><em>Head of Regulatory Affairs</em></p>
                <p>Jennifer brings 18 years of global regulatory experience across multiple therapeutic areas. She has successfully managed regulatory submissions in over 30 countries and specializes in navigating complex regulatory pathways.</p>
            </div>
        </div>
    </div>
    """, unsafe_allow_html=True)
    
    # Core values
    st.markdown("""
    <div class="card">
        <h3>Our Core Values</h3>
        <ul>
            <li><strong>Excellence:</strong> We are committed to delivering the highest quality services and solutions.</li>
            <li><strong>Integrity:</strong> We operate with honesty, transparency, and ethical standards in all that we do.</li>
            <li><strong>Innovation:</strong> We continuously seek new and improved ways to address challenges and enhance outcomes.</li>
            <li><strong>Collaboration:</strong> We work closely with our clients as true partners in their success.</li>
            <li><strong>Expertise:</strong> We maintain deep knowledge and stay current with evolving regulations and industry best practices.</li>
        </ul>
    </div>
    """, unsafe_allow_html=True)
    
    # Certifications and credentials
    st.markdown("""
    <div class="card">
        <h3>Certifications & Credentials</h3>
        <p>QRx maintains the highest professional standards and credentials in the industry:</p>
        <ul>
            <li>ISO 9001:2015 Certified</li>
            <li>Regulatory Affairs Professionals Society (RAPS) Corporate Member</li>
            <li>International Society for Pharmaceutical Engineering (ISPE) Member</li>
            <li>Parenteral Drug Association (PDA) Corporate Member</li>
            <li>International Council for Harmonisation of Technical Requirements for Pharmaceuticals for Human Use (ICH) Training Certified</li>
        </ul>
    </div>
    """, unsafe_allow_html=True)
    
    # Contact section
    st.markdown("""
    <div class="card" style="background: linear-gradient(135deg, #f0f9ff, #e0f2fe); text-align: center; padding: 2rem; margin-top: 2rem;">
        <h3>Partner with QRx for Your Pharmaceutical Quality & Regulatory Needs</h3>
        <p style="margin-bottom: 1.5rem;">Learn how our team can help you achieve and maintain compliance while optimizing your quality processes.</p>
    </div>
    """, unsafe_allow_html=True)
    
    if st.button("Contact Our Team", key="contact_btn_about"):
        change_page('contact')

# Regulatory Page
elif st.session_state.current_page == 'regulatory':
    st.markdown('<h2 class="section-title">Regulatory Compliance Services</h2>', unsafe_allow_html=True)
    
    # Introduction
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

# 📌 FORM PAGE
if st.session_state.page == "form":
    st.markdown('<div class="main-header"><h1>🧪 QAI Model AI-Powered Quality Assistance</h1><p> CREATED BY :- MEERA ACHARYA & RAJ PATEL</P><p>Enter details below to generate a comprehensive quality report</p></div>', unsafe_allow_html=True)

    # User Input Form in a card layout
    # st.markdown('<div class="card">', unsafe_allow_html=True)
    col1, col2 = st.columns(2)

    with col1:
        options["prodct_type"] = st.selectbox("🌎 Select Product Type", 
            ["API", "Tablets(com)", "Syrups", "Infusion", "Capsules", "Injectables","Other"])
        
    with col2:
        options["report_type"] = st.selectbox("🌎 Select Report Type", 
            ["Pathway", "List of license", "Detailed Information"])
        if options["report_type"] == "Detailed Information":
            options["resultsToCheck"] = st.text_area("🔍 Enter Your Results:", height=200, placeholder="Provide Licence you need info about here...", key="checkResults")
        options["regulatory"] = st.selectbox("🌎 Select Regulatory Authority", 
            ["CDSCO", "United States (FDA)", "European Union (EMA)","Brazil (ANVISA)", "Australia (TGA)"])
        
    # st.markdown('</div>', unsafe_allow_html=True)

    # Analysis Options in a separate card
    # st.markdown('<div class="card">', unsafe_allow_html=True)
   # Submit button with enhanced styling
    submit_button = st.button("🚀 Generate Report")
    if submit_button:
        if not all([options.get("prodct_type"), options.get("report_type"), options.get("regulatory")]):
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
    st.markdown(f"**💊 Product Type:** {st.session_state.prodct_type}",)
    st.markdown(f"**📦 Reoprt Type:** {st.session_state.report_type}")
    st.markdown(f"**⚡ Regulatory Authority:** {st.session_state.regulatory}")
    # st.markdown('</div>', unsafe_allow_ht ml=True)

    
    if st.session_state.api_response:
        components.html("<div class='table-container'>"+st.session_state.api_response+"</div>",height=800,width=1000,scrolling=True)
    else:
        st.warning("⚠️ No response received. Please try again.")

    
    if st.button("Contact Our Regulatory Team", key="contact_reg_btn"):
        change_page('contact')

# Quality Page
elif st.session_state.current_page == 'quality':
    st.markdown('<h2 class="section-title">Quality Assurance Services</h2>', unsafe_allow_html=True)
    
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
# Close main content wrapper
st.markdown('</div>', unsafe_allow_html=True)