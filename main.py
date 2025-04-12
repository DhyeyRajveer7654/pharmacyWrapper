import streamlit as st
import streamlit.components.v1 as components
from string import Template
import requests
import os

# Try to import RDKit with a fallback if it doesn't work
try:
    from rdkit import Chem
    from rdkit.Chem import Draw
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    print("RDKit not available. Chemical structure rendering will be disabled.")

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
    if not RDKIT_AVAILABLE:
        st.warning("RDKit is not available. Cannot render molecular structure.")
        return None
        
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
        text-align: center;
    }
    
    .footer a {
        color: #1e40af;
        text-decoration: none;
        transition: color 0.3s ease;
    }
    
    .footer a:hover {
        color: #3b82f6;
    }
    
    /* Button styling */
    .primary-button {
        background-color: #1e40af;
        color: white;
        padding: 0.7rem 1.5rem;
        border-radius: 5px;
        font-weight: 600;
        display: inline-block;
        transition: all 0.3s ease;
        text-decoration: none;
        border: none;
        cursor: pointer;
    }
    
    .primary-button:hover {
        background-color: #1e3a8a;
        transform: translateY(-2px);
    }
    
    /* Modal/popup styling */
    .popup-overlay {
        position: fixed;
        top: 0;
        left: 0;
        right: 0;
        bottom: 0;
        background-color: rgba(0, 0, 0, 0.5);
        z-index: 1000;
        display: flex;
        justify-content: center;
        align-items: center;
    }
    
    .popup-content {
        background-color: white;
        padding: 2rem;
        border-radius: 10px;
        max-width: 90%;
        width: 800px;
        max-height: 90vh;
        overflow-y: auto;
        position: relative;
    }
    
    .close-button {
        position: absolute;
        top: 1rem;
        right: 1rem;
        background: none;
        border: none;
        font-size: 1.5rem;
        cursor: pointer;
        color: #1e40af;
    }
</style>
""", unsafe_allow_html=True)

###############################################################################
# HOME PAGE
###############################################################################
def show_home_page():
    st.markdown('<div class="hero">', unsafe_allow_html=True)
    st.markdown('<h2>Your Partner in Pharmaceutical Quality & Regulatory Compliance</h2>', unsafe_allow_html=True)
    st.markdown('<p>QRx provides expert guidance and support for pharmaceutical companies navigating complex quality and regulatory requirements. Our specialized knowledge helps you bring safe, effective products to market faster.</p>', unsafe_allow_html=True)
    st.markdown('</div>', unsafe_allow_html=True)
    
    # Three columns for services highlights
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.markdown('<div class="card">', unsafe_allow_html=True)
        st.markdown('<h3>Quality Analysis</h3>', unsafe_allow_html=True)
        st.markdown('<p>Comprehensive quality testing and verification services to ensure your products meet all applicable standards.</p>', unsafe_allow_html=True)
        if st.button("Learn More", key="quality_learn_more"):
            st.session_state.current_page = 'quality'
            st.rerun()
        st.markdown('</div>', unsafe_allow_html=True)
        
    with col2:
        st.markdown('<div class="card">', unsafe_allow_html=True)
        st.markdown('<h3>Regulatory Support</h3>', unsafe_allow_html=True)
        st.markdown('<p>Navigate complex regulatory landscapes with our expert guidance on compliance requirements across multiple jurisdictions.</p>', unsafe_allow_html=True)
        if st.button("Learn More", key="regulatory_learn_more"):
            st.session_state.current_page = 'regulatory'
            st.rerun()
        st.markdown('</div>', unsafe_allow_html=True)
        
    with col3:
        st.markdown('<div class="card">', unsafe_allow_html=True)
        st.markdown('<h3>Consulting Services</h3>', unsafe_allow_html=True)
        st.markdown('<p>Strategic pharmaceutical consulting to optimize your development, manufacturing, and compliance processes.</p>', unsafe_allow_html=True)
        if st.button("Learn More", key="services_learn_more"):
            st.session_state.current_page = 'services'
            st.rerun()
        st.markdown('</div>', unsafe_allow_html=True)
    
    # About us section
    st.markdown('<h2 class="section-title">Why Choose QRx?</h2>', unsafe_allow_html=True)
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown('<div class="card">', unsafe_allow_html=True)
        st.markdown('<h3>Expert Team</h3>', unsafe_allow_html=True)
        st.markdown('<p>Our team of pharmaceutical scientists, regulatory specialists, and quality professionals brings decades of combined experience from leading global pharmaceutical companies.</p>', unsafe_allow_html=True)
        st.markdown('</div>', unsafe_allow_html=True)
        
    with col2:
        st.markdown('<div class="card">', unsafe_allow_html=True)
        st.markdown('<h3>Tailored Solutions</h3>', unsafe_allow_html=True)
        st.markdown('<p>We understand that every client has unique needs. Our services are customized to address your specific challenges and objectives.</p>', unsafe_allow_html=True)
        st.markdown('</div>', unsafe_allow_html=True)
    
    # CTA section
    st.markdown('<div style="text-align: center; margin: 3rem 0;">', unsafe_allow_html=True)
    st.markdown('<h2>Ready to streamline your pharmaceutical compliance?</h2>', unsafe_allow_html=True)
    if st.button("Contact Us Today", key="home_contact_button"):
        st.session_state.current_page = 'contact'
        st.rerun()
    st.markdown('</div>', unsafe_allow_html=True)
    
    # Footer
    st.markdown('<div class="footer">', unsafe_allow_html=True)
    st.markdown('<p>© 2025 QRx Pharmaceutical Consultants. All rights reserved.</p>', unsafe_allow_html=True)
    st.markdown('</div>', unsafe_allow_html=True)

###############################################################################
# SERVICES PAGE
###############################################################################
def show_services_page():
    st.markdown('<div class="hero">', unsafe_allow_html=True)
    st.markdown('<h2>Our Comprehensive Services</h2>', unsafe_allow_html=True)
    st.markdown('<p>QRx offers a wide range of pharmaceutical consulting services designed to help you navigate complex quality and regulatory requirements while optimizing your development and manufacturing processes.</p>', unsafe_allow_html=True)
    st.markdown('</div>', unsafe_allow_html=True)
    
    # Services overview
    st.markdown('<h2 class="section-title">Core Service Areas</h2>', unsafe_allow_html=True)
    
    # Create three columns for different service areas
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.markdown('<div class="card">', unsafe_allow_html=True)
        st.markdown('<h3>Research & Development</h3>', unsafe_allow_html=True)
        st.markdown('<ul>', unsafe_allow_html=True)
        st.markdown('<li>Formulation development support</li>', unsafe_allow_html=True)
        st.markdown('<li>Analytical method development</li>', unsafe_allow_html=True)
        st.markdown('<li>Stability studies design</li>', unsafe_allow_html=True)
        st.markdown('<li>Technical documentation</li>', unsafe_allow_html=True)
        st.markdown('</ul>', unsafe_allow_html=True)
        st.markdown('</div>', unsafe_allow_html=True)
        
    with col2:
        st.markdown('<div class="card">', unsafe_allow_html=True)
        st.markdown('<h3>Regulatory Affairs</h3>', unsafe_allow_html=True)
        st.markdown('<ul>', unsafe_allow_html=True)
        st.markdown('<li>Regulatory strategy development</li>', unsafe_allow_html=True)
        st.markdown('<li>Submission preparation and review</li>', unsafe_allow_html=True)
        st.markdown('<li>Regulatory intelligence</li>', unsafe_allow_html=True)
        st.markdown('<li>Post-approval maintenance</li>', unsafe_allow_html=True)
        st.markdown('</ul>', unsafe_allow_html=True)
        if st.button("See Regulatory Details", key="see_regulatory"):
            st.session_state.current_page = 'regulatory'
            st.rerun()
        st.markdown('</div>', unsafe_allow_html=True)
        
    with col3:
        st.markdown('<div class="card">', unsafe_allow_html=True)
        st.markdown('<h3>Quality Assurance</h3>', unsafe_allow_html=True)
        st.markdown('<ul>', unsafe_allow_html=True)
        st.markdown('<li>Quality management systems</li>', unsafe_allow_html=True)
        st.markdown('<li>Audit preparation and remediation</li>', unsafe_allow_html=True)
        st.markdown('<li>Quality risk management</li>', unsafe_allow_html=True)
        st.markdown('<li>Validation and qualification</li>', unsafe_allow_html=True)
        st.markdown('</ul>', unsafe_allow_html=True)
        if st.button("See Quality Details", key="see_quality"):
            st.session_state.current_page = 'quality'
            st.rerun()
        st.markdown('</div>', unsafe_allow_html=True)
    
    # Additional services section
    st.markdown('<h2 class="section-title">Specialized Services</h2>', unsafe_allow_html=True)
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown('<div class="card">', unsafe_allow_html=True)
        st.markdown('<h3>Training & Education</h3>', unsafe_allow_html=True)
        st.markdown('<p>We offer customized training programs for your staff on various aspects of pharmaceutical development, manufacturing, quality, and regulatory affairs.</p>', unsafe_allow_html=True)
        st.markdown('</div>', unsafe_allow_html=True)
        
    with col2:
        st.markdown('<div class="card">', unsafe_allow_html=True)
        st.markdown('<h3>Due Diligence</h3>', unsafe_allow_html=True)
        st.markdown('<p>Comprehensive evaluation of pharmaceutical products, facilities, and documentation for potential investments, mergers, or acquisitions.</p>', unsafe_allow_html=True)
        st.markdown('</div>', unsafe_allow_html=True)
    
    # CTA section
    st.markdown('<div style="text-align: center; margin: 3rem 0;">', unsafe_allow_html=True)
    st.markdown('<h2>Ready to optimize your pharmaceutical operations?</h2>', unsafe_allow_html=True)
    if st.button("Contact Us Today", key="services_contact_button"):
        st.session_state.current_page = 'contact'
        st.rerun()
    st.markdown('</div>', unsafe_allow_html=True)
    
    # Footer
    st.markdown('<div class="footer">', unsafe_allow_html=True)
    st.markdown('<p>© 2025 QRx Pharmaceutical Consultants. All rights reserved.</p>', unsafe_allow_html=True)
    st.markdown('</div>', unsafe_allow_html=True)

###############################################################################
# CONTACT PAGE
###############################################################################
def show_contact_page():
    st.markdown('<div class="hero">', unsafe_allow_html=True)
    st.markdown('<h2>Contact QRx</h2>', unsafe_allow_html=True)
    st.markdown('<p>We\'re here to answer your questions and discuss how we can support your pharmaceutical quality and regulatory compliance needs. Reach out to our team of experts today.</p>', unsafe_allow_html=True)
    st.markdown('</div>', unsafe_allow_html=True)
    
    # Contact form and information
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown('<div class="card">', unsafe_allow_html=True)
        st.markdown('<h3>Send Us a Message</h3>', unsafe_allow_html=True)
        
        # Contact form
        name = st.text_input("Name")
        email = st.text_input("Email")
        company = st.text_input("Company")
        subject = st.selectbox("Subject", ["General Inquiry", "Quality Services", "Regulatory Support", "Partnership Opportunities", "Other"])
        message = st.text_area("Message", height=150)
        
        if st.button("Submit"):
            st.success("Your message has been sent! Our team will get back to you within 24 hours.")
            
        st.markdown('</div>', unsafe_allow_html=True)
        
    with col2:
        st.markdown('<div class="card">', unsafe_allow_html=True)
        st.markdown('<h3>Contact Information</h3>', unsafe_allow_html=True)
        
        st.markdown('<p><strong>Email:</strong> redoxYlifecare@gmail.com</p>', unsafe_allow_html=True)
        st.markdown('<p><strong>Phone:</strong> +918849122744, +919723449306 </p>', unsafe_allow_html=True)
        st.markdown('<p><strong>Address:</strong> 2,3 medicare complex, old housing road <br> Surendranagar, 363001 <br> India </p>', unsafe_allow_html=True)
        
        st.markdown('<h4 style="margin-top: 2rem;">Office Hours</h4>', unsafe_allow_html=True)
        st.markdown('<p>Monday - Saturday: 10:00 AM - 8:00 PM IST</p>', unsafe_allow_html=True)
        
        st.markdown('</div>', unsafe_allow_html=True)
    
    # Footer
    st.markdown('<div class="footer">', unsafe_allow_html=True)
    st.markdown('<p>© 2025 QRx Pharmaceutical Consultants. All rights reserved.</p>', unsafe_allow_html=True)
    st.markdown('</div>', unsafe_allow_html=True)

###############################################################################
# ABOUT PAGE
###############################################################################
def show_about_page():
    st.markdown('<div class="hero">', unsafe_allow_html=True)
    st.markdown('<h2>About QRx</h2>', unsafe_allow_html=True)
    st.markdown('<p>QRx is a premier pharmaceutical consulting firm specializing in quality assurance and regulatory compliance. With decades of combined experience, our team of experts provides comprehensive solutions to help pharmaceutical companies navigate complex regulatory landscapes and maintain the highest quality standards.</p>', unsafe_allow_html=True)
    st.markdown('</div>', unsafe_allow_html=True)
    
    # Company overview
    st.markdown('<h2 class="section-title">Our Story</h2>', unsafe_allow_html=True)
    
    st.markdown('<div class="card">', unsafe_allow_html=True)
    st.markdown('<p>Founded in 2015 by a team of pharmaceutical industry veterans, QRx was established with a mission to simplify quality and regulatory compliance for pharmaceutical companies of all sizes. What began as a small consulting practice has grown into a trusted partner for dozens of companies across the pharmaceutical industry.</p>', unsafe_allow_html=True)
    st.markdown('<p>Our founders recognized that many pharmaceutical companies struggle with navigating the complex and ever-changing regulatory landscape while maintaining rigorous quality standards. Drawing on their experience from leading global pharmaceutical firms, they created QRx to provide expert guidance and practical solutions to these challenges.</p>', unsafe_allow_html=True)
    st.markdown('<p>Today, QRx serves clients ranging from emerging biotech startups to established pharmaceutical manufacturers, helping them bring safe and effective products to market efficiently while maintaining compliance with global regulations.</p>', unsafe_allow_html=True)
    st.markdown('</div>', unsafe_allow_html=True)
    
    # Core values section
    st.markdown('<h2 class="section-title">Our Core Values</h2>', unsafe_allow_html=True)
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.markdown('<div class="card">', unsafe_allow_html=True)
        st.markdown('<h3>Excellence</h3>', unsafe_allow_html=True)
        st.markdown('<p>We are committed to delivering the highest quality service and exceeding client expectations in everything we do.</p>', unsafe_allow_html=True)
        st.markdown('</div>', unsafe_allow_html=True)
        
    with col2:
        st.markdown('<div class="card">', unsafe_allow_html=True)
        st.markdown('<h3>Integrity</h3>', unsafe_allow_html=True)
        st.markdown('<p>We operate with honesty, transparency, and strong ethical principles in all our client relationships and business practices.</p>', unsafe_allow_html=True)
        st.markdown('</div>', unsafe_allow_html=True)
        
    with col3:
        st.markdown('<div class="card">', unsafe_allow_html=True)
        st.markdown('<h3>Innovation</h3>', unsafe_allow_html=True)
        st.markdown('<p>We continuously seek creative and efficient solutions to complex pharmaceutical quality and regulatory challenges.</p>', unsafe_allow_html=True)
        st.markdown('</div>', unsafe_allow_html=True)
    
    # Team section
    st.markdown('<h2 class="section-title">Our Leadership Team</h2>', unsafe_allow_html=True)
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown('<div class="card">', unsafe_allow_html=True)
        st.markdown('<h3>Dr. Sarah Johnson</h3>', unsafe_allow_html=True)
        st.markdown('<p><strong>Founder & CEO</strong></p>', unsafe_allow_html=True)
        st.markdown('<p>Dr. Johnson has over 20 years of experience in pharmaceutical development and regulatory affairs, previously serving as VP of Regulatory Affairs at a global pharmaceutical company. She holds a Ph.D. in Pharmaceutical Sciences and an MBA.</p>', unsafe_allow_html=True)
        st.markdown('</div>', unsafe_allow_html=True)
        
    with col2:
        st.markdown('<div class="card">', unsafe_allow_html=True)
        st.markdown('<h3>Dr. Michael Chen</h3>', unsafe_allow_html=True)
        st.markdown('<p><strong>Chief Scientific Officer</strong></p>', unsafe_allow_html=True)
        st.markdown('<p>With extensive experience in pharmaceutical quality systems and compliance, Dr. Chen leads our quality consulting practice. He previously held senior positions at FDA and major pharmaceutical manufacturers. He holds a Ph.D. in Chemistry.</p>', unsafe_allow_html=True)
        st.markdown('</div>', unsafe_allow_html=True)
    
    # Footer
    st.markdown('<div class="footer">', unsafe_allow_html=True)
    st.markdown('<p>© 2025 QRx Pharmaceutical Consultants. All rights reserved.</p>', unsafe_allow_html=True)
    st.markdown('</div>', unsafe_allow_html=True)

###############################################################################
# REGULATORY PAGE
###############################################################################
def show_regulatory_page():
    st.markdown('<div class="hero">', unsafe_allow_html=True)
    st.markdown('<h2>Regulatory Compliance Services</h2>', unsafe_allow_html=True)
    st.markdown('<p>QRx offers comprehensive regulatory compliance services to help pharmaceutical companies navigate complex regulatory landscapes across global markets. Our team of regulatory experts provides strategic guidance and practical support throughout the product lifecycle.</p>', unsafe_allow_html=True)
    st.markdown('</div>', unsafe_allow_html=True)
    
    # Regulatory services overview
    st.markdown('<h2 class="section-title">Our Regulatory Services</h2>', unsafe_allow_html=True)
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown('<div class="card">', unsafe_allow_html=True)
        st.markdown('<h3>Regulatory Strategy</h3>', unsafe_allow_html=True)
        st.markdown('<p>We develop comprehensive regulatory strategies tailored to your specific products and target markets. Our approach helps you navigate complex regulatory pathways efficiently, saving time and resources while maximizing chances for approval.</p>', unsafe_allow_html=True)
        st.markdown('<ul>', unsafe_allow_html=True)
        st.markdown('<li>Global regulatory pathway assessment</li>', unsafe_allow_html=True)
        st.markdown('<li>Regulatory risk evaluation</li>', unsafe_allow_html=True)
        st.markdown('<li>Strategic planning for submissions</li>', unsafe_allow_html=True)
        st.markdown('<li>Regulatory intelligence and monitoring</li>', unsafe_allow_html=True)
        st.markdown('</ul>', unsafe_allow_html=True)
        st.markdown('</div>', unsafe_allow_html=True)
        
    with col2:
        st.markdown('<div class="card">', unsafe_allow_html=True)
        st.markdown('<h3>Submission Support</h3>', unsafe_allow_html=True)
        st.markdown('<p>Our experts assist with preparation, review, and management of regulatory submissions across multiple jurisdictions, ensuring compliance with all requirements and standards.</p>', unsafe_allow_html=True)
        st.markdown('<ul>', unsafe_allow_html=True)
        st.markdown('<li>IND/CTA preparation and submissions</li>', unsafe_allow_html=True)
        st.markdown('<li>NDA/MAA preparation and submissions</li>', unsafe_allow_html=True)
        st.markdown('<li>DMF preparation and maintenance</li>', unsafe_allow_html=True)
        st.markdown('<li>Responses to regulatory authority queries</li>', unsafe_allow_html=True)
        st.markdown('</ul>', unsafe_allow_html=True)
        st.markdown('</div>', unsafe_allow_html=True)
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown('<div class="card">', unsafe_allow_html=True)
        st.markdown('<h3>Regulatory Compliance</h3>', unsafe_allow_html=True)
        st.markdown('<p>We help ensure ongoing compliance with regulatory requirements throughout the product lifecycle, from development through post-marketing phases.</p>', unsafe_allow_html=True)
        st.markdown('<ul>', unsafe_allow_html=True)
        st.markdown('<li>GMP compliance assessment and remediation</li>', unsafe_allow_html=True)
        st.markdown('<li>Regulatory inspection preparation</li>', unsafe_allow_html=True)
        st.markdown('<li>Post-approval change management</li>', unsafe_allow_html=True)
        st.markdown('<li>Pharmacovigilance and safety reporting</li>', unsafe_allow_html=True)
        st.markdown('</ul>', unsafe_allow_html=True)
        st.markdown('</div>', unsafe_allow_html=True)
        
    with col2:
        st.markdown('<div class="card">', unsafe_allow_html=True)
        st.markdown('<h3>Global Regulatory Support</h3>', unsafe_allow_html=True)
        st.markdown('<p>Our global regulatory expertise covers major and emerging markets worldwide, helping you navigate varied requirements across different regions.</p>', unsafe_allow_html=True)
        st.markdown('<ul>', unsafe_allow_html=True)
        st.markdown('<li>US FDA regulations compliance</li>', unsafe_allow_html=True)
        st.markdown('<li>European Medicines Agency (EMA) requirements</li>', unsafe_allow_html=True)
        st.markdown('<li>Health Canada submissions</li>', unsafe_allow_html=True)
        st.markdown('<li>Emerging markets regulatory strategy</li>', unsafe_allow_html=True)
        st.markdown('</ul>', unsafe_allow_html=True)
        st.markdown('</div>', unsafe_allow_html=True)
    
    # Regulatory expertise highlights
    st.markdown('<h2 class="section-title">Regulatory Areas of Expertise</h2>', unsafe_allow_html=True)
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.markdown('<div class="card">', unsafe_allow_html=True)
        st.markdown('<h3>Small Molecules</h3>', unsafe_allow_html=True)
        st.markdown('<p>Comprehensive regulatory support for small molecule drug products, including generics, new chemical entities, and complex formulations.</p>', unsafe_allow_html=True)
        st.markdown('</div>', unsafe_allow_html=True)
        
    with col2:
        st.markdown('<div class="card">', unsafe_allow_html=True)
        st.markdown('<h3>Biologics</h3>', unsafe_allow_html=True)
        st.markdown('<p>Specialized regulatory guidance for biological products, biosimilars, and advanced therapy medicinal products (ATMPs).</p>', unsafe_allow_html=True)
        st.markdown('</div>', unsafe_allow_html=True)
        
    with col3:
        st.markdown('<div class="card">', unsafe_allow_html=True)
        st.markdown('<h3>Combination Products</h3>', unsafe_allow_html=True)
        st.markdown('<p>Expert navigation of the complex regulatory landscape for drug-device combination products.</p>', unsafe_allow_html=True)
        st.markdown('</div>', unsafe_allow_html=True)
    
    # CTA section
    st.markdown('<div style="text-align: center; margin: 3rem 0;">', unsafe_allow_html=True)
    st.markdown('<h2>Need regulatory support for your pharmaceutical products?</h2>', unsafe_allow_html=True)
    if st.button("Contact Our Regulatory Experts", key="regulatory_contact_button"):
        st.session_state.current_page = 'contact'
        st.rerun()
    st.markdown('</div>', unsafe_allow_html=True)
    
    # Footer
    st.markdown('<div class="footer">', unsafe_allow_html=True)
    st.markdown('<p>© 2025 QRx Pharmaceutical Consultants. All rights reserved.</p>', unsafe_allow_html=True)
    st.markdown('</div>', unsafe_allow_html=True)

###############################################################################
# QUALITY PAGE
###############################################################################
def show_quality_page():
    st.markdown('<div class="hero">', unsafe_allow_html=True)
    st.markdown('<h2>Pharmaceutical Quality Services</h2>', unsafe_allow_html=True)
    st.markdown('<p>QRx delivers comprehensive quality analysis and assurance services to help pharmaceutical companies meet global quality standards. Our experienced team ensures your products are manufactured with the highest level of quality and consistency.</p>', unsafe_allow_html=True)
    st.markdown('</div>', unsafe_allow_html=True)
    
    # Example quality analysis form - ONLY SHOWN ON QUALITY PAGE
    st.markdown('<h2 class="section-title">Quality Analysis Capabilities</h2>', unsafe_allow_html=True)
    
    st.markdown('<div class="card">', unsafe_allow_html=True)
    st.markdown('<h3>Comprehensive Product Quality Assessment</h3>', unsafe_allow_html=True)
    
    # Sample product selection form (simplified for demo)
    st.markdown('<p>Explore our quality analysis capabilities with this interactive demo:</p>', unsafe_allow_html=True)
    
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
        
    st.markdown('</div>', unsafe_allow_html=True)
    
    # Quality services overview - moved after the analysis form
    st.markdown('<h2 class="section-title">Our Quality Services</h2>', unsafe_allow_html=True)
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown('<div class="card">', unsafe_allow_html=True)
        st.markdown('<h3>Quality Management Systems</h3>', unsafe_allow_html=True)
        st.markdown('<p>We help design, implement, and improve quality management systems that ensure compliance with global regulations while optimizing operational efficiency.</p>', unsafe_allow_html=True)
        st.markdown('<ul>', unsafe_allow_html=True)
        st.markdown('<li>QMS development and implementation</li>', unsafe_allow_html=True)
        st.markdown('<li>Gap analysis and remediation</li>', unsafe_allow_html=True)
        st.markdown('<li>SOP development and review</li>', unsafe_allow_html=True)
        st.markdown('<li>Quality metrics and performance monitoring</li>', unsafe_allow_html=True)
        st.markdown('</ul>', unsafe_allow_html=True)
        st.markdown('</div>', unsafe_allow_html=True)
        
    with col2:
        st.markdown('<div class="card">', unsafe_allow_html=True)
        st.markdown('<h3>Quality Assurance</h3>', unsafe_allow_html=True)
        st.markdown('<p>Our quality assurance services help maintain ongoing compliance and product quality throughout the manufacturing process.</p>', unsafe_allow_html=True)
        st.markdown('<ul>', unsafe_allow_html=True)
        st.markdown('<li>Batch record review</li>', unsafe_allow_html=True)
        st.markdown('<li>Change control management</li>', unsafe_allow_html=True)
        st.markdown('<li>Deviation investigation and CAPA implementation</li>', unsafe_allow_html=True)
        st.markdown('<li>Release testing oversight</li>', unsafe_allow_html=True)
        st.markdown('</ul>', unsafe_allow_html=True)
        st.markdown('</div>', unsafe_allow_html=True)
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown('<div class="card">', unsafe_allow_html=True)
        st.markdown('<h3>Quality Control</h3>', unsafe_allow_html=True)
        st.markdown('<p>We provide expert guidance on testing methods, specifications, and laboratory operations to ensure accurate and reliable product quality data.</p>', unsafe_allow_html=True)
        st.markdown('<ul>', unsafe_allow_html=True)
        st.markdown('<li>Analytical method validation</li>', unsafe_allow_html=True)
        st.markdown('<li>Laboratory compliance (GLP/GMP)</li>', unsafe_allow_html=True)
        st.markdown('<li>Stability testing programs</li>', unsafe_allow_html=True)
        st.markdown('<li>Specification development and justification</li>', unsafe_allow_html=True)
        st.markdown('</ul>', unsafe_allow_html=True)
        st.markdown('</div>', unsafe_allow_html=True)
        
    with col2:
        st.markdown('<div class="card">', unsafe_allow_html=True)
        st.markdown('<h3>Auditing & Inspection Readiness</h3>', unsafe_allow_html=True)
        st.markdown('<p>We help prepare your facilities, systems, and personnel for regulatory inspections and customer audits.</p>', unsafe_allow_html=True)
        st.markdown('<ul>', unsafe_allow_html=True)
        st.markdown('<li>Mock FDA/EMA inspections</li>', unsafe_allow_html=True)
        st.markdown('<li>Audit response support</li>', unsafe_allow_html=True)
        st.markdown('<li>Supplier qualification audits</li>', unsafe_allow_html=True)
        st.markdown('<li>GMP audit programs</li>', unsafe_allow_html=True)
        st.markdown('</ul>', unsafe_allow_html=True)
        st.markdown('</div>', unsafe_allow_html=True)
    
    # CTA section
    st.markdown('<div style="text-align: center; margin: 3rem 0;">', unsafe_allow_html=True)
    st.markdown('<h2>Ready to elevate your quality assurance program?</h2>', unsafe_allow_html=True)
    if st.button("Contact Our Quality Experts", key="quality_contact_button"):
        st.session_state.current_page = 'contact'
        st.rerun()
    st.markdown('</div>', unsafe_allow_html=True)
    
    # Footer
    st.markdown('<div class="footer">', unsafe_allow_html=True)
    st.markdown('<p>© 2025 QRx Pharmaceutical Consultants. All rights reserved.</p>', unsafe_allow_html=True)
    st.markdown('</div>', unsafe_allow_html=True)

###############################################################################
# MAIN PAGE RENDERING
###############################################################################

# Display the appropriate page based on session state
if st.session_state.current_page == 'home':
    show_home_page()
elif st.session_state.current_page == 'services':
    show_services_page()
elif st.session_state.current_page == 'contact':
    show_contact_page()
elif st.session_state.current_page == 'about':
    show_about_page()
elif st.session_state.current_page == 'regulatory':
    show_regulatory_page()
elif st.session_state.current_page == 'quality':
    show_quality_page()