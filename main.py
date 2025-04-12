import streamlit as st
import streamlit.components.v1 as components
import prompts
import chat_with_gpt
from string import Template
from rdkit import Chem
from rdkit.Chem import Draw
import requests
import os
import base64

def display_pdf(file_path):
    with open(file_path, "rb") as f:
        base64_pdf = base64.b64encode(f.read()).decode("utf-8")
    pdf_display = f'''
        <iframe src="data:application/pdf;base64,{base64_pdf}" width="100%" height="600" type="application/pdf">
        </iframe>
    '''
    return pdf_display
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
col1, col2, col3, col4, col5, col6 = st.columns(6)
with col1:
    st.markdown('<h1 style="color:#1e40af; margin:0; padding:0; font-size:2.5rem; font-weight:700;">QRx</h1>', unsafe_allow_html=True)
with col2:
    if st.button("HOME", key="nav_home", use_container_width=True, type="secondary", help="Go to home page"):
        st.session_state.current_page = 'home'
        st.rerun()
with col3:
    if st.button("CONTACT", key="nav_contact", use_container_width=True, type="secondary", help="Contact us"):
        st.session_state.current_page = 'contact'
        st.rerun()
with col4:
    if st.button("ABOUT", key="nav_about", use_container_width=True, type="secondary", help="About QRx"):
        st.session_state.current_page = 'about'
        st.rerun()
with col5:
    if st.button("REGULATORY", key="nav_regulatory", use_container_width=True, type="secondary", help="Regulatory Info"):
        st.session_state.current_page = 'regulatory'
        st.rerun()
with col6:
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
        height: 10px;
    }
    
    .header {
        display: flex;
        center-content: space-between;
        align-items: center;
        max-width: 1400px;
        margin: 0 auto;
        height: 50%;
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
        margin-top: 10px; /* Increased to ensure content doesn't hide under header */
        padding-top: 0.5rem;
    }
    
    /* Hero section */
    .hero {
        background: linear-gradient(135deg, #e0f2fe, #bfdbfe);
        color: #1e40af;
        padding: 1rem 1rem;
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
        padding: 1rem;
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
        background-color: #fffff;
        color: #1e3a8a;
        padding: 1rem;
        border-radius: 5px 5px 0 0;
        margin-top: 1rem;
    }
    
    .footer h3 {
        color: #ffffff;
        margin-bottom: 1rem;
    }
    
    .footer-bottom {
        text-align: center;
        margin-top: 1rem;
        padding-top: 1rem;
        border-top: 1px solid #ffffff;
        color: #ffffff;
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
        color: #000000;
        font-size: 3rem;
        margin-bottom: 0.5rem;
    }

    .main-header p {
        color: #000000;
        font-size: 1.1rem;
    }

    /* Results styling */
    .result-section {
        margin-top: 2rem;
    }

    .result-header {
        color: #000000;
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
        <p>QRx AI provides comprehensive quality assurance and regulatory compliance solutions for the pharmaceutical students. 
        With our models and advanced AI-powered tools, we ensure your products meet the highest standards at every stage.</p>
    </div>
    """, unsafe_allow_html=True)
    
    # Services section
    st.markdown('<h2 class="section-title">Our Services</h2>', unsafe_allow_html=True)
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("""
        <div class="card">
            <h3>QRx AI powered Regulatory Complaince</h3>
            <p>Navigate complex regulatory requirements with our comprehensive compliance model:</p>
            <ul>
                <li>Pathway for beginners</li>
                <li>Liciense check list</li>
                <li>Detailed Information on Particular forms</li>
            </ul>
        </div>
        """, unsafe_allow_html=True)
        
    with col2:
        st.markdown("""
        <div class="card">
            <h3>QAI AI Powered Quality Assurance</h3>
            <p>Ensure product safety, efficacy, and compliance with our quality assurance services:</p>
            <ul>
                <li>Pharmacopieal Complaince</li>
                <li>Deatiled Method of Preparation</li>
                <li>All the Evaluation Checklist</li>
                <li>Check your results with pharmacopiea</li>
                <li>FTIR Graphs of Anti-Cancer Drugs (IP)</li>
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
        if st.button("Learn More About us", key="home_About_us"):
            st.session_state.current_page = 'About us'
            st.rerun()
    col1, col_gap, col2 = st.columns([3, 1, 6])
    with col1:
        st.markdown("## 📄 ISO Certification")
        if os.path.exists("iso_preview.jpg"):
            st.image("iso_preview.jpg", caption="ISO Certificate Preview", width=400)
        else:
            st.warning("ISO preview image not found.")
        with open("iso_certificate.pdf", "rb") as f:
            st.download_button("📥 Download Full ISO Certificate", f, file_name="iso_certificate.pdf")
        
    with col2:
        st.markdown("## 🏆 QAI Model Award")
    
    # Optional: check if image exists
        if os.path.exists("ncip_award.jpg"):
            st.image("ncip_award.jpg", caption="Award Ceremony – Nirma University, 2025", use_column_width=True)
        else:
            st.warning("Award image not found. Please upload 'ncip_award.jpg'.")

    # Professional award description
    st.markdown("""
    **QAI MODEL** proudly secured **2nd Prize** at the **National Conference of Institute of Pharmacy (NCIP 2025)**  
    hosted by **Nirma University**. This recognition reflects the innovation, impact, and future potential of QRx AI in transforming pharmaceutical education and compliance.
    """)
    # Footer
    st.markdown("""
    <div class="footer">
        <p>© 2025 QRx AI. All rights reserved. Built with passion by pharmacy minds.</p>
    </div>
    """, unsafe_allow_html=True)


# Contact Page
elif st.session_state.current_page == 'contact':
    st.markdown("""
    <div class="hero">
        <h2>Get in Touch</h2>
        <p>Have questions about our services? Need expert assistance with your pharmaceutical quality and regulatory challenges? Contact our team and get the help you need.</p>
    </div>
    """, unsafe_allow_html=True)
    
    col1, col2 = st.columns(2)

    with col1:
        st.markdown("""
        <div class="card">
            <h3>Main Office</h3>
            <h4>Redoxy Lifecare</h4>
            <p>2,3 Medicare Complex,<br>
            Old Housing Road,<br>
            Surendranagar, 363001 <br>
            Gujarat, India</p>
            <p>Phone: +91-8849122744, +91-9723449306<br>
            Email: redoxylifecare@gmail.com</p>
        </div>
        """, unsafe_allow_html=True)

    with col2:
        st.image('LOGO.png', width=200)

        # Footer
    st.markdown("""
        <div class="footer">
            <p>© 2025 QRx Pharmaceutical Consultants. All rights reserved.</p>
        </div>
        """, unsafe_allow_html=True)

# About Page
elif st.session_state.current_page == 'about':
   
    st.markdown("""
    <div class="hero">
        <h2>Our Story</h2>
        <p>QRx AI was founded in 2024 by Meera Acharya, a passionate pharmacy student at A.P.M.C. College of Pharmacy and Research, with a dream to blend artificial intelligence and pharmaceutical science. 
        Her journey began in the 6th semester of her B.Pharm at C.U. Shah College of Pharmacy and Research, where she fell in love with AI. 
        Meera never hesitates to take up an opportunity in this field — her dedication led her to win 2nd prize for the same AI model at NCIP 2025.</p>
    </div>
    """, unsafe_allow_html=True)

    # Company Profile
    st.markdown("""
    <div class="card">
        <h3>Who We Are</h3>
        <p>QRx AI is a startup focused on empowering pharmacy students by simplifying regulatory compliance and enabling them to perform dosage form manufacturing as per pharmacopoeial standards in laboratory settings.</p>
        <p>Registered under Redocy Lifecare, QRx AI is co-founded by Raj H. Patel, with Meera Acharya as the founder. Our tools provide easy-to-use solutions and educational resources tailored for pharma students.</p>
    </div>
    """, unsafe_allow_html=True)

    # Mission and Values
    col1, col2 = st.columns(2)

    with col1:
        st.markdown("""
        <div class="card">
            <h3>Our Mission</h3>
            <p>To create something unique in the field of pharmacy by integrating artificial intelligence into pharmaceutical processes.</p>
            <p>We aim to bridge the gap between AI and the pharma industry while educating students and professionals alike.</p>
        </div>
        """, unsafe_allow_html=True)

    with col2:
        st.markdown("""
        <div class="card">
            <h3>Our Values</h3>
            <ul>
                <li><strong>Passion:</strong> Driven by curiosity and love for AI</li>
                <li><strong>Education:</strong> Focused on making learning easy and practical</li>
                <li><strong>Innovation:</strong> Building smart, student-focused tools</li>
                <li><strong>Teamwork:</strong> Collaborating across disciplines</li>
                <li><strong>Impact:</strong> Making real-world lab experiences more accessible</li>
            </ul>
        </div>
        """, unsafe_allow_html=True)

    # Team Section
    st.markdown('<h3 class="section-title">Our Leadership</h3>', unsafe_allow_html=True)

    col1, col2, col3 = st.columns(3)

    with col1:
        st.image("meera.png", width=120)
        st.markdown("""
            <h3>Meera Rahul Acharya</h3>
            <p><em>Founder & Tech Lead</em></p>
            <p>Innovator and visionary behind QRx AI. Obsessed with AI and passionate about transforming pharmacy education. She ideates and designs all project implementations.</p>
        """, unsafe_allow_html=True)

    with col2:
        st.image("raj.png", width=120)
        st.markdown("""
            <h3>Raj H Patel</h3>
            <p><em>Co-Founder & Tech Lead</em></p>
            <p>Owner of Redocy Lifecare and co-founder of QRx AI. Raj handles all the coding and back-end development, turning Meera's ideas into working AI solutions.</p>
        """, unsafe_allow_html=True)

    with col3:
        st.image("dhyey.png", width=120)
        st.markdown("""
            <h3>Dhyey Rajveer</h3>
            <p><em>Model Development Support</em></p>
            <p>Key contributor in the creation and fine-tuning of the AI models. His technical insights have been instrumental in QRx AI’s early success.</p>
        """, unsafe_allow_html=True)

        # Footer
        st.markdown("""
        <div class="footer">
            <p>© 2025 QRx AI. All rights reserved. Built with passion by pharmacy minds.</p>
        </div>
        """, unsafe_allow_html=True)

# Regulatory Page - ONLY SHOWS REGULATORY CONTENT
elif st.session_state.current_page == 'regulatory':
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
        st.markdown('<div class="main-header"><h1>🧪 QRx AI-Powered Regulatory Complaince</h1><p> CREATED BY :- MEERA ACHARYA & RAJ PATEL</P><p>Enter details below to generate a comprehensive quality report</p></div>', unsafe_allow_html=True)

        # User Input Form in a card layout
        # st.markdown('<div class="card">', unsafe_allow_html=True)
        col1, col2 = st.columns(2)

        with col1:
            options["prodct_type"] = st.selectbox("💊 Select Product Type", 
                ["Active Pharmaceutical Ingredient (API)", "Tablets(regular)", "Syrups", "Infusion", "Capsules", "Injectables","Other"])
            
        with col2:
            options["report_type"] = st.selectbox("📝 Select Report Type", 
                ["Pathway", "List of license", "Detailed Information"])
            if options["report_type"] == "Detailed Information":
                options["resultsToCheck"] = st.text_area("📝 Enter Your Results:", height=200, placeholder="Provide Licence you need info about here...", key="checkResults")
            options["regulatory"] = st.selectbox("🌎 Select Regulatory Authority", 
                ["CDSCO", "United States (FDA)", "European Union (EMA)","Brazil (ANVISA)", "Australia (TGA)"])
            
        # st.markdown('</div>', unsafe_allow_html=True)

        # Analysis Options in a separate card
        # st.markdown('<div class="card">', unsafe_allow_html=True)
    # Submit button with enhanced styling
        submit_button = st.button("🚀 Generate Regulatory Report")
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
        <p>© 2025 QRx AI. All rights reserved. Built with passion by pharmacy minds.</p>
    </div>
    """, unsafe_allow_html=True)

# Quality Page - ONLY SHOWS QUALITY CONTENT
elif st.session_state.current_page == 'quality':
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

    
    if st.button("Contact Our Quality Experts", key="quality_contact_button"):
        st.session_state.current_page = 'contact'
        st.rerun()
        
    # Footer
    st.markdown("""
    <div class="footer">
        <p>© 2025 QRx AI. All rights reserved. Built with passion by pharmacy minds.</p>
    </div>
    """, unsafe_allow_html=True)

# Close the main content div
st.markdown('</div>', unsafe_allow_html=True)