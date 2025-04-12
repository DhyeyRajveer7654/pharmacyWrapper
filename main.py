import streamlit as st
import streamlit.components.v1 as components
import prompts
import chat_with_gpt
from string import Template
from rdkit import Chem
from rdkit.Chem import Draw
import requests
import os

# Set flag for RDKit availability
RDKIT_AVAILABLE = True

###############################################################################
# UTILITY FUNCTIONS
###############################################################################

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
    col1, = st.columns(1)
    
    with col1:
        st.markdown("""
        <div class="card">
            <h3>Main Office</h3>
            <h4>Redoxy lifecare<h4>
            <p>2,3 medicare complex,<br>
            old housing road,<br>
            Surendranagar, 363001 <br>
            Gujarat, India</p>
            <p>Phone: +91-8849122744, +91-9723449306<br>
            Email: redoxylifecare@gmail.com</p>
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
    st.markdown('<h2 class="section-title">About Us</h2>', unsafe_allow_html=True)
    
    
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
        st.markdown("""
        <div class="card" style="text-align: center;">
            <img src="meera.png" alt="Meera Acharya" style="width:120px; height:120px; object-fit:cover; border-radius:50%; margin-bottom: 10px;">
            <h3>Meera Acharya</h3>
            <p><em>Founder & Tech Lead</em></p>
            <p>Innovator and visionary behind QRx AI. Obsessed with AI and passionate about transforming pharmacy education. She ideates and designs all project implementations.</p>
        </div>
        """, unsafe_allow_html=True)

    with col2:
        st.markdown("""
        <div class="card" style="text-align: center;">
            <img src="raj.png" alt="Raj H. Patel" style="width:120px; height:120px; object-fit:cover; border-radius:50%; margin-bottom: 10px;">
            <h3>Raj H. Patel</h3>
            <p><em>Co-Founder & Tech Lead</em></p>
            <p>Owner of Redocy Lifecare and co-founder of QRx AI. Raj handles all the coding and back-end development, turning Meera's ideas into working AI solutions.</p>
        </div>
        """, unsafe_allow_html=True)

    with col3:
        st.markdown("""
        <div class="card" style="text-align: center;">
            <img src="dhyey.png" alt="Dhyey Rajveer" style="width:120px; height:120px; object-fit:cover; border-radius:50%; margin-bottom: 10px;">
            <h3>Dhyey Rajveer</h3>
            <p><em>Model Development Support</em></p>
            <p>Key contributor in the creation and fine-tuning of the AI models. His technical insights have been instrumental in QRx AI’s early success.</p>
        </div>
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
        <p>© 2025 QRx Pharmaceutical Consultants. All rights reserved.</p>
    </div>
    """, unsafe_allow_html=True)

# Quality Page - ONLY SHOWS QUALITY CONTENT
elif st.session_state.current_page == 'quality':
    st.markdown('<h2 class="section-title">Pharmaceutical Quality Services</h2>', unsafe_allow_html=True)
    
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
        <p>© 2025 QRx Pharmaceutical Consultants. All rights reserved.</p>
    </div>
    """, unsafe_allow_html=True)

# Close the main content div
st.markdown('</div>', unsafe_allow_html=True)