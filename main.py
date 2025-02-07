import streamlit as st
import streamlit.components.v1 as components
import prompts
import chat_with_gpt

# Set Page Configuration
st.set_page_config(page_title="QAI Model", layout="wide", page_icon="🧪")

# Force Custom CSS Styling for UI
st.markdown(
    """
    <style>
        /* Set background color */
        body { background-color: #0e1117; color: white; font-family: 'Arial', sans-serif; }
        
        /* Style input fields */
        .stTextInput>div>div>input, .stSelectbox>div>div>select, .stTextArea>div>textarea { 
            background-color: #1e222a !important; 
            color: white !important; 
            border-radius: 10px !important; 
            padding: 10px;
        }

        /* Style buttons */
        .stButton>button { 
            background: linear-gradient(90deg, #007BFF, #00D4FF); 
            color: white; 
            border-radius: 10px; 
            font-size: 16px; 
            padding: 10px; 
            font-weight: bold; 
            border: none;
            transition: 0.3s;
        }
        
        .stButton>button:hover { 
            background: linear-gradient(90deg, #00D4FF, #007BFF);
            transform: scale(1.05);
        }

        /* Title Styling */
        .stMarkdown h1, .stMarkdown h2, .stMarkdown h3 { 
            color: #00D4FF; 
            text-align: center; 
        }

        /* Card Styling */
        .card {
            background-color: #1e222a; 
            padding: 20px; 
            border-radius: 12px; 
            box-shadow: 0px 4px 10px rgba(255, 255, 255, 0.1);
            margin: 20px;
        }
    </style>
    """,
    unsafe_allow_html=True
)

# Form Page UI
if "page" not in st.session_state:
    st.session_state.page = "form"
if "api_response" not in st.session_state:
    st.session_state.api_response = None

options = dict()

# Header
st.markdown("## 🧪 QAI Model - AI-Powered Quality Assurance")
st.markdown("#### 🔍 Enter the details below to analyze pharmaceutical quality.")

# User Input Form
with st.form("input_form"):
    col1, col2 = st.columns(2)

    with col1:
        options["product_name"] = st.text_input("📌 Product Name", placeholder="e.g., Paracetamol")
        options["powerOfDrug"] = st.text_input("⚡ Power of Drug", placeholder="e.g., 500 mg")
        
    with col2:
        options["quanOfMed"] = st.text_input("📦 Quantity of Medicine", placeholder="e.g., 1000 tablets")
        options["jurisdiction"] = st.selectbox("🌎 Select Jurisdiction", 
            ["INDIAN PHARMACOPIEA", "BRITISH PHARMACOPIEA", "UNITED STATES PHARMACOPOEIA", "COMPARE WITH ALL"])

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
        # Generate prompt
        prompt = prompts.getPromptForOptions(options)

        # Call GPT API
        with st.spinner("🛠️ Processing... Please wait"):
            api_response = chat_with_gpt.chatWithGpt(prompt)
            st.session_state.api_response = api_response

        # Save user input
        st.session_state.update(options)

        # Navigate to result page
        st.write("✅ Data submitted successfully! Click 'Submit' again to see results.")

# Display Result Page
if "api_response" in st.session_state and st.session_state.api_response:
    st.markdown("## 📋 Generated Report")
    components.html(st.session_state.api_response, height=1000, width=1000, scrolling=True)

    if st.session_state.ftir_required:
        with st.spinner("📡 Fetching FTIR Data..."):
            ftir_data = chat_with_gpt.get_ftir_from_gpt(st.session_state.product_name)
            st.markdown("### 🔬 FTIR Data from GPT-4")
            st.write(ftir_data)

    if st.button("🔙 Go Back to Form"):
        st.session_state.clear()  
        st.experimental_rerun()
