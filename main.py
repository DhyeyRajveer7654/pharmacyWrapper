import streamlit as st
import streamlit.components.v1 as components
import prompts
import chat_with_gpt

# Set Page Configuration
st.set_page_config(page_title="QAI Model", layout="wide", page_icon="🧪")

# ✅ **Modern Styling**
st.markdown("""
    <style>
        /* Transparent Background */
        .stApp {
            background: url('https://source.unsplash.com/1600x900/?science,technology') no-repeat center center fixed;
            background-size: cover;
        }

        /* Report Container */
        .report-container {
            background: rgba(255, 255, 255, 0.2);
            backdrop-filter: blur(10px);
            padding: 30px;
            border-radius: 12px;
            box-shadow: 0px 4px 10px rgba(255, 255, 255, 0.2);
        }

        /* Glowing Blue Headings */
        .title {
            text-align: center;
            font-size: 32px;
            font-weight: bold;
            color: #00BFFF;
            text-shadow: 0px 0px 10px #00BFFF;
        }

        /* Print Button */
        .print-button {
            position: absolute;
            top: 10px;
            right: 20px;
            background: linear-gradient(90deg, #007BFF, #00D4FF);
            color: white;
            padding: 10px 15px;
            border-radius: 10px;
            font-size: 16px;
            font-weight: bold;
            cursor: pointer;
            border: none;
            transition: 0.3s;
        }
        .print-button:hover {
            background: linear-gradient(90deg, #00D4FF, #007BFF);
            transform: scale(1.05);
        }

        /* Hide Print Button when printing */
        @media print {
            .print-button { display: none; }
        }
    </style>
""", unsafe_allow_html=True)

# Page Navigation
if "page" not in st.session_state:
    st.session_state.page = "form"
if "api_response" not in st.session_state:
    st.session_state.api_response = None

options = {}

# 📌 FORM PAGE
if st.session_state.page == "form":
    st.markdown('<div class="title">🧪 QAI Model - AI-Powered Quality Assurance</div>', unsafe_allow_html=True)

    with st.form("input_form"):
        col1, col2 = st.columns(2)

        with col1:
            options["product_name"] = st.text_input("💊 Product Name", placeholder="e.g., Paracetamol")
            options["powerOfDrug"] = st.text_input("⚡ Power of Drug", placeholder="e.g., 500 mg")

        with col2:
            options["quanOfMed"] = st.text_input("📦 Quantity of Medicine", placeholder="e.g., 1000 tablets")
            options["jurisdiction"] = st.selectbox("🌎 Select Jurisdiction", 
                ["INDIAN PHARMACOPIEA", "BRITISH PHARMACOPIEA", "UNITED STATES PHARMACOPOEIA", "COMPARE WITH ALL"])

        options["typeOfInfo"] = st.radio("📊 Select Information Required:", 
                ["METHOD OF PREPARATION", "CHARACTERIZATION/EVALUATION", "Both of above", "CHECK RESULTS"])

        if options["typeOfInfo"] == "CHECK RESULTS":
            options["resultsToCheck"] = st.text_area("🔍 Enter Your Results:", height=150, placeholder="Paste lab results here...")

        options["ftir_required"] = st.checkbox("📡 Retrieve FTIR Data")

        submit_button = st.form_submit_button("🚀 Generate Report")

    if submit_button:
        if not all([options["product_name"], options["quanOfMed"], options["powerOfDrug"]]):
            st.error("⚠️ Please fill in all required fields!")
        else:
            prompt = prompts.getPromptForOptions(options)
            with st.spinner("🛠️ Processing... Please wait"):
                api_response = chat_with_gpt.chatWithGpt(prompt)
                st.session_state.api_response = api_response

            st.session_state.update(options)
            st.session_state.page = "result"
            st.experimental_rerun()

# 📌 RESULT PAGE
elif st.session_state.page == "result":
    st.markdown('<div class="title">📑 Submission Summary</div>', unsafe_allow_html=True)

    st.markdown(f"**💊 Product Name:** {st.session_state.product_name}")
    st.markdown(f"**📦 Quantity of Medicine:** {st.session_state.quanOfMed}")
    st.markdown(f"**⚡ Power of Drug:** {st.session_state.powerOfDrug}")

    # ✅ **Working Print Button**
    print_js = """
        <script>
            function printReport() {
                var divContents = document.getElementById("report").innerHTML;
                var newWindow = window.open('', '', 'height=900, width=1200');
                newWindow.document.write('<html><head><title>Report</title></head><body>');
                newWindow.document.write(divContents);
                newWindow.document.write('</body></html>');
                newWindow.document.close();
                newWindow.print();
            }
        </script>
    """

    st.markdown(print_js, unsafe_allow_html=True)
    st.markdown('<button onclick="printReport()" class="print-button">🖨️ Print Report</button>', unsafe_allow_html=True)

    # ✅ **Display Report in a TABLE**
    st.markdown('<div id="report">', unsafe_allow_html=True)
    if st.session_state.api_response:
        st.markdown(prompts.TABLE_STYLE, unsafe_allow_html=True)
        components.html(st.session_state.api_response, height=1000, width=1000, scrolling=True)
    else:
        st.warning("⚠️ No response received from GPT API.")
    st.markdown('</div>', unsafe_allow_html=True)

    if st.button("🔙 Go Back to Form"):
        st.session_state.page = "form"
        st.experimental_rerun()
