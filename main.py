import streamlit as st
import streamlit.components.v1 as components
import prompts
import chat_with_gpt

# Set Page Configuration
st.set_page_config(page_title="QAI Model", layout="wide", page_icon="🧪")

# ✅ **Modern Styling**
st.markdown("""
    <style>
        .stApp {
            background: url('https://source.unsplash.com/1600x900/?science,technology') no-repeat center center fixed;
            background-size: cover;
        }
        .title {
            text-align: center;
            font-size: 32px;
            font-weight: bold;
            color: #00BFFF;
            text-shadow: 0px 0px 10px #00BFFF;
        }
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
        @media print {
            .print-button { display: none; }
        }
    </style>
""", unsafe_allow_html=True)

# ✅ **Fixing Print Button (Now Works Every Time)**
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

# Page Navigation
if "page" not in st.session_state:
    st.session_state.page = "form"
if "api_response" not in st.session_state:
    st.session_state.api_response = None

# 📌 FORM PAGE
if st.session_state.page == "form":
    st.markdown('<div class="title">🧪 QAI Model - AI-Powered Quality Assurance</div>', unsafe_allow_html=True)

    with st.form("input_form"):
        options = {}
        options["product_name"] = st.text_input("💊 Product Name")
        options["powerOfDrug"] = st.text_input("⚡ Power of Drug")
        options["quanOfMed"] = st.text_input("📦 Quantity of Medicine")
        options["typeOfInfo"] = st.radio("📊 Select Information Required:", 
                ["METHOD OF PREPARATION", "CHARACTERIZATION/EVALUATION", "Both of above", "CHECK RESULTS"])
        options["ftir_required"] = st.checkbox("📡 Retrieve FTIR Data")

        submit_button = st.form_submit_button("🚀 Generate Report")

    if submit_button:
        prompt = prompts.getPromptForOptions(options)
        with st.spinner("🛠️ Processing... Please wait"):
            api_response = chat_with_gpt.chatWithGpt(prompt)
            st.session_state.api_response = api_response

        st.session_state.page = "result"
        st.experimental_rerun()

# 📌 RESULT PAGE
elif st.session_state.page == "result":
    st.markdown(print_js, unsafe_allow_html=True)
    st.markdown('<button onclick="printReport()" class="print-button">🖨️ Print Report</button>', unsafe_allow_html=True)

    # ✅ **Displaying Report in a TABLE**
    st.markdown('<div id="report">', unsafe_allow_html=True)
    if st.session_state.api_response:
        st.markdown(prompts.TABLE_STYLE, unsafe_allow_html=True)
        components.html(st.session_state.api_response, height=800, width=1000, scrolling=True)
    else:
        st.warning("⚠️ No response received from GPT API.")
    st.markdown('</div>', unsafe_allow_html=True)

    if st.button("🔙 Go Back to Form"):
        st.session_state.page = "form"
        st.experimental_rerun()
