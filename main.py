import streamlit as st
import pandas as pd
import prompts
from chat_with_gpt import get_gpt_response

# Set page config
st.set_page_config(page_title="QAI Model", page_icon="🔬", layout="wide")

# Custom Styling
st.markdown(
    """
    <style>
        body {
            background-color: #E0F7FA;
        }
        .stTextInput, .stTextArea, .stSelectbox, .stNumberInput {
            border-radius: 10px;
            border: 1px solid #0078D4;
        }
        .report-container {
            background: white;
            padding: 15px;
            border-radius: 10px;
            box-shadow: 0px 4px 10px rgba(0,0,0,0.1);
        }
        .title {
            font-size: 24px;
            font-weight: bold;
            color: #0078D4;
            margin-bottom: 20px;
        }
    </style>
    """,
    unsafe_allow_html=True
)

st.title("🔬 QAI Model - Pharmaceutical Quality Analysis")

# User Inputs
st.sidebar.header("⚙️ Configuration")
option = st.sidebar.selectbox("Select Evaluation Type", ["Characterization", "Preparation Method", "FTIR Analysis"])

input_text = st.text_area("Enter your formulation details:")

if st.button("Generate Report"):
    if input_text:
        with st.spinner("Generating your report..."):
            prompt = prompts.getPromptForOptions(option)
            response = get_gpt_response(prompt + " " + input_text)
            
            # Format results in a table
            if option == "Characterization":
                df = pd.DataFrame([response.split("\n")], columns=["Property", "Value"])
                st.table(df)
            else:
                st.markdown("### 📜 Report")
                st.markdown(f"<div class='report-container'>{response}</div>", unsafe_allow_html=True)
    else:
        st.warning("⚠️ Please enter formulation details.")
