from openai import OpenAI
import streamlit as st

openai_api_key = st.secrets["api"]["key"]

def chatWithGpt(prompt):
    try:
        client = OpenAI(api_key=openai_api_key)
        response = client.chat.completions.create(
            messages=[{"role": "user", "content": prompt}],
            model="gpt-4o",
        )
        return str(response.choices[0].model_dump(mode='python')['message']['content'])
    except Exception as e:
        return f"Error: {str(e)}"

# ✅ New function to get FTIR data
def get_ftir_from_gpt(drug_name):
    prompt = f"Provide the standard FTIR peak values for {drug_name}, listing wavenumbers (cm⁻¹) and corresponding functional groups."
    return chatWithGpt(prompt)
