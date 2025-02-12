from openai import OpenAI
import streamlit as st

def chatWithGpt(prompt):
    try:
        # Initialize OpenAI client with error handling for missing API key
        if 'api' not in st.secrets or 'OPENAI_API_KEY' not in st.secrets['api']:
            st.error("⚠️ OpenAI API key not found. Please check your configuration.")
            return None

        client = OpenAI(api_key=st.secrets['api']['OPENAI_API_KEY'])

        # Make API call
        response = client.chat.completions.create(
            messages=[{"role": "user", "content": prompt}],
            model="gpt-4",  # Using standard GPT-4 model
            temperature=0.7,
            max_tokens=2000
        )
        return str(response.choices[0].message.content)
    except Exception as e:
        st.error(f"⚠️ Error: {str(e)}")
        return None

def get_ftir_from_gpt(drug_name):
    prompt = f"Provide the standard FTIR peak values for {drug_name}, listing wavenumbers (cm⁻¹) and corresponding functional groups in an HTML table format with a dark background and white text."
    ftir_data = chatWithGpt(prompt)
    if ftir_data:
        # Ensure the response is wrapped in the proper table styling
        styled_data = f"""
        <div class="table-container">
            {ftir_data}
        </div>
        """
        return styled_data
    return None