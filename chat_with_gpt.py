from openai import OpenAI
import streamlit as st

def chatWithGpt(prompt):
    """
    Send a prompt to OpenAI's API and get a response.
    """
    try:
        if not st.secrets["api"]["key"]:
            return "Error: OpenAI API key not configured. Please check your secrets."

        client = OpenAI(api_key=st.secrets["api"]["key"])
        response = client.chat.completions.create(
            messages=[{"role": "user", "content": prompt}],
            model="gpt-4o",  # Using the latest model as of May 13, 2024
            temperature=0.7,
            max_tokens=2000
        )
        return str(response.choices[0].message.content)
    except Exception as e:
        error_msg = str(e)
        if "Invalid API key" in error_msg:
            return "Error: Invalid OpenAI API key. Please check your configuration."
        elif "Rate limit" in error_msg:
            return "Error: API rate limit exceeded. Please try again later."
        else:
            return f"Error: {error_msg}"

def get_ftir_from_gpt(drug_name):
    """
    Get FTIR data for a specific drug using OpenAI API.
    """
    if not drug_name:
        return "Error: Please provide a valid drug name."

    try:
        prompt = f"""
        Provide a detailed FTIR spectral analysis for {drug_name} in an HTML table format.
        Include:
        - Wavenumber ranges (cm⁻¹)
        - Corresponding functional groups
        - Peak characteristics
        - Significance in drug identification

        Format the response as a professional-looking HTML table with proper styling.
        """
        return chatWithGpt(prompt)
    except Exception as e:
        return f"Error retrieving FTIR data: {str(e)}"