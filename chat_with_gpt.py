from openai import OpenAI
import streamlit as st
import requests

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

def fetch_image_from_pubchem(product_name):
    search_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{product_name}/PNG"
    response = requests.get(search_url)
    if response.status_code == 200:
        return search_url
    return None

def fetch_structure_image(product_name):
    pubchem_image = fetch_image_from_pubchem(product_name)
    if pubchem_image:
        return pubchem_image
    prompt = f"Find and provide a direct URL to the chemical structure image of {product_name} from trusted sources like PubChem or PubMed."
    return chatWithGpt(prompt)
