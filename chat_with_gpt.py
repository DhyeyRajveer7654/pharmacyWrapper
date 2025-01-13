from openai import OpenAI
import streamlit as st


openai_api_key = st.secrets["api"]["key"]


def chatWithGpt(prompt):
    try:
        # Use the newer method (hypothetical example, replace with actual updates)
        client = OpenAI(
            api_key=openai_api_key,  # This is the default and can be omitted
        )
        response = client.chat.completions.create(
            messages=[
                {
                    "role": "user",
                    "content": prompt
                }
            ],
            model="gpt-4o-mini",
        )
        # Extract the response content
        return str(response.choices[0].model_dump(mode='python')['message']['content'])
    except Exception as e:
        return f"Error: {str(e)}"