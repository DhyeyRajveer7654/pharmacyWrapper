from string import Template

# Styling for a centered table with white text
TABLE_STYLE = """
<style>
    .table-container {
        display: flex;
        justify-content: center;
        align-items: center;
        margin-top: 20px;
    }
    table {
        width: 80%;
        border-collapse: collapse;
        background-color: #1e1e1e;
        color: white;
        border-radius: 10px;
        font-size: 16px;
        border: 1px solid #444;
        text-align: center;
    }
    th {
        background-color: #007BFF;
        color: white;
        padding: 12px;
        text-align: center;
    }
    td {
        border: 1px solid #444;
        padding: 10px;
        text-align: center;
    }
    tr:nth-child(even) {
        background-color: #292b2c;
    }
    tr:hover {
        background-color: #007BFF;
        color: white;
    }
</style>
"""

# Function to format responses as a centered HTML table
def format_as_table(data):
    table_html = TABLE_STYLE + '<div class="table-container"><table><tr>'
    headers = data[0].keys()

    for header in headers:
        table_html += f"<th>{header}</th>"
    table_html += "</tr>"

    for row in data:
        table_html += "<tr>"
        for value in row.values():
            table_html += f"<td>{value}</td>"
        table_html += "</tr>"

    table_html += "</table></div>"
    return table_html

# GPT Prompts for different sections
METHOD_OF_PREPARATION_PROMPT = Template("""
Provide a detailed **method of preparation** for **$product_name** ($quanOfMed), each containing **$powerOfDrug** of the active ingredient, based on **$jurisdiction** standards.
Format the response as a **centered table** with the following columns:
- **Ingredient**
- **Quantity**
- **Purpose**
- **Preparation Steps**
Ensure the response is **only** in HTML table format with no extra text.
""")

CHARACTERIZATION_EVALUATION_PROMPT = Template("""
Provide a detailed **characterization and evaluation** of **$product_name** ($quanOfMed), each containing **$powerOfDrug**, based on **$jurisdiction** standards.
Format the response as a **centered table** with the following columns:
- **Test Name**
- **Procedure**
- **Acceptable Limits**
- **Expected Outcome**
Ensure the response is **only** in HTML table format with no extra text.
""")

COMBINED_PROMPT = Template("""
Provide a **comprehensive guide** for the formulation and testing of **$product_name** ($quanOfMed), each containing **$powerOfDrug**, based on **$jurisdiction** standards.
Format the response as a **centered table** with two sections:
1️⃣ **Formulation Process** (List ingredients, quantity, purpose).  
2️⃣ **Testing Process** (List test name, method, expected results).  
Ensure the response is **only** in HTML table format with no extra text.
""")

CHECK_RESULTS_PROMPT = Template("""
Compare the following **evaluation results** of **$product_name** ($powerOfDrug) for quantity **$quanOfMed** with the **$jurisdiction** standards.
Provide the response as a **centered table** with:
- **Test Parameter**
- **User Result**
- **Standard Requirement**
- **Pass/Fail**
Ensure the response is **only** in HTML table format with no extra text.
$resultsToCheck
""")

FTIR_PROMPT = Template("""
Provide the **standard FTIR spectrum data** for **$product_name**.
Format the response as a **centered table** with:
- **Wavenumber (cm⁻¹)**
- **Functional Group**
- **Peak Description**
Ensure the response is **only** in HTML table format with no extra text.
""")

# Function to generate GPT prompt
def getPromptForOptions(options):
    jurisdiction = options['jurisdiction']
    if jurisdiction == "COMPARE WITH ALL OF THEM":
        jurisdiction = "INDIAN PHARMACOPIEA, BRITISH PHARMACOPIEA, and UNITED STATES PHARMACOPOEIA"

    if options['typeOfInfo'] == "METHOD OF PREPARATION":
        prompt = METHOD_OF_PREPARATION_PROMPT.substitute(options)
    elif options['typeOfInfo'] == "CHARACTARIZATION/EVALUATION":
        prompt = CHARACTERIZATION_EVALUATION_PROMPT.substitute(options)
    elif options['typeOfInfo'] == "Both of above":
        prompt = COMBINED_PROMPT.substitute(options)
    elif options['typeOfInfo'] == "CHECK RESULTS":
        prompt = CHECK_RESULTS_PROMPT.substitute(options)
    else:
        prompt = ""

    return prompt
