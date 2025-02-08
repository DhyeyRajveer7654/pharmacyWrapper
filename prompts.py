from string import Template

# Styling for a centered, left-aligned table with white text
TABLE_STYLE = """
<style>
    .table-container {
        display: flex;
        justify-content: center;
        align-items: center;
        margin-top: 20px;
    }
    table {
        width: 85%;
        border-collapse: collapse;
        background-color: #1e1e1e;
        color: white;
        border-radius: 10px;
        font-size: 16px;
        border: 1px solid #444;
        text-align: left;
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
        text-align: left;
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

# Function to format responses as a detailed HTML table
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

# 📌 **Method of Preparation (Highly Detailed)**
METHOD_OF_PREPARATION_PROMPT = Template("""
Provide a **highly detailed** step-by-step **method of preparation** for **$product_name** ($quanOfMed), each containing **$powerOfDrug** of the active ingredient, based on **$jurisdiction** standards.
The response should be a **centered table** with the following columns:
- **Step Number**
- **Description**
- **Equipment Required**
- **Time Duration**
- **Critical Observations**

Ensure that the response is in a **table format only**, with all text **left-aligned** and no extra text outside the table.
""")

# 📌 **Characterization & Evaluation (Highly Detailed)**
CHARACTERIZATION_EVALUATION_PROMPT = Template("""
Provide a **comprehensive, detailed** characterization and evaluation of **$product_name** ($quanOfMed), each containing **$powerOfDrug**, based on **$jurisdiction** standards.
The response should be a **centered table** with the following columns:
- **Test Name**
- **Objective**
- **Procedure**
- **Equipment Used**
- **Acceptable Limits**
- **Expected Outcome**
- **Possible Errors & Solutions**

Ensure that the response is in a **table format only**, with all text **left-aligned** and no extra text outside the table.
""")

# 📌 **Combined Formulation & Testing Process (Highly Detailed)**
COMBINED_PROMPT = Template("""
Provide a **fully detailed** combined method for the **formulation and testing** of **$product_name** ($quanOfMed), each containing **$powerOfDrug**, based on **$jurisdiction** standards.
The response should include **two separate sections in table format**:
1️⃣ **Formulation Process**:
   - **Ingredient**
   - **Quantity**
   - **Purpose**
   - **Mixing Procedure**
   - **Processing Parameters**
2️⃣ **Testing Process**:
   - **Test Name**
   - **Methodology**
   - **Equipment Required**
   - **Expected Limits**
   - **Deviation Handling**

Ensure the response is in **table format only**, with all text **left-aligned** and no extra text outside the table.
""")

# 📌 **Quality Control & Results Checking (Highly Detailed)**
CHECK_RESULTS_PROMPT = Template("""
Compare the following **quality control evaluation results** of **$product_name** ($powerOfDrug) for quantity **$quanOfMed** with the **$jurisdiction** standards.
The response should be a **centered table** with:
- **Test Parameter**
- **User Result**
- **Standard Requirement**
- **Deviation (if any)**
- **Corrective Action**
- **Pass/Fail Status**

Ensure the response is in **table format only**, with all text **left-aligned** and no extra text outside the table.
$resultsToCheck
""")

# 📌 **FTIR Spectrum Analysis (Highly Detailed)**
FTIR_PROMPT = Template("""
Provide a **detailed standard FTIR spectrum data** for **$product_name**.
The response should be a **centered table** with the following columns:
- **Wavenumber (cm⁻¹)**
- **Functional Group**
- **Peak Description**
- **Significance in Drug Identification**
- **Possible Interferences**

Ensure the response is in **table format only**, with all text **left-aligned** and no extra text outside the table.
""")

# 📌 **Dissolution & Stability Studies (Highly Detailed)**
DISSOLUTION_STABILITY_PROMPT = Template("""
Provide a **highly detailed** report on the **dissolution and stability studies** of **$product_name** ($quanOfMed), each containing **$powerOfDrug**, based on **$jurisdiction** standards.
The response should be a **centered table** with:
- **Study Type (Dissolution/Stability)**
- **Test Conditions**
- **Sampling Time Points**
- **Equipment Used**
- **Acceptable Limits**
- **Results Interpretation**
- **Corrective Measures if Deviations Occur**

Ensure the response is in **table format only**, with all text **left-aligned** and no extra text outside the table.
""")

# Function to generate GPT prompt
def getPromptForOptions(options):
    jurisdiction = options['jurisdiction']
    if jurisdiction == "COMPARE WITH ALL OF THEM":
        jurisdiction = "INDIAN PHARMACOPIEA, BRITISH PHARMACOPIEA, MARTINDALE-EXTRA PHAMRACOPIEA and UNITED STATES PHARMACOPOEIA"

    if options['typeOfInfo'] == "METHOD OF PREPARATION":
        prompt = METHOD_OF_PREPARATION_PROMPT.substitute(options)
    elif options['typeOfInfo'] == "CHARACTARIZATION/EVALUATION":
        prompt = CHARACTERIZATION_EVALUATION_PROMPT.substitute(options)
    elif options['typeOfInfo'] == "Both of above":
        prompt = COMBINED_PROMPT.substitute(options)
    elif options['typeOfInfo'] == "CHECK RESULTS":
        prompt = CHECK_RESULTS_PROMPT.substitute(options)
    elif options['typeOfInfo'] == "DISSOLUTION & STABILITY":
        prompt = DISSOLUTION_STABILITY_PROMPT.substitute(options)
    elif options['typeOfInfo'] == "FTIR ANALYSIS":
        prompt = FTIR_PROMPT.substitute(options)
    else:
        prompt = ""

    return prompt
