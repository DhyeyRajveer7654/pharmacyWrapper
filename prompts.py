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

# 📌 **Highly Detailed Method of Preparation**
METHOD_OF_PREPARATION_PROMPT = Template("""
Provide a **highly detailed, step-by-step** **method of preparation** for **$product_name** ($quanOfMed), each containing **$powerOfDrug** of the active ingredient, based on **$jurisdiction** standards.

Ensure the response is a **centered HTML table** covering:
- **Step Number**
- **Step Description**
- **Equipment Required**
- **Time Duration**
- **Critical Observations**
- **Regulatory Considerations**

Each step must be detailed with **scientific justification**, including:
- How **ingredients are selected and handled**.
- Precautions to **avoid errors** during mixing, drying, compression, and packaging.
- How to ensure **uniformity, stability, and compliance** with pharmacopeial standards.

The response **must be in table format only** with **white text inside a dark background**, all text **left-aligned**, and no extra text outside the table.
""")

# 📌 **Highly Detailed Characterization & Evaluation**
CHARACTERIZATION_EVALUATION_PROMPT = Template("""
Provide a **comprehensive and detailed** characterization and evaluation of **$product_name** ($quanOfMed), each containing **$powerOfDrug**, based on **$jurisdiction** standards.

Ensure the response is a **centered HTML table** covering:
- **Test Name**
- **Objective**
- **Procedure (Step-by-step)**
- **Equipment Required**
- **Acceptable Limits**
- **Expected Outcomes**
- **Common Errors & Corrective Actions**
- **Regulatory Considerations**

Every test should include **scientific justification** explaining:
- Why the test is critical for **quality assurance**.
- What **deviations indicate** about product failure.
- **How to correct failures** to meet pharmacopeial standards.

The response **must be in table format only**, with **white text inside a dark background**, text **left-aligned**, and no extra text outside the table.
""")

# 📌 **Highly Detailed Combined Formulation & Testing**
COMBINED_PROMPT = Template("""
Provide a **detailed** combined **formulation and testing** report for **$product_name** ($quanOfMed), each containing **$powerOfDrug**, based on **$jurisdiction** standards.

The response should include **two separate centered tables**:
1️⃣ **Formulation Process**:
   - **Ingredient**
   - **Quantity**
   - **Purpose**
   - **Mixing & Processing Steps**
   - **Critical Processing Parameters**
   - **Possible Risks & Precautions**

2️⃣ **Testing & Quality Control**:
   - **Test Name**
   - **Testing Procedure**
   - **Equipment Used**
   - **Acceptance Criteria**
   - **Deviation Handling**
   - **Regulatory Considerations**

The response **must be in table format only**, with **white text inside a dark background**, text **left-aligned**, and no extra text outside the table.
""")

# 📌 **Highly Detailed Quality Control & Results Checking**
CHECK_RESULTS_PROMPT = Template("""
Compare the **quality control evaluation results** of **$product_name** ($powerOfDrug) for **$quanOfMed** with the **$jurisdiction** standards.

Ensure the response is a **centered HTML table** covering:
- **Test Parameter**
- **User Result**
- **Pharmacopeial Standard Requirement**
- **Deviation Analysis**
- **Corrective Action Plan**
- **Pass/Fail Status**

Each parameter must be explained in **scientific depth**, including:
- Why the parameter is **critical for drug quality**.
- What **failures indicate** about formulation issues.
- **How to correct issues** based on pharmacopeial standards.

The response **must be in table format only**, with **white text inside a dark background**, text **left-aligned**, and no extra text outside the table.

$resultsToCheck
""")

# 📌 **Highly Detailed FTIR Spectrum Analysis**
FTIR_PROMPT = Template("""
Provide a **detailed FTIR spectrum analysis** for **$product_name**.

Ensure the response is a **centered HTML table** covering:
- **Wavenumber (cm⁻¹)**
- **Functional Group**
- **Peak Description**
- **Significance in Drug Identification**
- **Potential Interferences**
- **Regulatory Considerations**

Explain:
- How FTIR confirms **drug identity**.
- What **peak deviations** indicate about formulation errors.
- How to **ensure FTIR compliance** with pharmacopeial standards.

The response **must be in table format only**, with **white text inside a dark background**, text **left-aligned**, and no extra text outside the table.
""")

# 📌 **Highly Detailed Dissolution & Stability Studies**
DISSOLUTION_STABILITY_PROMPT = Template("""
Provide a **comprehensive dissolution and stability study** for **$product_name** ($quanOfMed), each containing **$powerOfDrug**, based on **$jurisdiction** standards.

Ensure the response is a **centered HTML table** covering:
- **Study Type (Dissolution/Stability)**
- **Test Conditions**
- **Sampling Time Points**
- **Equipment Used**
- **Acceptance Limits**
- **Stability Period**
- **Corrective Actions for Failures**
- **Regulatory Considerations**

The response **must be in table format only**, with **white text inside a dark background**, text **left-aligned**, and no extra text outside the table.
""")

# 📌 **GPT Prompt Selection**
def getPromptForOptions(options):
    if options['typeOfInfo'] == "METHOD OF PREPARATION":
        return METHOD_OF_PREPARATION_PROMPT.substitute(options)
    elif options['typeOfInfo'] == "CHARACTARIZATION/EVALUATION":
        return CHARACTERIZATION_EVALUATION_PROMPT.substitute(options)
    elif options['typeOfInfo'] == "Both of above":
        return COMBINED_PROMPT.substitute(options)
    elif options['typeOfInfo'] == "CHECK RESULTS":
        return CHECK_RESULTS_PROMPT.substitute(options)
    elif options['typeOfInfo'] == "DISSOLUTION & STABILITY":
        return DISSOLUTION_STABILITY_PROMPT.substitute(options)
    elif options['typeOfInfo'] == "FTIR ANALYSIS":
        return FTIR_PROMPT.substitute(options)
    return ""
