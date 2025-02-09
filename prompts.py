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

# 📌 **Highly Detailed Method of Preparation with Excipients Quantity**
METHOD_OF_PREPARATION_PROMPT = Template("""
Provide a **highly detailed, step-by-step** **method of preparation** for **$product_name** ($quanOfMed), each containing **$powerOfDrug** of the active ingredient, based on **$jurisdiction** standards.

Ensure the response is a **centered HTML table** covering:
- **Step Number**
- **Step Description**
- **Equipment Required**
- **Time Duration**
- **Critical Observations**
- **Regulatory Considerations**

Additionally, provide a **reference table** showing the **exact quantity** of excipients required based on **$quanOfMed**.  
This table should include:
- **Ingredient Name (API & Excipients)**
- **Required Quantity per Dosage Unit**
- **Total Quantity Required for $quanOfMed**
- **Function in Formulation**
- **Solubility & Stability Considerations**

Each step must include **scientific justification**, including:
- How **ingredients are selected and handled**.
- Precautions to **avoid errors** during mixing, drying, compression, and packaging.
- How to ensure **uniformity, stability, and compliance** with pharmacopeial standards.

The response **must be in HTML table format only**, with **white text inside a dark background**, text **left-aligned**, and no extra text outside the table.
""")

# 📌 **Highly Detailed Combined Formulation & Testing with Excipients Quantity**
COMBINED_PROMPT = Template("""
Provide a **fully detailed** combined **formulation and testing** report for **$product_name** ($quanOfMed), each containing **$powerOfDrug**, based on **$jurisdiction** standards.

The response should include **two separate centered tables**:
1️⃣ **Formulation Process**:
   - **Ingredient**
   - **Quantity per Unit**
   - **Total Quantity for $quanOfMed**
   - **Purpose**
   - **Mixing & Processing Steps**
   - **Critical Processing Parameters**
   - **Possible Risks & Precautions**

2️⃣ **Testing & Quality Control**:
   - **Test Name**
   - **Testing Procedure**
   - **Equipment Used with steps of using equipment**
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

STRUCTURE_PROMPT = Template("""
Provide the **canonical SMILES notation** for the drug $product_name based on PubChem's database. Ensure that the SMILES code is accurate and matches PubChem's standard molecular structure for the drug. Return only the canonical SMILES code as provided by PubChem, and no other extra text. If the drug name is not valid, return only "NO DRUG FOUND".
""")

# 📌 **Fix: Force Output into Proper Table Format**
def format_table(response):
    """Ensure the GPT response is returned as a proper HTML table format."""
    if "<table>" not in response:  # Prevents raw HTML output issue
        rows = response.strip().split("\n")
        formatted_rows = ["<tr><td>" + "</td><td>".join(row.split(":")) + "</td></tr>" for row in rows if ":" in row]
        return f"<table border='1' style='border-collapse:collapse;width:100%'><tr><th>Parameter</th><th>Value</th></tr>{''.join(formatted_rows)}</table>"
    return response  # Return as-is if already a valid table

# 📌 **GPT Prompt Selection**
def getPromptForOptions(options):
    if options['typeOfInfo'] == "METHOD OF PREPARATION":
        return METHOD_OF_PREPARATION_PROMPT.substitute(options)
    elif options['typeOfInfo'] == "Both of above":
        return COMBINED_PROMPT.substitute(options)
    elif options['typeOfInfo'] == "CHECK RESULTS":
        return CHECK_RESULTS_PROMPT.substitute(options)
    elif options['typeOfInfo'] == "DISSOLUTION & STABILITY" or options['typeOfInfo'] == "CHARACTARIZATION/EVALUATION":
        return DISSOLUTION_STABILITY_PROMPT.substitute(options)
    elif options['typeOfInfo'] == "FTIR ANALYSIS":
        return FTIR_PROMPT.substitute(options)
    return ""

