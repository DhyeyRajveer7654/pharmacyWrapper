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

# Characterization/Evaluation Prompt
CHARACTERIZATION_PROMPT = Template("""
Provide a **comprehensive characterization and evaluation** report for **$product_name** ($quanOfMed), each containing **$powerOfDrug**, based on **$jurisdiction** standards.

Ensure the response is in **table format only** covering:
- **Characterization Parameter**
- **Testing Methodology**
- **Equipment Used**
- **Standard Requirements**
- **Observed Results**
- **Deviation & Correction Measures**
""")

# Keep all the prompts as they are
METHOD_OF_PREPARATION_PROMPT = Template("""
Provide a **highly detailed, step-by-step** **method of preparation** for **$product_name** ($quanOfMed), each containing **$powerOfDrug** of the active ingredient, based on **$jurisdiction** standards.

Ensure the response is a **table format only** covering:
- **Step Number**
- **Step Description**
- **Equipment Required**
- **Time Duration**
- **Critical Observations**
- **Regulatory Considerations**

Additionally, include a **reference table** for excipients:
- **Ingredient Name**
- **Required Quantity per Dosage Unit**
- **Total Quantity Required for $quanOfMed**
- **Function in Formulation**
- **Solubility & Stability Considerations**
""")

COMBINED_PROMPT = Template("""
Provide a **detailed** combined **formulation and testing** report for **$product_name** ($quanOfMed), each containing **$powerOfDrug**, based on **$jurisdiction** standards.

Include **two separate tables**:
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
   - **Equipment Used**
   - **Acceptance Criteria**
   - **Deviation Handling**
   - **Regulatory Considerations**
""")

CHECK_RESULTS_PROMPT = Template("""
Compare the **quality control evaluation results** of **$product_name** ($powerOfDrug) for **$quanOfMed** with the **$jurisdiction** standards.

Ensure the response is in **table format only** covering:
- **Test Parameter**
- **User Result**
- **Pharmacopeial Standard Requirement**
- **Deviation Analysis**
- **Corrective Action Plan**
- **Pass/Fail Status**

$resultsToCheck
""")

FTIR_PROMPT = Template("""
Provide a **detailed FTIR spectrum analysis** for **$product_name**.

Ensure the response is in **table format only** covering:
- **Wavenumber (cm⁻¹)**
- **Functional Group**
- **Peak Description**
- **Significance in Drug Identification**
- **Potential Interferences**
- **Regulatory Considerations**
""")

DISSOLUTION_STABILITY_PROMPT = Template("""
Provide a **comprehensive dissolution and stability study** for **$product_name** ($quanOfMed), each containing **$powerOfDrug**, based on **$jurisdiction** standards.

Ensure the response is in **table format only** covering:
- **Study Type (Dissolution/Stability)**
- **Test Conditions**
- **Sampling Time Points**
- **Equipment Used**
- **Acceptance Limits**
- **Stability Period**
- **Corrective Actions for Failures**
- **Regulatory Considerations**
""")

STRUCTURE_PROMPT = Template("""
Provide the **canonical SMILES notation** for **$product_name** based on PubChem's database.

Return **only** the canonical SMILES code. If the drug name is not valid, return **"NO DRUG FOUND"**.
""")

# GPT Prompt Selection
def getPromptForOptions(options):
    if options['typeOfInfo'] == "METHOD OF PREPARATION":
        return METHOD_OF_PREPARATION_PROMPT.substitute(options)
    elif options['typeOfInfo'] == "Both of above":
        return COMBINED_PROMPT.substitute(options)
    elif options['typeOfInfo'] == "CHECK RESULTS":
        return CHECK_RESULTS_PROMPT.substitute(options)
    elif options['typeOfInfo'] == "DISSOLUTION & STABILITY":
        return DISSOLUTION_STABILITY_PROMPT.substitute(options)
    elif options['typeOfInfo'] == "CHARACTERIZATION/EVALUATION":
        return CHARACTERIZATION_PROMPT.substitute(options)
    elif options['typeOfInfo'] == "FTIR ANALYSIS":
        return FTIR_PROMPT.substitute(options)
    return ""
