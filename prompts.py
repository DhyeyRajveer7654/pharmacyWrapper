from string import Template

# Styling for a professional pharma-based UI
TABLE_STYLE = """
<style>
    .table-container {
        display: flex;
        justify-content: center;
        align-items: center;
        margin-top: 20px;
    }
    table {
        width: 90%;
        border-collapse: collapse;
        background-color: #f4f8fb;
        color: #0b3d91;
        border-radius: 10px;
        font-size: 16px;
        border: 1px solid #004080;
        text-align: left;
    }
    th {
        background-color: #007BFF;
        color: white;
        padding: 12px;
        text-align: center;
    }
    td {
        border: 1px solid #004080;
        padding: 10px;
        text-align: left;
    }
    tr:nth-child(even) {
        background-color: #e6f0ff;
    }
    tr:hover {
        background-color: #cce0ff;
    }
</style>
"""

# 📌 **Highly Detailed Method of Preparation**
METHOD_OF_PREPARATION_PROMPT = Template("""
Provide a **detailed step-by-step** **method of preparation** for **$product_name** ($quanOfMed), each containing **$powerOfDrug** based on **$jurisdiction** standards.

Ensure the response covers:
- **Step Number, Description, Equipment Required, Time Duration, Critical Observations, Regulatory Considerations**
- **A table with exact excipient quantities**
- **Justification for ingredient selection and handling**
- **Precautions for uniformity, stability, and pharmacopeial compliance**

$TABLE_STYLE
""")

# 📌 **Characterization & Evaluation**
CHARACTERIZATION_EVALUATION_PROMPT = Template("""
Provide a **detailed characterization and evaluation** for **$product_name** ($quanOfMed), each containing **$powerOfDrug**, based on **$jurisdiction** standards.

Include:
- **Physical characteristics, Identification tests (IR, UV, HPLC), Hardness, Friability, Disintegration, Dissolution, Assay of content**
- **Step-by-step testing instructions, required equipment, SOPs, acceptance limits**
$TABLE_STYLE
""")

# 📌 **Check Results**
CHECK_RESULTS_PROMPT = Template("""
Compare the **quality control results** for **$product_name** ($quanOfMed) with the **$jurisdiction** standards.

Include:
- **Test Parameter, User Result, Pharmacopeial Standard Requirement, Deviation Analysis, Corrective Action Plan, Pass/Fail Status**
- **Scientific explanations on parameter importance, failures, and corrective actions**
$TABLE_STYLE
""")

# 📌 **FTIR Spectrum Analysis (Kept from Latest File)**
FTIR_PROMPT = Template("""
Provide an **FTIR spectrum analysis** for **$product_name**.

Include:
- **Wavenumber (cm⁻¹), Functional Group, Peak Description, Significance in Drug Identification, Potential Interferences, Regulatory Considerations**
$TABLE_STYLE
""")

# 📌 **SMILES Notation for Structure Feature**
STRUCTURE_PROMPT = Template("""
Provide the **canonical SMILES notation** for **$product_name** from PubChem.
""")

# 📌 **Prompt Selection**
def getPromptForOptions(options):
    if options['typeOfInfo'] == "METHOD OF PREPARATION":
        return METHOD_OF_PREPARATION_PROMPT.substitute(options)
    elif options['typeOfInfo'] == "CHARACTERIZATION/EVALUATION":
        return CHARACTERIZATION_EVALUATION_PROMPT.substitute(options)
    elif options['typeOfInfo'] == "CHECK RESULTS":
        return CHECK_RESULTS_PROMPT.substitute(options)
    elif options['typeOfInfo'] == "FTIR ANALYSIS":
        return FTIR_PROMPT.substitute(options)
    return ""
