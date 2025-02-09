from string import Template

# 📌 **Table Styling for Pharma-Themed Reports**
TABLE_STYLE = """
<style>
    .table-container { display: flex; justify-content: center; align-items: center; margin-top: 20px; }
    table { width: 90%; border-collapse: collapse; background-color: #f4f8fb; color: #0b3d91; border-radius: 10px; font-size: 16px; border: 1px solid #004080; text-align: left; }
    th { background-color: #007BFF; color: white; padding: 12px; text-align: center; }
    td { border: 1px solid #004080; padding: 10px; text-align: left; }
    tr:nth-child(even) { background-color: #e6f0ff; }
    tr:hover { background-color: #cce0ff; }
</style>
"""

# 📌 **Method of Preparation (Forces Table Output)**
METHOD_OF_PREPARATION_PROMPT = Template("""
Generate a **detailed step-by-step method of preparation** for **$product_name** ($quanOfMed), each containing **$powerOfDrug** based on **$jurisdiction** standards.

Ensure the response follows this **strict table format**:
$TABLE_STYLE
<table>
    <tr><th>Step Number</th><th>Step Description</th><th>Equipment Required</th><th>Time Duration</th><th>Critical Observations</th><th>Regulatory Considerations</th></tr>
</table>

Additionally, generate a **Reference Table** with excipient details in the same format.
""")

# 📌 **Characterization & Evaluation (Forces Table Output)**
CHARACTERIZATION_EVALUATION_PROMPT = Template("""
Generate a **detailed characterization and evaluation** for **$product_name** ($quanOfMed), each containing **$powerOfDrug**, based on **$jurisdiction** standards.

Ensure the response follows this **strict table format**:
$TABLE_STYLE
<table>
    <tr><th>Test Name</th><th>Testing Procedure</th><th>Equipment Used</th><th>Acceptance Criteria</th><th>Deviation Handling</th><th>Regulatory Considerations</th></tr>
</table>
""")

# 📌 **Check Results (Forces Table Output)**
CHECK_RESULTS_PROMPT = Template("""
Compare the **quality control results** for **$product_name** ($quanOfMed) with the **$jurisdiction** standards.

Ensure the response follows this **strict table format**:
$TABLE_STYLE
<table>
    <tr><th>Test Parameter</th><th>User Result</th><th>Pharmacopeial Standard</th><th>Deviation Analysis</th><th>Corrective Action Plan</th><th>Pass/Fail Status</th></tr>
</table>
""")

# 📌 **FTIR Spectrum Analysis (Forces Table Output)**
FTIR_PROMPT = Template("""
Generate an **FTIR spectrum analysis** for **$product_name**.

Ensure the response follows this **strict table format**:
$TABLE_STYLE
<table>
    <tr><th>Wavenumber (cm⁻¹)</th><th>Functional Group</th><th>Peak Description</th><th>Significance in Drug Identification</th><th>Potential Interferences</th><th>Regulatory Considerations</th></tr>
</table>
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
