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

# 📌 **Method of Preparation (Ensures Table Output)**
METHOD_OF_PREPARATION_PROMPT = Template("""
$TABLE_STYLE
<table>
    <tr><th>Step Number</th><th>Step Description</th><th>Equipment Required</th><th>Time Duration</th><th>Critical Observations</th><th>Regulatory Considerations</th></tr>
</table>
Generate a **detailed step-by-step method of preparation** for **${product_name}** (${quanOfMed}), each containing **${powerOfDrug}** based on **${jurisdiction}** standards.
""")

# 📌 **Characterization & Evaluation**
CHARACTERIZATION_EVALUATION_PROMPT = Template("""
$TABLE_STYLE
<table>
    <tr><th>Test Name</th><th>Testing Procedure</th><th>Equipment Used</th><th>Acceptance Criteria</th><th>Deviation Handling</th><th>Regulatory Considerations</th></tr>
</table>
Generate a **detailed characterization and evaluation** for **${product_name}** (${quanOfMed}), each containing **${powerOfDrug}**, based on **${jurisdiction}** standards.
""")

# 📌 **Check Results**
CHECK_RESULTS_PROMPT = Template("""
$TABLE_STYLE
<table>
    <tr><th>Test Parameter</th><th>User Result</th><th>Pharmacopeial Standard</th><th>Deviation Analysis</th><th>Corrective Action Plan</th><th>Pass/Fail Status</th></tr>
</table>
Compare the **quality control results** for **${product_name}** (${quanOfMed}) with the **${jurisdiction}** standards.

${resultsToCheck}
""")

# 📌 **FTIR Spectrum Analysis**
FTIR_PROMPT = Template("""
$TABLE_STYLE
<table>
    <tr><th>Wavenumber (cm⁻¹)</th><th>Functional Group</th><th>Peak Description</th><th>Significance in Drug Identification</th><th>Potential Interferences</th><th>Regulatory Considerations</th></tr>
</table>
Generate an **FTIR spectrum analysis** for **${product_name}**.
""")

# 📌 **Structure Feature**
STRUCTURE_PROMPT = Template("""
Provide the **canonical SMILES notation** for **${product_name}** from PubChem.
""")

# 📌 **Fix: Ensure All Required Keys Exist**
def getPromptForOptions(options):
    """Ensures all required keys exist before calling .substitute()"""
    required_keys = ["product_name", "quanOfMed", "powerOfDrug", "jurisdiction", "resultsToCheck"]
    
    # Set default values for missing keys
    for key in required_keys:
        options[key] = options.get(key, "N/A")

    # Debugging: Print the options dictionary (REMOVE AFTER TESTING)
    print("DEBUG: Options dictionary before calling substitute():", options)

    # Select the appropriate prompt
    if options['typeOfInfo'] == "METHOD OF PREPARATION":
        return METHOD_OF_PREPARATION_PROMPT.substitute(options)
    elif options['typeOfInfo'] == "CHARACTERIZATION/EVALUATION":
        return CHARACTERIZATION_EVALUATION_PROMPT.substitute(options)
    elif options['typeOfInfo'] == "CHECK RESULTS":
        return CHECK_RESULTS_PROMPT.substitute(options)
    elif options['typeOfInfo'] == "FTIR ANALYSIS":
        return FTIR_PROMPT.substitute(options)
    
    return "Error: Invalid typeOfInfo selected."
