from string import Template

# 📌 **Highly Detailed Method of Preparation with Excipients Quantity**
METHOD_OF_PREPARATION_PROMPT = Template("""
Provide a detailed method for preparing **$product_name** ($quanOfMed), each containing **$powerOfDrug** of the active ingredient, based on **$jurisdiction** standards. 

Ensure the response includes:
- **A list of materials**, including API and excipients with quantities.
- **Purpose of each material** in the formulation.
- **Step-by-step preparation instructions**, covering:
  - Mixing, granulation, drying, lubrication, and compression.
- **Excludes** evaluation, quality control, or testing procedures.

Ensure the response is formatted as **a clean, minimal HTML table** with a **black background, white text, and properly spaced borders**.
""")

# 📌 **Characterization & Evaluation**
CHARACTERIZATION_EVALUATION_PROMPT = Template("""
Provide a detailed **characterization and evaluation** report for **$product_name** ($quanOfMed), each containing **$powerOfDrug**, based on **$jurisdiction** standards.

The response should include:
- **Physical characteristics** (size, shape, color, appearance).
- **Identification tests** (IR, UV, HPLC).
- **Weight variation, hardness, friability, disintegration, dissolution**.
- **Assay of content** for content uniformity.
- **Step-by-step instructions** for each test.
- **Equipment and SOPs** as per **$jurisdiction**.

The response **must be in table format only**, with **black background, white text, and proper formatting**.
""")

# 📌 **Combined Formulation & Testing**
COMBINED_PROMPT = Template("""
Provide a **comprehensive guide** for the **formulation and testing** of **$product_name** ($quanOfMed), each containing **$powerOfDrug**, based on **$jurisdiction** standards.

The response should cover:
1️⃣ **Formulation Process**
   - List of ingredients with **roles and quantities**.
   - **Step-by-step formulation process**.
   
2️⃣ **Testing & Quality Control**
   - **Testing criteria** and required procedures.
   - **Standards & expected results** as per **$jurisdiction**.

Ensure the response **remains in a table format with a dark theme**.
""")

# 📌 **Check Results Against Standards**
CHECK_RESULTS_PROMPT = Template("""
Compare the following **evaluation results** of **$powerOfDrug $product_name** ($quanOfMed) with the **$jurisdiction** standards:

$resultsToCheck

Ensure the response is **formatted as an HTML table only**, with black background and white text.
""")

# 📌 **FTIR Spectrum Analysis (Kept from Newer File)**
FTIR_PROMPT = Template("""
Provide a **detailed FTIR spectrum analysis** for **$product_name**.

Ensure the response includes:
- **Wavenumber (cm⁻¹)**
- **Functional Group**
- **Peak Description**
- **Significance in Drug Identification**
- **Potential Interferences**
- **Regulatory Considerations**

Ensure the response **remains in a dark-themed HTML table**.
""")

# 📌 **Prompt Selection Function**
def getPromptForOptions(options):
    if options['typeOfInfo'] == "METHOD OF PREPARATION":
        return METHOD_OF_PREPARATION_PROMPT.substitute(options)
    elif options['typeOfInfo'] == "CHARACTERIZATION/EVALUATION":
        return CHARACTERIZATION_EVALUATION_PROMPT.substitute(options)
    elif options['typeOfInfo'] == "Both of above":
        return COMBINED_PROMPT.substitute(options)
    elif options['typeOfInfo'] == "CHECK RESULTS":
        return CHECK_RESULTS_PROMPT.substitute(options)
    elif options['typeOfInfo'] == "FTIR ANALYSIS":
        return FTIR_PROMPT.substitute(options)
    return ""
