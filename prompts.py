from string import Template

# 🌟 Modernized Table Styling with Better Contrast
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
        background-color: rgba(255, 255, 255, 0.95); /* Light background for easy reading */
        color: #222; /* Dark text for strong contrast */
        border-radius: 10px;
        font-size: 16px;
        text-align: left;
        box-shadow: 0px 4px 10px rgba(0, 0, 0, 0.2);
    }
    th {
        background: linear-gradient(90deg, #ff758c, #ff7eb3); /* Attractive gradient */
        color: white;
        padding: 12px;
        text-align: center;
        font-weight: bold;
        border-radius: 8px;
    }
    td {
        border: 1px solid #ddd;
        padding: 10px;
        text-align: left;
    }
    tr:nth-child(even) {
        background-color: #f8f8f8; /* Soft gray for better readability */
    }
    tr:hover {
        background: #ff7eb3;
        color: white;
    }
</style>
"""

# 📌 **Highly Detailed Method of Preparation**
METHOD_OF_PREPARATION_PROMPT = Template("""
Provide a **step-by-step, highly detailed method of preparation** for **$product_name** ($quanOfMed), each containing **$powerOfDrug**, following **$jurisdiction** standards.

Ensure the response follows a **structured manufacturing process** with these steps:

1️⃣ **Ingredient Selection & Pre-processing**  
   - List all active & inactive ingredients with functions.  
   - Pre-treatment or purification methods.  
   
2️⃣ **Pre-formulation Studies**  
   - Solubility, pH adjustments, compatibility testing.  

3️⃣ **Manufacturing Process (Detailed Steps)**  
   - Exact **processing conditions** (temperature, pressure, mixing time).  
   - **Critical parameters** & precautions.  
   - **Scientific justifications** for each step.  

4️⃣ **Final Processing & Quality Control Checks**  
   - **Final drying, compression, or filling** procedures.  
   - **Quality control tests** before release.  

🔹 **Ensure clarity with scientific justifications and GMP-compliant precautions.**
""")

# 📌 **Detailed Formulation & Testing Report**
COMBINED_PROMPT = Template("""
Generate a **highly structured** combined **formulation and testing** report for **$product_name** ($quanOfMed), based on **$jurisdiction** standards.

The response should include **two structured tables**:
1️⃣ **Formulation Process**  
   - **Ingredient** | **Quantity per Unit** | **Total Quantity for $quanOfMed**  
   - **Purpose & Role in Formulation**  
   - **Mixing Steps & Processing Parameters**  

2️⃣ **Testing & Quality Control**  
   - **Test Name** | **Testing Procedure** | **Equipment Used**  
   - **Acceptance Criteria & Deviation Handling**  
   - **Regulatory Compliance & Stability Considerations**  

🔹 **Ensure precise scientific justifications for ingredient selection and process control.**
""")

# 📌 **Enhanced Quality Control & Results Checking**
CHECK_RESULTS_PROMPT = Template("""
Compare the **quality control evaluation results** of **$product_name** ($powerOfDrug) for **$quanOfMed** with the **$jurisdiction** standards.

Ensure the response is a **well-structured HTML table** covering:
- **Test Parameter**
- **User Result**
- **Pharmacopeial Standard Requirement**
- **Deviation Analysis**
- **Corrective Action Plan**
- **Pass/Fail Status**

🔹 **For each parameter, include scientific justifications on its importance, failure impact, and corrective actions.**
""")

# 📌 **Comprehensive FTIR Spectrum Analysis**
FTIR_PROMPT = Template("""
Provide a **detailed FTIR spectrum analysis** for **$product_name**.

Ensure the response is a **structured HTML table** covering:
- **Wavenumber (cm⁻¹)**
- **Functional Group**
- **Peak Description**
- **Significance in Drug Identification**
- **Potential Interferences**
- **Regulatory Compliance Considerations**

🔹 **Explain how deviations in peak values affect formulation quality and corrective actions to address issues.**
""")

# 📌 **Advanced Dissolution & Stability Studies**
DISSOLUTION_STABILITY_PROMPT = Template("""
Provide a **comprehensive dissolution and stability study** for **$product_name** ($quanOfMed), following **$jurisdiction** standards.

Ensure the response is structured into **two key sections**:
1️⃣ **Dissolution Study**  
   - **Study Type (Single-point or Multi-point Analysis)**  
   - **Test Conditions (pH, buffer system, agitation speed, temperature)**  
   - **Sampling Time Points & Equipment Used**  
   - **Acceptance Limits & Corrective Actions for Failures**  

2️⃣ **Stability Study**  
   - **Storage Conditions (Temperature, Humidity, Light Sensitivity)**  
   - **Stability Period & Observations**  
   - **Impact on Drug Potency & Safety**  
   - **GMP Compliance Considerations**  

🔹 **Ensure highly structured formatting with scientific justifications for each decision.**
""")

# 📌 **GPT Prompt Selection**
def getPromptForOptions(options):
    if options['typeOfInfo'] == "METHOD OF PREPARATION":
        return METHOD_OF_PREPARATION_PROMPT.substitute(options)
    elif options['typeOfInfo'] == "Both of above":
        return COMBINED_PROMPT.substitute(options)
    elif options['typeOfInfo'] == "CHECK RESULTS":
        return CHECK_RESULTS_PROMPT.substitute(options)
    elif options['typeOfInfo'] == "DISSOLUTION & STABILITY":
        return DISSOLUTION_STABILITY_PROMPT.substitute(options)
    elif options['typeOfInfo'] == "FTIR ANALYSIS":
        return FTIR_PROMPT.substitute(options)
    return ""
