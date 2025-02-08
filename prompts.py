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

# 📌 **Highly Detailed Method of Preparation with Excipients Quantity**
METHOD_OF_PREPARATION_PROMPT = Template("""
Provide a **highly detailed, step-by-step** **method of preparation** for **$product_name** ($quanOfMed), each containing **$powerOfDrug** of the active ingredient, based on **$jurisdiction** standards.

Ensure the response is a **well-formatted HTML table** covering:
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

The response **must be in an easy-to-read HTML table**, with **hover effects and high contrast**.
""")

# 📌 **Highly Detailed Combined Formulation & Testing with Excipients Quantity**
COMBINED_PROMPT = Template("""
Provide a **fully detailed** combined **formulation and testing** report for **$product_name** ($quanOfMed), each containing **$powerOfDrug**, based on **$jurisdiction** standards.

The response should include **two separate tables**:
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

The response **must be in a visually appealing table format**, with **modern styling and high readability**.
""")

# 📌 **Highly Detailed Quality Control & Results Checking**
CHECK_RESULTS_PROMPT = Template("""
Compare the **quality control evaluation results** of **$product_name** ($powerOfDrug) for **$quanOfMed** with the **$jurisdiction** standards.

Ensure the response is a **cleanly formatted HTML table** covering:
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

The response **must be in an easy-to-read HTML table**, ensuring **high contrast and clarity**.
""")

# 📌 **Highly Detailed FTIR Spectrum Analysis**
FTIR_PROMPT = Template("""
Provide a **detailed FTIR spectrum analysis** for **$product_name**.

Ensure the response is a **clear, formatted HTML table** covering:
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

The response **must be in an easy-to-read table**, with **alternating row colors and gradient headers**.
""")

# 📌 **Highly Detailed Dissolution & Stability Studies**
DISSOLUTION_STABILITY_PROMPT = Template("""
Provide a **comprehensive dissolution and stability study** for **$product_name** ($quanOfMed), each containing **$powerOfDrug**, based on **$jurisdiction** standards.

Ensure the response is a **clear and structured HTML table** covering:
- **Study Type (Dissolution/Stability)**
- **Test Conditions**
- **Sampling Time Points**
- **Equipment Used**
- **Acceptance Limits**
- **Stability Period**
- **Corrective Actions for Failures**
- **Regulatory Considerations**

The response **must be in an easy-to-read table format**, with **hover effects and clear contrast**.
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
