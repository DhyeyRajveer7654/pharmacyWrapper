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

FTIR_PROMPT = Template("""
Provide a **detailed FTIR spectrum analysis** for **$product_name**.

Ensure the response is a **professionally styled HTML table** covering:
- **Wavenumber (cm⁻¹)**
- **Functional Group**
- **Peak Description**
- **Significance in Drug Identification**
- **Potential Interferences**
- **Regulatory Considerations**

### **💠 Table Formatting Requirements**:
- Use a **dark blue header** (`#0B3D91`) with **white text**.
- Apply **alternating row colors**: 
  - Even rows: **Light blue (`#E3F2FD`)**
  - Odd rows: **White**
- Ensure **left-aligned text** for readability.
- Include **hover effect** (`#CFE2FF`) for better visibility.
- **No extra text** outside the table.

### **💡 Key Considerations for FTIR Analysis**:
- Explain how FTIR confirms **drug identity**.
- Identify what **peak deviations** indicate about formulation errors.
- Provide guidance on ensuring **FTIR compliance** with pharmacopeial standards.

### **⚡ Response Format (STRICTLY TABLE ONLY, NO EXTRA TEXT)**:
""")  # 🛠️ Fixed: Properly closed the Template string

# ✅ Fixed formatting issues in the response table
HTML_TABLE_FORMAT = """
<style>
    table {
        width: 100%;
        border-collapse: collapse;
        font-family: 'Arial', sans-serif;
        font-size: 16px;
        border: 1px solid black;
    }
    th {
        background-color: #0B3D91; /* Dark Blue Header */
        color: white;
        font-weight: bold;
        font-style: italic;
        padding: 12px;
        text-align: left;
        border: 1px solid black;
    }
    td {
        padding: 12px;
        text-align: left;
        border: 1px solid black;
        color: black;
    }
    tr:nth-child(even) {
        background-color: #E3F2FD; /* Light Blue */
    }
    tr:nth-child(odd) {
        background-color: white;
    }
    tr:hover {
        background-color: #CFE2FF; /* Slightly darker blue hover effect */
    }
</style>
<table>
    <tr>
        <th>Wavenumber (cm⁻¹)</th>
        <th>Functional Group</th>
        <th>Peak Description</th>
        <th>Significance in Drug Identification</th>
        <th>Potential Interferences</th>
        <th>Regulatory Considerations</th>
    </tr>
    <tr>
        <td>Example 1</td>
        <td>Functional Group 1</td>
        <td>Description 1</td>
        <td>Significance 1</td>
        <td>Interferences 1</td>
        <td>Regulatory 1</td>
    </tr>
    <tr>
        <td>Example 2</td>
        <td>Functional Group 2</td>
        <td>Description 2</td>
        <td>Significance 2</td>
        <td>Interferences 2</td>
        <td>Regulatory 2</td>
    </tr>
</table>
"""  # 🛠️ Fixed: Ensured proper closure of multi-line string

# ✅ Corrected Python syntax errors by separating the definitions properly

