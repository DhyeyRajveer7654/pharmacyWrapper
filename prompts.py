from string import Template

# ✅ **Table Styling for Better Readability**
TABLE_STYLE = """
<style>
    table {
        width: 100%;
        border-collapse: collapse;
        background-color: white;
        color: black;
        border-radius: 10px;
        font-size: 16px;
        text-align: left;
        box-shadow: 0px 4px 10px rgba(0, 0, 0, 0.2);
    }
    th {
        background: #007BFF;
        color: white;
        padding: 12px;
        text-align: center;
        font-weight: bold;
    }
    td {
        border: 1px solid #ddd;
        padding: 10px;
        text-align: left;
    }
    tr:nth-child(even) {
        background-color: #f8f8f8;
    }
    tr:hover {
        background: #007BFF;
        color: white;
    }
</style>
"""

# 📌 **Highly Detailed Method of Preparation**
METHOD_OF_PREPARATION_PROMPT = Template("""
### **Step-by-Step Manufacturing Process for $product_name ($quanOfMed)**

The following table provides a structured breakdown of the preparation steps:

<table>
    <tr>
        <th>Step</th>
        <th>Description</th>
        <th>Equipment</th>
        <th>Time</th>
        <th>Critical Observations</th>
    </tr>
    <tr>
        <td>1</td>
        <td>Weigh all active and inactive ingredients precisely.</td>
        <td>Weighing Scale</td>
        <td>5 min</td>
        <td>Ensure accuracy in measurement.</td>
    </tr>
    <tr>
        <td>2</td>
        <td>Mix ingredients uniformly.</td>
        <td>High-speed Mixer</td>
        <td>15 min</td>
        <td>Ensure homogenous blending.</td>
    </tr>
    <tr>
        <td>3</td>
        <td>Granulation process.</td>
        <td>Fluid Bed Granulator</td>
        <td>30 min</td>
        <td>Check granule size and moisture content.</td>
    </tr>
</table>
""")

# 📌 **Detailed Characterization Process (Step-by-Step for Each Test)**
CHARACTERIZATION_PROMPT = Template("""
### **Step-by-Step Characterization for $product_name**

#### **1️⃣ Particle Size Analysis**
<table>
    <tr>
        <th>Step</th>
        <th>Description</th>
        <th>Equipment</th>
        <th>Time</th>
        <th>Critical Observations</th>
    </tr>
    <tr>
        <td>1</td>
        <td>Load sample into particle analyzer.</td>
        <td>Laser Diffraction Analyzer</td>
        <td>10 min</td>
        <td>Avoid sample contamination.</td>
    </tr>
</table>

#### **2️⃣ Surface Morphology (SEM)**
<table>
    <tr>
        <th>Step</th>
        <th>Description</th>
        <th>Equipment</th>
        <th>Time</th>
        <th>Critical Observations</th>
    </tr>
    <tr>
        <td>1</td>
        <td>Place sample on SEM stub.</td>
        <td>Scanning Electron Microscope</td>
        <td>15 min</td>
        <td>Ensure uniform distribution.</td>
    </tr>
</table>
""")

# 📌 **Formulation & Testing**
FORMULATION_TESTING_PROMPT = Template("""
### **Formulation & Testing Report for $product_name ($quanOfMed)**

#### **Formulation Process**
<table>
    <tr>
        <th>Ingredient</th>
        <th>Quantity per Unit</th>
        <th>Total Quantity</th>
        <th>Purpose</th>
    </tr>
    <tr>
        <td>API (Active Ingredient)</td>
        <td>$powerOfDrug</td>
        <td>Calculated</td>
        <td>Therapeutic effect</td>
    </tr>
</table>

#### **Testing & Quality Control**
<table>
    <tr>
        <th>Test Name</th>
        <th>Procedure</th>
        <th>Equipment Used</th>
        <th>Acceptance Criteria</th>
    </tr>
    <tr>
        <td>Disintegration Test</td>
        <td>Check breakdown in water</td>
        <td>USP Disintegration Tester</td>
        <td>< 15 min</td>
    </tr>
</table>
""")

# 📌 **FTIR Spectrum Analysis**
FTIR_PROMPT = Template("""
### **FTIR Spectrum Analysis for $product_name**

<table>
    <tr>
        <th>Wavenumber (cm⁻¹)</th>
        <th>Functional Group</th>
        <th>Peak Description</th>
        <th>Significance</th>
    </tr>
    <tr>
        <td>3200-3500</td>
        <td>O-H Stretch</td>
        <td>Strong, Broad</td>
        <td>Indicates presence of alcohol</td>
    </tr>
</table>
""")

# 📌 **Dissolution & Stability Studies**
DISSOLUTION_STABILITY_PROMPT = Template("""
### **Dissolution & Stability Study for $product_name**

#### **Dissolution Study**
<table>
    <tr>
        <th>Condition</th>
        <th>pH</th>
        <th>Time Point</th>
        <th>Result</th>
    </tr>
    <tr>
        <td>Simulated Gastric Fluid</td>
        <td>1.2</td>
        <td>30 min</td>
        <td>80% release</td>
    </tr>
</table>

#### **Stability Study**
<table>
    <tr>
        <th>Storage Condition</th>
        <th>Duration</th>
        <th>Observed Changes</th>
    </tr>
    <tr>
        <td>25°C, 60% RH</td>
        <td>3 months</td>
        <td>No significant change</td>
    </tr>
</table>
""")

# 📌 **GPT Prompt Selection**
def getPromptForOptions(options):
    if options['typeOfInfo'] == "METHOD OF PREPARATION":
        return METHOD_OF_PREPARATION_PROMPT.substitute(options)
    elif options['typeOfInfo'] == "CHARACTERIZATION/EVALUATION":
        return CHARACTERIZATION_PROMPT.substitute(options)
    elif options['typeOfInfo'] == "Both of above":
        return FORMULATION_TESTING_PROMPT.substitute(options)
    elif options['typeOfInfo'] == "CHECK RESULTS":
        return DISSOLUTION_STABILITY_PROMPT.substitute(options)
    elif options['typeOfInfo'] == "FTIR ANALYSIS":
        return FTIR_PROMPT.substitute(options)
    return ""
