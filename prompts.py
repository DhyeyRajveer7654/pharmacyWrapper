from string import Template

# Styling for a centered, left-aligned table with a **transparent background**
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
        background-color: transparent;
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
        background-color: rgba(255, 255, 255, 0.1);
    }
    tr:hover {
        background-color: #007BFF;
        color: white;
    }
</style>
"""

# 📌 **Highly Detailed FTIR Spectrum Analysis (Table Format)**
FTIR_PROMPT = Template("""
Provide a **detailed FTIR spectrum analysis** for **$product_name**.

Ensure the response is a **centered HTML table** covering:
- **Wavenumber (cm⁻¹)**
- **Functional Group**
- **Peak Description**
- **Significance in Drug Identification**
- **Potential Interferences**
- **Regulatory Considerations**

The response **must be in table format only**, with **transparent background**, text **left-aligned**, and no extra text outside the table.
""")

# 📌 **GPT Prompt Selection**
def getPromptForOptions(options):
    if options['typeOfInfo'] == "FTIR ANALYSIS":
        return FTIR_PROMPT.substitute(options)
    return ""
