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

METHOD_OF_PREPARATION_PROMPT = Template("""
Provide a detailed method for preparing $product_name $quanOfMed, each containing $powerOfDrug of the active ingredient, based on $jurisdiction standards. The information should include: A list of all materials required, including the active pharmaceutical ingredient (API) and excipients, with their specific quantities for $quanOfMed. The purpose of each material used in the formulation. Step-by-step instructions in all the details needed for the preparation process, including quantities and methods for mixing, granulation, drying, lubrication, and compression. Do not include any information related to evaluation, quality control, or testing procedures. The focus should solely be on the formulation and preparation steps.                                        
                                        """)
CHARACTARIZATION_EVALUATION_PROMPT = Template("""
Provide a detailed characterization of $quanOfMed of  $product_name, each containing $powerOfDrug of the active ingredient, based on $jurisdiction standards. The response should include: List of characterization parameters for the tablets required as per $jurisdiction to test the standards of the prepared formulation, including: Physical characteristics (e.g., size, shape, color, appearance) Identification tests for $product_name (e.g., IR, UV, HPLC) Weight variation Hardness test Friability test Disintegration test Dissolution test Assay of $product_name content (for content uniformity) Related substances (if applicable) Step-by-step instructions for performing each of the characterization tests, including: Specific methods or apparatus required for each test Standard operating procedures (SOPs) for each test as per $jurisdiction standards Acceptable limits and expected outcomes for each test. The response should focus solely on characterization parameters and testing procedures, excluding any information related to preparation, formulation, or quality control processes.           
                                        """)
COMBINED_PROMPT = Template("""
Provide a comprehensive guide for the formulation and testing of $product_name $quanOfMed, each containing $powerOfDrug of the active ingredient, based on $jurisdiction standards. The response should cover all necessary processes as outlined in the $jurisdiction guidelines. 1. Formulation Process: Materials and Quantities: List all the required substances, including the active ingredient and excipients, with their specific quantities for preparing $quanOfMed. For each ingredient, explain its role in the final product. (TABLE FORM) Detailed Instructions: Outline the complete procedure for preparing the tablets, including the precise quantities of materials to be used. This should cover all necessary stages of the formulation, ensuring that each stage is in accordance with the standards prescribed in the $jurisdiction. 2. Testing and Analysis: Testing Criteria: Provide a full description of the criteria for testing the final product, in line with the $jurisdiction. For each aspect of testing, describe the methods, tools, and procedures required to ensure that the final product meets the prescribed specifications. Standards and Expectations: For each test, detail the acceptable limits and expected results, as stated in the $jurisdiction guidelines, to ensure the formulation meets the required specifications for quality and efficacy. The response should focus solely on the formulation and testing processes, with no inclusion of information related to evaluation, quality control, or unrelated procedures.                                        
                                        """)

CHECK_RESULTS_PROMPT = Template("""
Compare the following evaluation results of my $powerOfDrug $product_name for quantity $quanOfMed with the $jurisdiction standards:
$resultsToCheck
Please compare these results with the $jurisdiction standards and assess whether they meet the required specifications.""")

STRUCTURE_PROMPT = Template("""
Provide the **canonical SMILES notation** for the drug $product_name based on PubChem's database. Ensure that the SMILES code is accurate and matches PubChem's standard molecular structure for the drug. Return only the canonical SMILES code as provided by PubChem, and no other extra text. If the drug name is not valid, return only "NO DRUG FOUND".
""")

# 📌 **GPT Prompt Selection**
def getPromptForOptions(options):
    if options['typeOfInfo'] == "METHOD OF PREPARATION":
        return METHOD_OF_PREPARATION_PROMPT.substitute(options)
    elif options['typeOfInfo'] == "Both of above":
        return COMBINED_PROMPT.substitute(options)
    elif options['typeOfInfo'] == "CHECK RESULTS":
        return CHECK_RESULTS_PROMPT.substitute(options)
    elif options['typeOfInfo'] == "DISSOLUTION & STABILITY" or options['typeOfInfo'] == "CHARACTARIZATION/EVALUATION":
        return DISSOLUTION_STABILITY_PROMPT.substitute(options)
    elif options['typeOfInfo'] == "FTIR ANALYSIS":
        return FTIR_PROMPT.substitute(options)
    return ""
