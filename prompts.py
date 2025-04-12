from string import Template

TABLE_STYLE = """
<table border="1" style="border-collapse: collapse; text-align: center; width: 100%;">
    <thead>
        <tr>
            <th>Material</th>
            <th>Quantity</th>
            <th>Purpose</th>
            <th>Instructions</th>
        </tr>
    </thead>
    <tbody>
"""

TABLE_END = "</tbody></table>"

METHOD_OF_PREPARATION_PROMPT = Template("""
Provide a detailed comprehensive well structured table of method for preparing $product_name $quanOfMed, each containing $powerOfDrug of the active ingredient, based on $jurisdiction standards. The information should include:
1. A list of all materials required, including the active pharmaceutical ingredient (API) and excipients, with their specific quantities, including for $quanOfMed. For example, if there are $quanOfMed tablets, calculate quantities for the API and excipients for this $quanOfMed in the table as well.
2. The purpose of each material used in the formulation.
3. Step-by-step comprehensive well-structured table of instructions, including mixing, granulation, drying, lubrication, and compression methods. Include specific quantities and procedures.
Please ensure all information is strictly formatted in the following structured HTML table format:

""" + TABLE_STYLE)

CHARACTARIZATION_EVALUATION_PROMPT = Template("""
Provide a detailed comprehensive well structured table of characterization of $quanOfMed of $product_name, each containing $powerOfDrug of the active ingredient, based on $jurisdiction standards. The response should include:
1. List of characterization parameters for the tablets required as per $jurisdiction, including:
    - Physical characteristics (size, shape, color, appearance)
    - Identification tests (IR, UV, HPLC)
    - Weight variation, hardness test, friability test, disintegration test, dissolution test
    - Assay of $product_name content, related substances (if applicable)
2. Step-by-step comprehensive table of instructions for performing each of the characterization tests, including:
    - Methods or apparatus required for each test
    - Standard operating procedures (SOPs) for each test
    - Acceptable limits and expected outcomes for each test
Ensure that the response strictly follows the HTML table format outlined:

""" + TABLE_STYLE)

COMBINED_PROMPT = Template("""
Provide a comprehensive well structured table of guide for the formulation and testing of $product_name $quanOfMed, each containing $powerOfDrug of the active ingredient, based on $jurisdiction standards. The response should cover:
1. Formulation Process:
    - Materials and Quantities: List all required substances, including the active ingredient and excipients, with specific quantities for $quanOfMed.
    - Detailed Instructions: Outline the complete preparation procedure with quantities for each material.
2. Testing and Analysis:
    - Testing Criteria: Provide a full description of testing criteria for the final product in line with $jurisdiction.
    - Standards and Expectations: Detail acceptable limits and expected results for each test.
The response should only include output in HTML table format as shown below:

""" + TABLE_STYLE)

CHECK_RESULTS_PROMPT = Template("""
Compare the following evaluation results in comprehensive well structured table form of my $powerOfDrug $product_name for quantity $quanOfMed with the $jurisdiction standards:
$resultsToCheck
Please compare these results with the $jurisdiction standards and assess whether they meet the required specifications in the following structured table format:

""" + TABLE_STYLE)

FTIR_PROMPT = Template("""
Provide a detailed FTIR spectrum analysis for $product_name.

Ensure the response is formatted strictly in HTML table format as shown below:
<table border="1" style="border-collapse: collapse; text-align: center; width: 100%;">
  <tr>
    <th>Wavenumber (cm⁻¹)</th>
    <th>Functional Group</th>
    <th>Peak Description</th>
    <th>Significance in Drug Identification</th>
    <th>Potential Interferences</th>
    <th>Regulatory Considerations</th>
  </tr>
  <tr>
    <td>3300-3500</td>
    <td>N-H Stretching</td>
    <td>Strong, broad peak</td>
    <td>Confirms presence of amide group</td>
    <td>Water, amines</td>
    <td>Should match reference spectrum</td>
  </tr>
  <tr>
    <td>3150-3450</td>
    <td>O-H Stretching</td>
    <td>Broad peak</td>
    <td>Indicates phenol presence</td>
    <td>Alcohols, moisture</td>
    <td>Must align with pharmacopeial data</td>
  </tr>
</table>
After the table, provide the following explanations in a structured, point-wise format:

How FTIR Confirms Drug Identity
What Peak Deviations Indicate About Formulation Errors
How to Ensure FTIR Compliance with Pharmacopeial Standards
""")

STRUCTURE_PROMPT = Template("""
Provide the **canonical SMILES notation** for the drug $product_name based on PubChem's database. Ensure that the SMILES code is accurate and matches PubChem's standard molecular structure for the drug. Return only the canonical SMILES code as provided by PubChem, and no other extra text. If the drug name is not valid, return only "NO DRUG FOUND".
""")

PATHWAY_PROMPT = Template("""
I am a registered pharmacist in Australia who has completed a Bachelor's degree in Pharmacy. I am now planning to start my own pharmaceutical manufacturing company, specifically focused on the production of $product_type . I want a complete, step-by-step $report_type that outlines all the processes, legal requirements, licenses, documents, and compliance steps that I need to follow to start and operate a pharmaceutical company in accordance with $regulatory_authorities guidelines.

Please provide:

1. Provide a detailed comprehensive well structured table of $report_type outlining all the steps from company registration to manufacturing and marketing of pharmaceutical $product_type.
2. A list of all required licenses and approvals at each stage, including but not limited to manufacturing license, product registration as per $regulatory_authorities guidelines.
3. Direct links to each license application form and the relevant guidelines from $regulatory_authorities .
4. A detailed breakdown of $regulatory_authorities requirements for $product_type manufacturing, including links to the official data pdfs.
5. A dedicated section that covers the $regulatory_authorities Manufacturing License application process specifically:
    - What is the process?
    - Which documents are required?
    - Who is eligible to apply?
    - What technical and quality-related systems must be in place?
6. A comprehensive table that lists all documents required to be submitted with the $regulatory_authorities Manufacturing License application, categorized appropriately (e.g., company info, premises, quality systems, personnel, validation, etc.).
7. If applicable, include any requirements or processes related to pharmacovigilance, post-market surveillance, and export licensing or certification under $regulatory_authorities.
8. The goal is to have a single, well-structured document that can guide a beginner pharmacist like me through the entire journey of starting a $product_type manufacturing company in full compliance with $regulatory_authorities standards.
""")


LIST_OF_LICENSES_PROMPT = Template("""
Provide a comprehensive list of all necessary licenses and approvals required to manufacture and distribute $product_type in accordance with $regulatory_authorities guidelines. This should include:
1. Manufacturing License
2. Product Registration
3. Distribution License
4. Export/Import License (if applicable)
5. Any additional specialized licenses as required by $regulatory_authorities
The response should be formatted in a well-structured HTML table, showing:
    - License Name
    - License Description
    - Application Process
    - Required Documentation
    - Links to Application Forms and Guidelines
""" + TABLE_STYLE)


# Define the getPromptForOptions function here to avoid circular imports
def getPromptForOptions(options):
    # Check for report_type based prompts
    if options.get('report_type') == "Pathway":
        return PATHWAY_PROMPT.substitute(
            product_type=options.get('product_type', 'default_product'),
            report_type=options.get('report_type', 'default_report'),
            regulatory_authorities=options.get('regulatory', 'default_regulatory')
        )

    elif options.get('report_type') == "List of License":
        return LIST_OF_LICENSES_PROMPT.substitute(
            product_type=options.get('product_type', 'default_product'),
            regulatory_authorities=options.get('regulatory', 'default_regulatory')
        )

    # Default: Handle typeOfInfo prompts
    jurisdiction = options.get('jurisdiction', 'INDIAN PHARMACOPIEA, BRITISH PHARMACOPIEA, UNITED STATES PHARMACOPOEIA AND MARTINDALE-EXTRA PHARMACOPIEA')

    prompt_template = Template("")

    if options.get('typeOfInfo') == "METHOD OF PREPARATION":
        prompt_template = METHOD_OF_PREPARATION_PROMPT
    elif options.get('typeOfInfo') == "CHARACTARIZATION/EVALUATION":
        prompt_template = CHARACTARIZATION_EVALUATION_PROMPT
    elif options.get('typeOfInfo') == "Both of above":
        prompt_template = COMBINED_PROMPT
    elif options.get('typeOfInfo') == "CHECK RESULTS":
        return CHECK_RESULTS_PROMPT.substitute(
            product_name=options.get('product_name', ''),
            quanOfMed=options.get('quanOfMed', ''),
            powerOfDrug=options.get('powerOfDrug', ''),
            jurisdiction=jurisdiction,
            resultsToCheck=options.get('resultsToCheck', '')
        )
    elif options.get('typeOfInfo') == "FTIR ANALYSIS":
        return FTIR_PROMPT.substitute(options)
    elif options.get('typeOfInfo') == "SMILES STRUCTURE":
        return STRUCTURE_PROMPT.substitute(product_name=options.get('product_name', ''))

    # Default return for typeOfInfo prompts
    final_prompt = prompt_template.substitute(
        product_name=options.get('product_name', ''),
        quanOfMed=options.get('quanOfMed', ''),
        powerOfDrug=options.get('powerOfDrug', ''),
        jurisdiction=jurisdiction
    )

    return final_prompt