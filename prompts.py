from string import Template

# Table Style Constant
TABLE_STYLE = "RETURN ONLY HTML TABLE FORMAT STRICTLY AND ONLY GIVE RESULT IN STRUCTURED **TABLE** FORM."

# Prompt Templates
METHOD_OF_PREPARATION_PROMPT = Template("""
Provide a detailed comprehensive well structured table of method for preparing $product_name $quanOfMed, each containing $powerOfDrug of the active ingredient, based on $jurisdiction standards. The information should include: 
- A list of all materials required, including the active pharmaceutical ingredient (API) and excipients, with their specific quantities for $quanOfMed.
- The purpose of each material used in the formulation.
- Step-by-step comprehensive well structured table of instructions for the preparation process, including quantities and methods for mixing, granulation, drying, lubrication, and compression.
Do not include any information related to evaluation, quality control, or testing procedures.
""" + TABLE_STYLE)

CHARACTERIZATION_EVALUATION_PROMPT = Template("""
Provide a detailed comprehensive well structured table of characterization of $quanOfMed of $product_name, each containing $powerOfDrug of the active ingredient, based on $jurisdiction standards. The response should include:
- List of characterization parameters: Physical characteristics, identification tests, weight variation, hardness, friability, disintegration, dissolution, assay, related substances.
- Step-by-step comprehensive well structured table for each characterization test: method/apparatus, SOPs, acceptable limits, and expected outcomes.
Exclude information unrelated to characterization or testing procedures.
""" + TABLE_STYLE)

COMBINED_PROMPT = Template("""
Provide a comprehensive well structured table of guide for the formulation and testing of $product_name $quanOfMed, each containing $powerOfDrug of the active ingredient, based on $jurisdiction standards.

1. **Formulation Process**:
- Materials and Quantities: List required substances, including API and excipients, their roles, and specific quantities for $quanOfMed.
- Detailed Instructions: Step-by-step table covering preparation stages as per $jurisdiction standards.

2. **Testing and Analysis**:
- Testing Criteria: Describe methods/tools/procedures for quality check.
- Standards and Expectations: Detail acceptable limits and results based on $jurisdiction guidelines.

Exclude any unrelated content.
""" + TABLE_STYLE)

CHECK_RESULTS_PROMPT = Template("""
Compare the following evaluation results in comprehensive well structured table form of my $powerOfDrug $product_name for quantity $quanOfMed with the $jurisdiction standards:

$resultsToCheck

Please compare these results with the $jurisdiction standards and assess whether they meet the required specifications.
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

Also provide:
- How FTIR Confirms Drug Identity
- What Peak Deviations Indicate About Formulation Errors
- How to Ensure FTIR Compliance with Pharmacopeial Standards
""")

STRUCTURE_PROMPT = Template("""
Provide the **canonical SMILES notation** for the drug $product_name based on PubChem's database. 
Return only the canonical SMILES code and no other text. 
If the drug is not found, return "NO DRUG FOUND".
""")

PATHWAY_PROMPT = Template("""
You are a regulatory expert tasked with guiding a registered pharmacist in India who has completed a Bachelor's degree in Pharmacy and now plans to start a pharmaceutical manufacturing company focused on $product_type.

Please provide a **strictly HTML table-only output** that outlines all necessary steps, licenses, documentation, and compliance requirements according to the $regulatory_authorities guidelines.

Include the following columns in the HTML table:
- Stage Number
- Process/Activity
- Description of Step
- Required Licenses or Approvals
- Responsible Authority
- Relevant Links

DO NOT include any text before or after the table. DO NOT add headings, explanations, or markdown formatting. RETURN ONLY THE TABLE STRICTLY IN HTML TABLE FORMAT WITH HEADINGS. """ + TABLE_STYLE)

# Function to generate prompt
def getPromptForOptions(options):
    print("Options received:", options)
    
    jurisdiction = options.get('jurisdiction', '')
    if jurisdiction == "COMPARE WITH ALL OF THEM":
        jurisdiction = "INDIAN PHARMACOPIEA, BRITISH PHARMACOPIEA, UNITED STATES PHARMACOPOEIA AND MARTINDALE-EXTRA PHARMACOPIEA"

    if options.get('report_type') == "Pathway":
        return PATHWAY_PROMPT.substitute(
            product_type=options.get('product_type', 'default_product'),
            report_type=options.get('report_type', 'default_report'),
            regulatory_authorities=options.get('regulatory', 'default_regulatory')
        )

    prompt_type = options.get('typeOfInfo')
    if prompt_type == "METHOD OF PREPARATION":
        return METHOD_OF_PREPARATION_PROMPT.substitute(**options, jurisdiction=jurisdiction)
    elif prompt_type == "CHARACTARIZATION/EVALUATION":
        return CHARACTERIZATION_EVALUATION_PROMPT.substitute(**options, jurisdiction=jurisdiction)
    elif prompt_type == "Both of above":
        return COMBINED_PROMPT.substitute(**options, jurisdiction=jurisdiction)
    elif prompt_type == "CHECK RESULTS":
        return CHECK_RESULTS_PROMPT.substitute(**options, jurisdiction=jurisdiction)
    elif prompt_type == "FTIR ANALYSIS":
        return FTIR_PROMPT.substitute(**options)
    elif prompt_type == "SMILES STRUCTURE":
        return STRUCTURE_PROMPT.substitute(product_name=options.get('product_name', ''))
    
    return "Invalid option selected."
