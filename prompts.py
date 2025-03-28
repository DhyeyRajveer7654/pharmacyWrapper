from string import Template

# template = Template("My name is $name and I am $age years old.")
# print(template.substitute(name="Alice", age=25))

# product_name

# quanOfMed

# powerOfDrug

# jurisdiction

TABLE_STYLE = "RETURN ONLY HTML TABLE FORMAT STRICTLY AND ONLY GIVE RESULT IN STRUTURED **TABLE** FORM."

METHOD_OF_PREPARATION_PROMPT = Template("""
Provide a detailed comprehensive well structured table of  method for preparing $product_name $quanOfMed, each containing $powerOfDrug of the active ingredient, based on $jurisdiction standards. The information should include: A list of all materials required, including the active pharmaceutical ingredient (API) and excipients, with their specific quantities as well as for $quanOfMed which means if there are for example $quanOfMed tablets then caluclate qualities of api and excipients for this $quanOfMed also in same table. The purpose of each material used in the formulation. Step-by-step comprehensive well structured table of instructions in all the details needed for the preparation process, including quantities and methods for mixing, granulation, drying, lubrication, and compression. Do not include any information related to evaluation, quality control, or testing procedures. The focus should solely be on the formulation and preparation steps in comprehensive well structured table form and only give output of prompt no other things and """ + TABLE_STYLE)
CHARACTARIZATION_EVALUATION_PROMPT = Template("""
Provide a detailed comprehensive well structured table of characterization of $quanOfMed of  $product_name, each containing $powerOfDrug of the active ingredient, based on $jurisdiction standards. The response should include: List of characterization parameters for the tablets required as per $jurisdiction to test the standards of the prepared formulation, including: Physical characteristics (e.g., size, shape, color, appearance) Identification tests for $product_name (e.g., IR, UV, HPLC) Weight variation Hardness test Friability test Disintegration test Dissolution test Assay of $product_name content (for content uniformity) Related substances (if applicable) Step-by-step comprehensive well structured table of instructions for performing each of the characterization tests, including: Specific methods or apparatus required for each test Standard operating procedures (SOPs) for each test as per $jurisdiction standards Acceptable limits and expected outcomes for each test. The response should focus solely on characterization parameters and testing procedures, excluding any information related to preparation, formulation, or quality control processes in comprehensive well structured table form and only give output of prompt no other things and """+TABLE_STYLE)
COMBINED_PROMPT = Template("""
Provide a comprehensive well structured table of guide for the formulation and testing of $product_name $quanOfMed, each containing $powerOfDrug of the active ingredient, based on $jurisdiction standards. The response should cover all necessary processes as outlined in the $jurisdiction guidelines. 1. Formulation Process: Materials and Quantities: List all the required substances, including the active ingredient and excipients, with their specific quantities for preparing $quanOfMed. For each ingredient, explain its role in the final product. (TABLE FORM) Detailed Instructions: Outline the complete procedure for preparing the tablets, including the precise quantities of materials to be used. This should cover all necessary stages of the formulation, ensuring that each stage is in accordance with the standards prescribed in the $jurisdiction. 2. Testing and Analysis: Testing Criteria: Provide a full description of the criteria for testing the final product, in line with the $jurisdiction. For each aspect of testing, describe the methods, tools, and procedures required to ensure that the final product meets the prescribed specifications. Standards and Expectations: For each test, detail the acceptable limits and expected results, as stated in the $jurisdiction guidelines, to ensure the formulation meets the required specifications for quality and efficacy. The response should focus solely on the formulation and testing processes, with no inclusion of information related to evaluation, quality control, or unrelated procedures in comprehensive well structured table formand only give output of prompt no other things and """+TABLE_STYLE)

CHECK_RESULTS_PROMPT = Template("""
Compare the following evaluation results in comprehensive well structured table form of my $powerOfDrug $product_name for quantity $quanOfMed with the $jurisdiction standards:
$resultsToCheck
Please compare these results with the $jurisdiction standards and assess whether they meet the required specifications in """+TABLE_STYLE)

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

def getPromptForOptions(options):
    jurisdiction = ""
    if options['jurisdiction'] == "COMPARE WITH ALL OF THEM":
        jurisdiction = "INDIAN PHARMACOPIEA, BRITISH PHARMACOPIEA, UNITED STATES PHARMACOPOEIA AND MARTINDALE-EXTRA PHARMACOPIEA"
    prompt_template = Template("")
    if options['typeOfInfo'] == "METHOD OF PREPARATION":
        prompt_template = METHOD_OF_PREPARATION_PROMPT
    elif options['typeOfInfo'] == "CHARACTARIZATION/EVALUATION":
        prompt_template = CHARACTARIZATION_EVALUATION_PROMPT
    elif options['typeOfInfo'] == "Both of above":
        prompt_template = COMBINED_PROMPT
    elif options['typeOfInfo'] == "CHECK RESULTS":
        prompt_template = CHECK_RESULTS_PROMPT
        final_prompt = prompt_template.substitute(product_name=options['product_name'], quanOfMed=options['quanOfMed'], powerOfDrug=options['powerOfDrug'], jurisdiction=jurisdiction, resultsToCheck=options['resultsToCheck'])
        return CHECK_RESULTS_PROMPT.substitute(options)
    elif options['typeOfInfo'] == "FTIR ANALYSIS":
        return FTIR_PROMPT.substitute(options)
    
    if options['jurisdiction'] == "COMPARE WITH ALL OF THEM":
        jurisdiction = "COMAPRE IN TABLE ALL 4 INDIAN PHARMACOPIEA, BRITISH PHARMACOPIEA, UNITED STATES PHARMACOPOEIA AND MARTINDALE-EXTRA PHARMACOPIEA"
    final_prompt = prompt_template.substitute(product_name=options['product_name'], quanOfMed=options['quanOfMed'], powerOfDrug=options['powerOfDrug'], jurisdiction=jurisdiction)
    print("Options Dictionary:", options.keys())  # Debugging output
    print(final_prompt)
    return final_prompt
