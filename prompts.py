from string import Template

# Template for Table Style
TABLE_STYLE = """
<table border="1" style="border-collapse: collapse; text-align: center; width: 100%;">
    <thead>
        <tr>
            <th>Parameter</th>
            <th>Description</th>
            <th>Details</th>
        </tr>
    </thead>
    <tbody>
        $content
    </tbody>
</table>
"""

# METHOD OF PREPARATION PROMPT
METHOD_OF_PREPARATION_PROMPT = Template("""
Provide a detailed comprehensive well-structured table of method for preparing $product_name $quanOfMed, each containing $powerOfDrug of the active ingredient, based on $jurisdiction standards. The information should include: 
1. A list of all materials required, including the active pharmaceutical ingredient (API) and excipients, with their specific quantities as well as for $quanOfMed (e.g., if there are $quanOfMed tablets, calculate quantities of API and excipients accordingly).
2. A table detailing the purpose of each material used in the formulation.
3. A step-by-step table for preparation instructions, covering mixing, granulation, drying, lubrication, and compression.
4. ONLY output in HTML table format.
""")

# CHARACTERIZATION/EVALUATION PROMPT
CHARACTARIZATION_EVALUATION_PROMPT = Template("""
Provide a detailed comprehensive well-structured table of characterization of $quanOfMed of $product_name, each containing $powerOfDrug of the active ingredient, based on $jurisdiction standards. The response should include:
1. Characterization parameters for the tablets, including:
   - Physical characteristics
   - Identification tests (e.g., IR, UV, HPLC)
   - Weight variation
   - Hardness test
   - Friability test
   - Disintegration test
   - Dissolution test
   - Assay of $product_name content
2. A step-by-step table of instructions for each test (including methods, apparatus, SOPs, and acceptable limits).
3. ONLY output in HTML table format.
""")

# COMBINED PROMPT (Formulation + Testing)
COMBINED_PROMPT = Template("""
Provide a comprehensive well-structured table of guide for the formulation and testing of $product_name $quanOfMed, each containing $powerOfDrug of the active ingredient, based on $jurisdiction standards. The response should cover:
1. **Formulation Process**:
   - Materials and Quantities
   - Detailed Instructions for Preparation
2. **Testing and Analysis**:
   - Testing Criteria
   - Standards and Expectations
3. ONLY output in HTML table format.
""")

# CHECK RESULTS PROMPT
CHECK_RESULTS_PROMPT = Template("""
Compare the following evaluation results of $product_name for quantity $quanOfMed with the $jurisdiction standards:
$content
ONLY output in HTML table format.
""")

# FTIR ANALYSIS PROMPT
FTIR_PROMPT = Template("""
Provide a detailed FTIR spectrum analysis for $product_name. Ensure the response is formatted strictly in HTML table format as shown below:
<table border="1" style="border-collapse: collapse; text-align: center; width: 100%;">
  <tr>
    <th>Wavenumber (cm⁻¹)</th>
    <th>Functional Group</th>
    <th>Peak Description</th>
    <th>Significance in Drug Identification</th>
    <th>Potential Interferences</th>
    <th>Regulatory Considerations</th>
  </tr>
  $spectrum_data
</table>
After the table, provide the following explanations:
- How FTIR Confirms Drug Identity
- What Peak Deviations Indicate About Formulation Errors
- How to Ensure FTIR Compliance with Pharmacopeial Standards
""")

# SMILES STRUCTURE PROMPT
STRUCTURE_PROMPT = Template("""
Provide the **canonical SMILES notation** for the drug $product_name based on PubChem's database. Ensure that the SMILES code is accurate and matches PubChem's standard molecular structure for the drug. Return only the canonical SMILES code as provided by PubChem, and no other extra text. If the drug name is not valid, return only "NO DRUG FOUND".
""")

# PATHWAY PROMPT (For starting a pharmaceutical manufacturing company)
PATHWAY_PROMPT = Template("""
I am a registered pharmacist in India and I am starting a pharmaceutical manufacturing company for $product_type. Provide a comprehensive and detailed HTML table with all the following:
and manufacturing plant is to be Established in India so only give  data for manufacturing in India but as per &regulatory_authorities.
1. **All steps** required from company registration to manufacturing and marketing of $product_type.
2. **Documents required** at each stage.
3. **Regulatory bodies** involved at each stage.
4. **Eligibility criteria** for each process.
5. **Links to the application forms** for each process. If the exact link is unavailable, provide the **official webpage link** or guidance on how to access it.

The output should ONLY be in a structured HTML table format. Each row should represent a specific process in the pathway with the following columns:

- Stage: The stage in the process.
- Step: The specific action to be taken in this stage.
- Description: A detailed description of the step.
- Documents Required: List all documents necessary for the step.
- Regulatory Body: The body responsible for that step.
- Application Link: The direct application link (or official webpage link if the specific link is unavailable).

If a step requires accessing a specific document or application form, provide a direct link to the official website of $regulatory_authorities or relevant data website or pdf links.

<table border="1" style="border-collapse: collapse; text-align: center; width: 100%;">
    <thead>
        <tr>
            <th>Stage</th>
            <th>Step</th>
            <th>Description</th>
            <th>Documents Required</th>
            <th>Regulatory Body</th>
            <th>Application Link</th>
        </tr>
    </thead>
    <tbody>
        <!-- Fill table rows with each step from registration to marketing -->
    </tbody>
</table>
""")

LIST_OF_LICENSES_PROMPT = Template("""
Provide a detailed and comprehensive HTML table of licenses required for manufacturing, marketing, and distribution of $product_type under the guidelines of $regulatory_authorities. The table should include the following:

1. A list of all required licenses.
2. For each license:
   - **Application process**: A description of how to apply for each license all applicationi and manufacturing plant is to be Established in India so only give  data for manufacturing in India but as per &regulatory_authorities.
   - **Documents required**: List the exact documents necessary for the application.
   - **Eligibility criteria**: Describe who is eligible to apply for each license.
   - **Issuing authority**: The official body responsible for issuing the license.
   - **Official application link**: Provide the direct link to the official application page, or provide guidance on where to find it on the relevant official website.

The output should be in a structured HTML table format with these columns:

- License Name: The name of the required license.
- Application Process: Description of the application process for the license.
- Documents Required: List the specific documents needed.
- Eligibility Criteria: Who is eligible to apply for this license.
- Issuing Body: The official authority issuing the license.
- Application Link: The link to the official application page (or a reference to the official website).

If you can't find the exact link, provide the **best possible link** to the related **regulatory page or guidance** for obtaining the license. and do not provide any other response other than table nothing above table or below table

<table border="1" style="border-collapse: collapse; text-align: center; width: 100%;">
    <thead>
        <tr>
            <th>License Name</th>
            <th>Application Process</th>
            <th>Documents Required</th>
            <th>Eligibility Criteria</th>
            <th>Issuing Body</th>
            <th>Application Link</th>
        </tr>
    </thead>
    <tbody>
        <!-- Fill table rows with each required license details -->
    </tbody>
</table>
""")
DETAILED_INFORMATION_PROMPT = Template("""
Provide detailed and professional information on the "$detailed_information" document, which is required for manufacturing "$product_type" in India as per the guidelines of "$regulatory_authorities". 

This response should provide **only detailed information** for the "$detailed_information" document, and it should be structured in HTML format. Structure the output as follows:

1. **Document Overview**:
   - Describe the purpose and importance of the "$detailed_information" in pharmaceutical manufacturing.
   - Explain its role in compliance with regulations set by "$regulatory_authorities".

2. **Key Components**:
    - List and explain the essential sections of the "$detailed_information":
     - Facility Information
     - Personnel and Organizational Structure
     - Manufacturing Processes
     - Equipment and Facilities
     - Quality Control Systems
     - Safety and Environmental Controls
     - Documentation and Record Keeping

3. **Regulatory Compliance**:
   - Provide details on how the "$detailed_information" helps ensure compliance with industry standards, such as GMP (Good Manufacturing Practices) and other relevant regulations.
   
4. **Submission and Updates**:
   - Detail the process for submitting the "$detailed_information" to regulatory authorities.
   - Explain the frequency and conditions under which the "$detailed_information" must be updated.

5. **Application for the Document**:
   - **How to Apply**: Explain the step-by-step process for applying for the "$detailed_information" document.
   - **Application Time**: Mention how much time it typically takes to apply for and get approval for the "$detailed_information" document.
   - **Cost**: Provide the cost of applying for the "$detailed_information" (if applicable).
   - **Necessary Documents**: List the documents required to apply for the "$detailed_information" and also provide link of that documents which can give good infomration about neccessary documents.
   - **Application Link**: Provide a direct link or guidance on where to apply for the "$detailed_information". Include the link to the official application page, if available.
   - **PDF of Document**: If the "$detailed_information" has an official PDF or downloadable document, provide a link to that document.

6. **Additional Resources**:
   - Provide **<a href="URL">links</a>** to any relevant official resources or documents for reference.
   
7. **Formatting**:
   - Use **<h3>** for document title and section headers, and **<h4>** for subheaders.
   - Use **<ul>** and **<li>** to list items and explanations.
   - Use **<p>** for detailed paragraphs and explanations.
   - Include **<a href="URL">link</a>** to any external references or related documents (if applicable).
   - Avoid using Markdown or plain text formatting. 
   - **Only return valid HTML content**.

Ensure the content is **detailed**, **professional**, and focused specifically on the "$detailed_information", providing comprehensive regulatory and process-oriented information.
""")
def getPromptForOptions(options):
    report_type = options.get('report_type', '').strip().lower()
    type_of_info = options.get('typeOfInfo', '').strip().upper()
    jurisdiction = options.get(
        'jurisdiction', 
        'INDIAN PHARMACOPIEA, BRITISH PHARMACOPIEA, UNITED STATES PHARMACOPOEIA AND MARTINDALE-EXTRA PHARMACOPIEA'
    )

    # --- Handle report_type-based prompts ---
    if report_type == "pathway":
        return PATHWAY_PROMPT.substitute(
            product_type=options.get('product_type', ''),
            report_type=report_type,
            regulatory_authorities=options.get('regulatory', '')
        )

    elif report_type == "list of license":
        return LIST_OF_LICENSES_PROMPT.substitute(
            product_type=options.get('product_type', ''),
            regulatory_authorities=options.get('regulatory', '')
        )

    elif report_type == "detailed information":
        return DETAILED_INFORMATION_PROMPT.substitute(
            product_type=options.get('product_type', ''),
            report_type=report_type,
            regulatory_authorities=options.get('regulatory', ''),
            detailed_information=options.get('detailed_information', '')
        )

    # --- Handle typeOfInfo-based prompts ---
    if type_of_info == "METHOD OF PREPARATION":
        return METHOD_OF_PREPARATION_PROMPT.substitute(
            product_name=options.get('product_name', ''),
            quanOfMed=options.get('quanOfMed', ''),
            powerOfDrug=options.get('powerOfDrug', ''),
            jurisdiction=jurisdiction
        )

    elif type_of_info == "CHARACTARIZATION/EVALUATION":
        return CHARACTARIZATION_EVALUATION_PROMPT.substitute(
            product_name=options.get('product_name', ''),
            quanOfMed=options.get('quanOfMed', ''),
            powerOfDrug=options.get('powerOfDrug', ''),
            jurisdiction=jurisdiction
        )

    elif type_of_info == "BOTH OF ABOVE":
        return COMBINED_PROMPT.substitute(
            product_name=options.get('product_name', ''),
            quanOfMed=options.get('quanOfMed', ''),
            powerOfDrug=options.get('powerOfDrug', ''),
            jurisdiction=jurisdiction
        )

    elif type_of_info == "CHECK RESULTS":
        return CHECK_RESULTS_PROMPT.substitute(
            product_name=options.get('product_name', ''),
            quanOfMed=options.get('quanOfMed', ''),
            powerOfDrug=options.get('powerOfDrug', ''),
            jurisdiction=jurisdiction,
            content=options.get('resultsToCheck', '')
        )

    elif type_of_info == "FTIR ANALYSIS":
        return FTIR_PROMPT.substitute(options)

    elif type_of_info == "SMILES STRUCTURE":
        return STRUCTURE_PROMPT.substitute(
            product_name=options.get('product_name', '')
        )

    # --- If neither matched ---
    print("⚠️ Unrecognized report_type or typeOfInfo:", options)
    return "ERROR: Unrecognized report_type or typeOfInfo"