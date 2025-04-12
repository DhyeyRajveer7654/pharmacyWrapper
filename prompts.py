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
You are to provide ONLY a clean HTML table showing the full pathway to starting a $product_type manufacturing company under $regulatory_authorities. Do NOT include explanations, introductions, disclaimers, or text outside the table. No apologies, no extra context.

<table border="1" style="border-collapse: collapse; text-align: center; width: 100%;">
    <thead>
        <tr>
            <th>Stage</th>
            <th>Step</th>
            <th>Description</th>
            <th>Required Documents</th>
            <th>Issuing Body</th>
            <th>Official Link</th>
        </tr>
    </thead>
    <tbody>
        <!-- Fill table rows with each step from registration to marketing -->
    </tbody>
</table>
""")


# LIST OF LICENSES PROMPT (To be added)
LIST_OF_LICENSES_PROMPT = Template("""
You are to provide ONLY a clean HTML table, without any explanations, introductions, or disclaimers. Do NOT provide any text outside the table. Do not say "I cannot provide...", "I'm sorry...", or any statements before or after.

<table border="1" style="border-collapse: collapse; text-align: center; width: 100%;">
    <thead>
        <tr>
            <th>License Name</th>
            <th>Application Process</th>
            <th>Documents Required</th>
            <th>Eligibility Criteria</th>
            <th>Issuing Authority</th>
            <th>Application Link</th>
        </tr>
    </thead>
    <tbody>
        <!-- Fill table rows with licenses for $product_type under $regulatory_authorities -->
    </tbody>
</table>
""")



def getPromptForOptions(options):
    report_type = options.get('report_type', '').strip().lower()

    if report_type == "pathway":
        return PATHWAY_PROMPT.substitute(
            product_type=options.get('product_type', 'default_product'),
            report_type=options.get('report_type', 'default_report'),
            regulatory_authorities=options.get('regulatory', 'default_regulatory')
        )

    elif report_type == "list of license":
        return LIST_OF_LICENSES_PROMPT.substitute(
            product_type=options.get('product_type', 'default_product'),
            regulatory_authorities=options.get('regulatory', 'default_regulatory')
        )

    # Rest of typeOfInfo handling remains unchanged...

    # Fallback for unrecognized input
    print("⚠️ Unrecognized report_type or typeOfInfo:", options)
    return "ERROR: Unrecognized report_type or typeOfInfo"


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