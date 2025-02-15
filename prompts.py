from string import Template

# template = Template("My name is $name and I am $age years old.")
# print(template.substitute(name="Alice", age=25))

# product_name

# quanOfMed

# powerOfDrug

# jurisdiction

from string import Template

METHOD_OF_PREPARATION_PROMPT = Template("""
Provide a comprehensive and structured table detailing the method for preparing $product_name $quanOfMed, each containing $powerOfDrug of the active ingredient, based on $jurisdiction standards. The table should include all materials required, such as the active pharmaceutical ingredient (API) and excipients, with specific quantities, including calculations for $quanOfMed (e.g., if preparing $quanOfMed tablets, provide corresponding API and excipient amounts). Each material’s purpose in the formulation must be included. Additionally, provide a structured table with step-by-step instructions covering all preparation details, including mixing, granulation, drying, lubrication, and compression. The output should strictly be in a well-structured table format, excluding any details related to evaluation, quality control, or testing procedures. The response should consist only of the requested table without any additional text.
""")

CHARACTERIZATION_EVALUATION_PROMPT = Template("""
Generate a comprehensive and well-structured table outlining the characterization of $quanOfMed of $product_name, each containing $powerOfDrug of the active ingredient, as per $jurisdiction standards. The table should list all necessary characterization parameters for evaluating the formulation’s compliance, including physical characteristics (size, shape, color, and appearance), identification tests (IR, UV, HPLC), weight variation, hardness test, friability test, disintegration test, dissolution test, assay for $product_name content (content uniformity), and related substances (if applicable). Additionally, the table should detail the step-by-step procedures for each test, including required methods, apparatus, standard operating procedures (SOPs), and jurisdiction-based acceptance criteria. The output should be exclusively in a structured table format, without any accompanying explanatory text or unrelated information.
""")

COMBINED_PROMPT = Template("""
Generate a structured and comprehensive table that provides a detailed guide for the formulation and testing of $product_name $quanOfMed, each containing $powerOfDrug of the active ingredient, in compliance with $jurisdiction standards. The response should be divided into two sections:

1. **Formulation Process:**  
   - A table listing all required materials, including the API and excipients, along with their specific quantities for preparing $quanOfMed. Each ingredient’s role in the formulation must also be included.  
   - A detailed procedural table outlining the complete formulation steps, ensuring all necessary stages, including mixing, granulation, drying, lubrication, and compression, adhere to $jurisdiction standards.  

2. **Testing and Analysis:**  
   - A table specifying all required testing parameters, methods, and apparatus to ensure the final product meets $jurisdiction specifications.  
   - Expected results and acceptable limits for each test must be included in the table as per $jurisdiction guidelines.  

The output must strictly be in a structured table format, containing no additional text or explanations beyond the required details.
""")

CHECK_RESULTS_PROMPT = Template("""
Compare the following evaluation results for my $powerOfDrug $product_name (quantity: $quanOfMed) against the $jurisdiction standards in a structured and comprehensive table format:

$resultsToCheck

The table should clearly display the comparison between my results and the $jurisdiction standards, highlighting whether each parameter meets the required specifications. The output must strictly be in table form without any additional text or explanations.
""")

FTIR_PROMPT = Template("""
Provide a **detailed FTIR spectrum analysis** for **$product_name**.

Ensure the response is a **centered colorful blue table** covering:
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

The response **must be in colorful blue table format only**, with **white text**, text **left-aligned**, and no extra text outside the table.
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
    final_prompt = prompt_template.substitute(
    product_name=options['product_name'], 
    quanOfMed=options['quanOfMed'], 
    powerOfDrug=options['powerOfDrug'], 
    jurisdiction=jurisdiction
)
    print("Options Dictionary:", options.keys())  # Debugging output
    print(final_prompt)
    return final_prompt