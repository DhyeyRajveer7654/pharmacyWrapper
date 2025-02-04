from string import Template

# Templates for various types of pharmaceutical information
METHOD_OF_PREPARATION_PROMPT = Template("""
Provide a detailed method for preparing $product_name $quanOfMed, each containing $powerOfDrug of the active ingredient, based on $jurisdiction standards. 
The information should include: 
- A list of all materials required, including the active pharmaceutical ingredient (API) and excipients, with their specific quantities for $quanOfMed. 
- The purpose of each material used in the formulation. 
- Step-by-step instructions in all the details needed for the preparation process, including quantities and methods for mixing, granulation, drying, lubrication, and compression. 

Do not include any information related to evaluation, quality control, or testing procedures. The focus should solely be on the formulation and preparation steps.
""")

CHARACTARIZATION_EVALUATION_PROMPT = Template("""
Provide a detailed characterization of $quanOfMed of $product_name, each containing $powerOfDrug of the active ingredient, based on $jurisdiction standards. 
The response should include: 
- A list of characterization parameters for the tablets required as per $jurisdiction to test the standards of the prepared formulation, including: 
  - Physical characteristics (e.g., size, shape, color, appearance) 
  - Identification tests for $product_name (e.g., IR, UV, HPLC) 
  - Weight variation 
  - Hardness test 
  - Friability test 
  - Disintegration test 
  - Dissolution test 
  - Assay of $product_name content (for content uniformity) 
  - Related substances (if applicable) 

Step-by-step instructions for performing each of the characterization tests, including:
- Specific methods or apparatus required for each test 
- Standard operating procedures (SOPs) for each test as per $jurisdiction standards 
- Acceptable limits and expected outcomes for each test. 

The response should focus solely on characterization parameters and testing procedures, excluding any information related to preparation, formulation, or quality control processes.
""")

COMBINED_PROMPT = Template("""
Provide a comprehensive guide for the formulation and testing of $product_name $quanOfMed, each containing $powerOfDrug of the active ingredient, based on $jurisdiction standards. 
The response should cover all necessary processes as outlined in the $jurisdiction guidelines.

1. **Formulation Process**:
   - Materials and Quantities: List all the required substances, including the active ingredient and excipients, with their specific quantities for preparing $quanOfMed. 
   - For each ingredient, explain its role in the final product. (TABLE FORM)
   - Detailed Instructions: Outline the complete procedure for preparing the tablets, including the precise quantities of materials to be used. 
   - This should cover all necessary stages of the formulation, ensuring that each stage is in accordance with the standards prescribed in the $jurisdiction. 

2. **Testing and Analysis**:
   - Testing Criteria: Provide a full description of the criteria for testing the final product, in line with the $jurisdiction. 
   - For each aspect of testing, describe the methods, tools, and procedures required to ensure the final product meets the prescribed specifications. 
   - Standards and Expectations: For each test, detail the acceptable limits and expected results, as stated in the $jurisdiction guidelines, to ensure the formulation meets the required specifications for quality and efficacy. 

The response should focus solely on the formulation and testing processes, with no inclusion of information related to quality control or unrelated procedures.
""")

CHECK_RESULTS_PROMPT = Template("""
Compare the following evaluation results of my $powerOfDrug $product_name for quantity $quanOfMed with the $jurisdiction standards:
$resultsToCheck

Please compare these results with the $jurisdiction standards and assess whether they meet the required specifications.
""")

# **New prompt for retrieving the chemical structure of a drug**
CHEMICAL_STRUCTURE_PROMPT = Template("""
Provide a high-quality image (PNG, JPEG, or JPG) of the chemical structure of $product_name. 
The structure must be sourced from trusted scientific databases such as:
- **Pharmacopoeias**
- **PubChem**
- **PubMed**
- **NCBI**
- **Peer-reviewed scholarly articles** 

Additionally, include the **IUPAC name** and **molecular formula** alongside the structure if available. 
Ensure the structure matches the standard representation from these authoritative sources.
""")

# Function to ensure output is only an HTML table
def addTextToReturnOnlyHtmlTable(promptWithoutHtml):
    promptToOnlyGetHtmlTable = "\nEnsure the response is only a valid HTML table with no additional text or explanation. The table should have a black background (background-color: black) and all text should be white (color: white) with white borders and proper spacing. Format the response as clean, minimal HTML."
    return promptWithoutHtml + promptToOnlyGetHtmlTable

# Function to generate the required prompt based on user-selected options
def getPromptForOptions(options):
    jurisdiction = ""
    if options['jurisdiction'] == "COMPARE WITH ALL OF THEM":
        jurisdiction = "INDIAN PHARMACOPIEA, BRITISH PHARMACOPIEA, and UNITED STATES PHARMACOPOEIA"
    
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
        return addTextToReturnOnlyHtmlTable(final_prompt)
    
    if options['jurisdiction'] == "COMPARE WITH ALL OF THEM":
        jurisdiction = "Show different results according to all of these jurisdictions: INDIAN PHARMACOPIEA, BRITISH PHARMACOPIEA, and UNITED STATES PHARMACOPOEIA"
    
    final_prompt = prompt_template.substitute(product_name=options['product_name'], quanOfMed=options['quanOfMed'], powerOfDrug=options['powerOfDrug'], jurisdiction=jurisdiction)
    final_prompt = addTextToReturnOnlyHtmlTable(final_prompt)
    return final_prompt

# Function to get the prompt for only the chemical structure
def getPromptForChemicalStructure(product_name):
    return CHEMICAL_STRUCTURE_PROMPT.substitute(product_name=product_name)
