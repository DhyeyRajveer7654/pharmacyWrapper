from string import Template
from rdkit import Chem

METHOD_OF_PREPARATION_PROMPT = Template("""
Provide a detailed method for preparing $product_name $quanOfMed, each containing $powerOfDrug of the active ingredient, based on $jurisdiction standards.
""")

CHARACTERIZATION_EVALUATION_PROMPT = Template("""
Provide a detailed characterization of $quanOfMed of $product_name, each containing $powerOfDrug of the active ingredient, based on $jurisdiction standards.
""")

COMBINED_PROMPT = Template("""
Provide a comprehensive guide for the formulation and testing of $product_name $quanOfMed, each containing $powerOfDrug of the active ingredient, based on $jurisdiction standards.
""")

CHECK_RESULTS_PROMPT = Template("""
Compare the following evaluation results of my $powerOfDrug $product_name for quantity $quanOfMed with the $jurisdiction standards:
$resultsToCheck
""")

STRUCTURE_GENERATION_PROMPT = Template("""
Generate a valid SMILES representation for $product_name.
""")

def getPromptForOptions(options):
    if options['typeOfInfo'] == "METHOD OF PREPARATION":
        return METHOD_OF_PREPARATION_PROMPT.substitute(options)
    elif options['typeOfInfo'] == "CHARACTERIZATION/EVALUATION":
        return CHARACTERIZATION_EVALUATION_PROMPT.substitute(options)
    elif options['typeOfInfo'] == "Both of above":
        return COMBINED_PROMPT.substitute(options)
    elif options['typeOfInfo'] == "CHECK RESULTS":
        return CHECK_RESULTS_PROMPT.substitute(options)

def get_smiles_for_product(product_name):
    """Generate a SMILES code based on the product name."""
    # Placeholder logic: A real implementation would use a proper lookup or ML model
    smiles_dict = {
        "Paracetamol": "CC(=O)Nc1ccc(O)cc1",
        "Aspirin": "CC(=O)Oc1ccccc1C(=O)O",
        "Caffeine": "Cn1cnc2c1c(=O)n(c(=O)n2C)C",
    }
    return smiles_dict.get(product_name, "C")  # Default to carbon if unknown
