from string import Template

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

def getPromptForOptions(options):
    if options['typeOfInfo'] == "METHOD OF PREPARATION":
        return METHOD_OF_PREPARATION_PROMPT.substitute(options)
    elif options['typeOfInfo'] == "CHARACTARIZATION/EVALUATION":
        return CHARACTERIZATION_EVALUATION_PROMPT.substitute(options)
    elif options['typeOfInfo'] == "Both of above":
        return COMBINED_PROMPT.substitute(options)
    elif options['typeOfInfo'] == "CHECK RESULTS":
        return CHECK_RESULTS_PROMPT.substitute(options)
    return ""
