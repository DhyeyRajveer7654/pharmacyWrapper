from string import Template

TABLE_STYLE = """
<style>
    table {
        width: 100%;
        border-collapse: collapse;
        background-color: white;
        color: black;
        border-radius: 10px;
        font-size: 16px;
        border: 1px solid #ddd;
    }
    th {
        background-color: #007BFF;
        color: white;
        padding: 12px;
        text-align: center;
    }
    td {
        border: 1px solid #ddd;
        padding: 10px;
        text-align: left;
    }
    tr:nth-child(even) {
        background-color: #f2f2f2;
    }
    tr:hover {
        background-color: #007BFF;
        color: white;
    }
</style>
"""

def format_as_table(data):
    table_html = TABLE_STYLE + "<table><tr>"
    headers = data[0].keys()

    for header in headers:
        table_html += f"<th>{header}</th>"
    table_html += "</tr>"

    for row in data:
        table_html += "<tr>"
        for value in row.values():
            table_html += f"<td>{value}</td>"
        table_html += "</tr>"

    table_html += "</table>"
    return table_html

METHOD_OF_PREPARATION_PROMPT = Template("""
Provide a detailed method for preparing $product_name ($quanOfMed), each containing $powerOfDrug of the active ingredient, based on $jurisdiction standards.
Ensure the response is **only** in HTML table format with no extra text.
""")

CHARACTERIZATION_EVALUATION_PROMPT = Template("""
Provide a detailed characterization of $product_name ($quanOfMed), each containing $powerOfDrug, based on $jurisdiction standards.
Ensure the response is **only** in HTML table format with no extra text.
""")

CHECK_RESULTS_PROMPT = Template("""
Compare evaluation results of $product_name ($powerOfDrug) for quantity $quanOfMed with $jurisdiction standards.
Ensure the response is **only** in HTML table format with no extra text.
$resultsToCheck
""")

def getPromptForOptions(options):
    if options['typeOfInfo'] == "METHOD OF PREPARATION":
        return METHOD_OF_PREPARATION_PROMPT.substitute(options)
    elif options['typeOfInfo'] == "CHARACTARIZATION/EVALUATION":
        return CHARACTERIZATION_EVALUATION_PROMPT.substitute(options)
    elif options['typeOfInfo'] == "CHECK RESULTS":
        return CHECK_RESULTS_PROMPT.substitute(options)
    return ""
