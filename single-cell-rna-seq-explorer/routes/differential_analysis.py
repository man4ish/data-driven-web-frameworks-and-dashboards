import pandas as pd
import numpy as np
from flask import request, render_template

def perform_differential_analysis(file, genes=None):
    """
    Perform differential expression analysis on the uploaded file.
    
    :param file: The uploaded file (CSV/Excel)
    :param genes: A comma-separated string of genes to filter (optional)
    :return: A list of top differentially expressed genes with log-fold changes
    """
    # Load the uploaded CSV or Excel file (you can modify this as per your needs)
    data = pd.read_csv(file)

    # Example: Columns in the file - ['gene', 'condition_1', 'condition_2']
    if genes:
        # If genes are provided, filter the data
        genes = genes.split(',')  # User input as comma-separated list
        data = data[data['gene'].isin(genes)]  # Filter for specified genes
    
    # Compute the log fold change (example differential analysis)
    data['log_fold_change'] = np.log2(data['condition_2'] / data['condition_1'])

    # Sort by log fold change and return the top 5 differentially expressed genes
    top_genes = data.sort_values(by='log_fold_change', ascending=False).head(5)

    # Convert results to a dictionary for rendering
    result = top_genes[['gene', 'log_fold_change']].to_dict(orient='records')
    
    return result

def differential_analysis():
    """
    Route handler for differential expression analysis.
    
    :return: Rendered template with analysis results
    """
    if request.method == 'POST':
        # Retrieve file and genes from the form
        file = request.files['file']
        genes = request.form.get('genes')

        # Perform differential analysis with the uploaded file and gene input
        results = perform_differential_analysis(file, genes)

        # Pass the results to the differential analysis template
        return render_template('differential.html', results=results)
    
    return render_template('home.html')  # If GET request, show home page for file upload
