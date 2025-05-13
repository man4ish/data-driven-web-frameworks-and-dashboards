from flask import render_template, request
import pandas as pd
from .utils.gsea import run_gsea, plot_gsea_results, prepare_ranked_list

def gsea():
    if request.method == 'POST':
        file = request.files['file']
        gene_set = request.form['gene_set']  # Gene set chosen by the user
        
        # Load data and preprocess it
        df = pd.read_csv(file)
        
        # Prepare ranked gene list based on fold change
        ranked_genes = prepare_ranked_list(df)
        
        # Run GSEA
        gsea_results = run_gsea(ranked_genes, gene_sets=gene_set)
        
        # Generate plots
        plot_gsea_results(gsea_results)
        
        # Return results to be rendered in the template
        results = gsea_results[['GeneSet', 'EnrichmentScore', 'PValue']].values.tolist()
        return render_template('gsea_results.html', results=results)

    return render_template('upload_gsea.html')
