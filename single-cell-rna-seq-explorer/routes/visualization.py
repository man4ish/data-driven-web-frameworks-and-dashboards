from flask import request, render_template
from .utils.preprocessing import load_and_preprocess
from .utils.plotting import generate_all_plots
import os

def visualize():
    filepath = request.args.get('filepath')
    genes = request.args.get('genes')

    if not filepath:
        return "Error: 'filepath' parameter is required in the URL.", 400

    selected_genes = [gene.strip() for gene in genes.split(',')] if genes else None

    results = load_and_preprocess(filepath)

    # If there's an error in preprocessing, return it with a 400 status
    if 'error' in results:
        return results['error'], 400

    # Save the output.csv file in the 'results' folder after processing
    os.makedirs('results', exist_ok=True)  # Ensure the results folder exists
    output_file_path = os.path.join('results', 'output.csv')
    results['df'].to_csv(output_file_path, index=False)  # Save dataframe as CSV

    plots = generate_all_plots(results['adata'], results['df'], selected_genes)

    return render_template(
        'visualize.html',
        pca_plot=plots['pca_plot'],
        umap_plot=plots['umap_plot'],
        ml_cluster_plot=plots['ml_cluster_plot'],
        autoencoder_umap_plot=plots['autoencoder_umap_plot'],
        heatmap_plot=plots.get('heatmap_plot')
    )
