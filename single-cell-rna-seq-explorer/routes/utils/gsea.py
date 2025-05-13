import gseapy as gp
import matplotlib.pyplot as plt
import pandas as pd

def run_gsea(ranked_genes, gene_sets='MSigDB', output_dir='gsea_results'):
    """
    Perform Gene Set Enrichment Analysis (GSEA).

    ranked_genes: DataFrame or Series with genes and their associated statistic (e.g., fold change).
    gene_sets: Gene set collection (e.g., 'MSigDB', or path to a custom gene set).
    output_dir: Directory to save GSEA results.

    Returns: GSEA results object.
    """
    # Perform the GSEA
    gsea = gp.prerank(rnk=ranked_genes, gene_sets=gene_sets, outdir=output_dir)
    return gsea

def prepare_ranked_list(df):
    """
    Prepare ranked list of genes based on fold change or other statistics.

    df: DataFrame with gene expression data.
    """
    ranked_genes = df[['gene', 'fold_change']].sort_values(by='fold_change', ascending=False)
    ranked_genes.set_index('gene', inplace=True)
    return ranked_genes

import matplotlib.pyplot as plt

def plot_gsea_results(gsea_results):
    # Assuming gsea_results is a pandas DataFrame with relevant columns
    plt.figure(figsize=(10, 6))
    plt.plot(gsea_results['Rank'], gsea_results['EnrichmentScore'])
    plt.title('GSEA Enrichment Plot')
    plt.xlabel('Rank')
    plt.ylabel('Enrichment Score')

    # Save the plot as a PNG file in the static folder
    plt.savefig('static/gsea_plot.png')
    plt.close()

