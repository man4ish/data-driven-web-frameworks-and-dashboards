import os
import uuid
import pandas as pd
import scanpy as sc  # For scRNA-seq analysis
import plotly.express as px
from flask import Flask, request, render_template, redirect, url_for

app = Flask(__name__)
UPLOAD_FOLDER = 'uploads'
os.makedirs(UPLOAD_FOLDER, exist_ok=True)

@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        file = request.files['file']
        if file:
            filepath = os.path.join(UPLOAD_FOLDER, f"{uuid.uuid4()}_{file.filename}")
            file.save(filepath)
            return redirect(url_for('visualize', filepath=filepath))
    return render_template('index.html')

@app.route('/visualize')
def visualize():
    filepath = request.args.get('filepath')

    # Read the uploaded CSV file
    df = pd.read_csv(filepath)

    # Ensure the dataframe contains only numeric data
    df = df.apply(pd.to_numeric, errors='coerce')  # Convert everything to numeric, coerce errors to NaN
    df = df.dropna(axis=0, how='any')  # Drop rows with NaN values
    
    # Perform single-cell RNA-seq analysis
    adata = sc.AnnData(df)  # Assuming the input data is already a gene expression matrix

    # Preprocessing steps (e.g., normalization, log transformation)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Perform PCA
    sc.tl.pca(adata, svd_solver='arpack', n_comps=2)  # Set n_comps to 2 or the number of components you need

    # Create a PCA plot
    pca_fig = px.scatter(x=adata.obsm['X_pca'][:, 0], y=adata.obsm['X_pca'][:, 1], title='PCA Plot')
    pca_plot = pca_fig.to_html(full_html=False)

    # Perform UMAP (after PCA)
    sc.pp.neighbors(adata, n_neighbors=10, use_rep='X_pca')
    sc.tl.umap(adata)
    umap_fig = px.scatter(x=adata.obsm['X_umap'][:, 0], y=adata.obsm['X_umap'][:, 1], title='UMAP Plot')
    umap_plot = umap_fig.to_html(full_html=False)

    # Perform clustering (e.g., Louvain)
    sc.tl.louvain(adata)
    adata.obs['louvain'] = adata.obs['louvain'].astype('category')

    # Plot UMAP with clusters
    cluster_fig = px.scatter(x=adata.obsm['X_umap'][:, 0], y=adata.obsm['X_umap'][:, 1], color=adata.obs['louvain'], title='UMAP with Clusters')
    cluster_plot = cluster_fig.to_html(full_html=False)

    return render_template('visualize.html', pca_plot=pca_plot, umap_plot=umap_plot, cluster_plot=cluster_plot)

if __name__ == '__main__':
    app.run(debug=True)
