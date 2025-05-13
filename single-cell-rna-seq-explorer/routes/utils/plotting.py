import plotly.express as px
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from ..models.autoencoder import build_autoencoder
import scanpy as sc
import plotly.graph_objects as go

def generate_all_plots(adata, df, selected_genes=None):
    plots = {}

    try:
        # PCA plot
        sc.tl.pca(adata, svd_solver='arpack', n_comps=2)
        plots['pca_plot'] = px.scatter(
            x=adata.obsm['X_pca'][:, 0],
            y=adata.obsm['X_pca'][:, 1],
            title='PCA Plot'
        ).to_html(full_html=False)

        # UMAP plot
        sc.pp.neighbors(adata, n_neighbors=10, use_rep='X_pca')
        sc.tl.umap(adata)
        plots['umap_plot'] = px.scatter(
            x=adata.obsm['X_umap'][:, 0],
            y=adata.obsm['X_umap'][:, 1],
            title='UMAP Plot'
        ).to_html(full_html=False)

        # KMeans Clustering on UMAP
        X_scaled = StandardScaler().fit_transform(df)
        adata.obs['kmeans_cluster'] = KMeans(n_clusters=3).fit_predict(X_scaled)
        plots['ml_cluster_plot'] = px.scatter(
            x=adata.obsm['X_umap'][:, 0],
            y=adata.obsm['X_umap'][:, 1],
            color=adata.obs['kmeans_cluster'],
            title='UMAP with KMeans Clusters'
        ).to_html(full_html=False)

        # Autoencoder + UMAP
        reduced = build_autoencoder(X_scaled).predict(X_scaled)
        adata.obsm['X_auto'] = reduced  # Store reduced data in `X_auto`
        sc.pp.neighbors(adata, use_rep='X_auto')
        sc.tl.umap(adata)
        plots['autoencoder_umap_plot'] = px.scatter(
            x=adata.obsm['X_umap'][:, 0],
            y=adata.obsm['X_umap'][:, 1],
            title='UMAP of Autoencoder-Reduced Data'
        ).to_html(full_html=False)

        # Gene Expression Heatmap (if selected_genes are provided)
        
        if selected_genes:
            # Filter data for selected genes
            genes_to_plot = [gene for gene in selected_genes if gene in adata.var_names]
            if genes_to_plot:
                sc.pl.heatmap(
                    adata, 
                    var_names=genes_to_plot, 
                    groupby='kmeans_cluster', 
                    show=False
                )
                # Save heatmap to a file (PNG)
                heatmap_file = 'static/heatmap.png'
                sc.pl.heatmap(
                    adata, 
                    var_names=genes_to_plot, 
                    groupby='kmeans_cluster', 
                    save=heatmap_file,
                    show=False
                )
                # Convert heatmap to Plotly figure
                heatmap_fig = go.Figure(go.Image(z=heatmap_file))
                plots['heatmap_plot'] = heatmap_fig.to_html(full_html=False)
            else:
                plots['heatmap_plot'] = None
        else:
            plots['heatmap_plot'] = None

    except Exception as e:
        return {"error": str(e)}

    return plots
