# routes/go_enrichment.py

from flask import render_template, request
import pandas as pd
from goatools import obo_parser
from goatools.go_enrichment import GOEnrichmentStudy
import matplotlib.pyplot as plt
import io
import base64
import plotly.express as px

# Load GO terms from the OBO file (we'll use a public GO OBO file here)
def load_go_terms():
    go_terms = obo_parser.GODag("http://purl.obolibrary.org/obo/go.obo")
    return go_terms

# GO Enrichment Analysis function
def go_enrichment():
    if request.method == 'POST':
        # Get uploaded DEG file (CSV with gene ids)
        file = request.files['deg_file']
        deg_df = pd.read_csv(file)  # assign to deg_df, not eg_df
        deg_genes = deg_df['Gene'].tolist()  # Use 'Gene' since your CSV has that column
        
        print("Columns in DataFrame:", deg_df.columns.tolist())
        # Load GO annotations
        go_terms = load_go_terms()

        # Example background genes (you can replace this with a list of all genes in your dataset)
        background_genes = ['Gene1', 'Gene2', 'Gene3', 'Gene4', 'Gene5']  # Replace with real gene list

        # Perform GO enrichment analysis using GOEnrichmentStudy
        go_enrichment = GOEnrichmentStudy(
            background_genes,
            go_terms,
            # Replace with actual gene ontology annotations
            go_annotations={'Gene1': ['GO:0008150'], 'Gene2': ['GO:0008150'], 'Gene3': ['GO:0003674']}, 
            methods=['fdr_bh']
        )
        go_results = go_enrichment.run_study(deg_genes)

        # Plot GO terms enrichment
        go_terms_enrichment = [(term, result.p_fdr_bh) for term, result in zip(go_terms, go_results)]
        go_terms_enrichment_df = pd.DataFrame(go_terms_enrichment, columns=['GO Term', 'FDR'])

        # Plotting FDR values
        fig = plt.figure()
        go_terms_enrichment_df.plot(kind='barh', x='GO Term', y='FDR', ax=fig.gca())
        plt.title("GO Term Enrichment (FDR)")
        
        # Convert plot to base64 to embed in HTML
        img = io.BytesIO()
        fig.savefig(img, format='png')
        img.seek(0)
        img_data = base64.b64encode(img.read()).decode('utf-8')

        # Display interactive chart with Plotly
        fig = px.bar(go_terms_enrichment_df, x="GO Term", y="FDR", title="GO Term Enrichment (FDR)")
        graph_html = fig.to_html(full_html=False)

        return render_template('go_enrichment_results.html', plot_url=img_data, graph_html=graph_html, results=go_terms_enrichment_df)

    return render_template('go_enrichment_form.html')
