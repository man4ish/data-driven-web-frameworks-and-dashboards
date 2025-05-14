from flask import render_template, request
import pandas as pd
import matplotlib.pyplot as plt
import io
import base64
import plotly.express as px
import os

from goatools.obo_parser import GODag
from goatools.anno.genetogo_reader import Gene2GoReader
from goatools.go_enrichment import GOEnrichmentStudy
import shutil
import gzip
import urllib.request
from mygene import MyGeneInfo

def get_gene_info(gene_ids):
    mg = MyGeneInfo()
    
    # Query multiple gene IDs at once (assumes they are Entrez IDs)
    gene_info = mg.querymany(gene_ids, scopes='entrezgene', fields='symbol,name,chromosome,type_of_gene', species='human')

    # Convert to a dict: {gene_id: {symbol, name, chromosome, type_of_gene}}
    gene_metadata = {}
    for entry in gene_info:
        if not entry.get('notfound'):
            gene_metadata[entry['query']] = {
                'symbol': entry.get('symbol', ''),
                'name': entry.get('name', ''),
                'chromosome': entry.get('chromosome', ''),
                'type': entry.get('type_of_gene', '')
            }
    return gene_metadata


def load_go_terms():
    obo_file = "go-basic.obo"
    if not os.path.exists(obo_file):
        print("Downloading GO OBO file...")
        url = "http://geneontology.org/ontology/go-basic.obo"
        urllib.request.urlretrieve(url, obo_file)
    return GODag(obo_file)


# Load GO annotations using NCBI gene2go
def load_go_resources():
    obo_file = "go-basic.obo"
    gene2go_gz = "gene2go.gz"
    gene2go_file = "gene2go"

    # Download go-basic.obo if not present
    if not os.path.exists(obo_file):
        print("Downloading go-basic.obo...")
        urllib.request.urlretrieve("http://purl.obolibrary.org/obo/go/go-basic.obo", obo_file)
    
    '''
    # Download and unzip gene2go.gz if not present
    if not os.path.exists(gene2go_file):
        print("Downloading gene2go.gz...")
        urllib.request.urlretrieve("https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz", gene2go_gz)

        print("Unzipping gene2go.gz...")
        with gzip.open(gene2go_gz, 'rb') as f_in, open(gene2go_file, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    '''

    # Load GO DAG and gene2go associations
    go_dag = GODag(obo_file)
    gene2go = Gene2GoReader(gene2go_file, taxids=[9606])  # Tax ID 9606 = Homo sapiens
    ns2assoc = gene2go.get_ns2assc()

    return go_dag, ns2assoc

# GO Enrichment Analysis function
def go_enrichment():
    if request.method == 'POST':
        file = request.files['deg_file']
        deg_df = pd.read_csv(file)
        deg_genes = deg_df['Gene'].tolist()
        
        print(deg_genes)

        # Load data
        go_terms = load_go_terms()
        # Load GO DAG and annotation
        go_dag, ns2assoc = load_go_resources()

        # Combine all namespaces (BP, MF, CC) into one gene2go dict
        go_annotations = {}
        for ns, assoc in ns2assoc.items():
            for gene_id, go_ids in assoc.items():
                go_annotations.setdefault(str(gene_id), set()).update(go_ids)

        # Background genes
        background_genes = list(go_annotations.keys())

        # Run enrichment
        go_enrichment = GOEnrichmentStudy(
            background_genes,
            go_annotations,
            go_dag,
            methods=["fdr_bh"]
        )
        go_results = go_enrichment.run_study(deg_genes)

        # Format results with significant FDR < 0.05
        go_terms_enrichment = [
            (f"{res.GO} ({go_dag[res.GO].name})", res.p_fdr_bh)
            for res in go_results if res.p_fdr_bh < 0.05
        ]
        go_terms_enrichment_df = pd.DataFrame(go_terms_enrichment, columns=["GO Term", "FDR"])

        # Fetch gene metadata
        gene_metadata = get_gene_info(deg_genes)

        deg_info = pd.DataFrame([
            {
                'Entrez ID': gene_id,
                'Symbol': gene_metadata[gene_id]['symbol'],
                'Description': gene_metadata[gene_id]['name'],
                'Chromosome': gene_metadata[gene_id]['chromosome'],
                'Type': gene_metadata[gene_id]['type'],
            }
            for gene_id in deg_genes if gene_id in gene_metadata
        ])

        # Plotly bar plot
        fig = px.bar(go_terms_enrichment_df, x="FDR", y="GO Term", orientation="h",
                     title="GO Term Enrichment (FDR < 0.05)", height=600)
        graph_html = fig.to_html(full_html=False)

        # Render results in the template
        return render_template('go_enrichment_results.html',
                               graph_html=graph_html,
                               results=go_terms_enrichment_df.to_dict(orient='records'),
                               deg_metadata=deg_info.to_dict(orient='records'))

    return render_template('go_enrichment_form.html')