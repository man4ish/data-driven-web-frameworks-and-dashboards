import csv
import requests
from django.shortcuts import render
from django.contrib import messages
from .forms import GeneListUploadForm

from django.shortcuts import render

def home(request):
    # Render a simple template, e.g. 'pathway_enrichment/home.html'
    return render(request, 'pathway_enrichment/upload.html')

def parse_gene_pvalue(file):
    genes = []
    decoded_file = file.read().decode('utf-8').splitlines()
    reader = csv.reader(decoded_file, delimiter='\t')
    for row in reader:
        if len(row) < 2:
            continue
        gene, pvalue = row[0].strip(), row[1].strip()
        try:
            pval = float(pvalue)
            genes.append((gene, pval))
        except ValueError:
            continue
    return genes


def run_enrichr(gene_list, library='KEGG_2019_Human'):
    ENRICHR_ADDLIST_URL = 'https://maayanlab.cloud/Enrichr/addList'
    ENRICHR_ENRICH_URL = 'https://maayanlab.cloud/Enrichr/enrich'

    genes_str = '\n'.join(gene_list)
    payload = {'list': (None, genes_str), 'description': (None, 'Uploaded gene list')}
    response = requests.post(ENRICHR_ADDLIST_URL, files=payload)
    if not response.ok:
        return None

    data = response.json()
    user_list_id = data.get('userListId')
    if not user_list_id:
        return None

    params = {'userListId': user_list_id, 'backgroundType': library}
    response = requests.get(ENRICHR_ENRICH_URL, params=params)
    if not response.ok:
        return None

    return response.json()


def run_gsea(genes_pvalues):
    ranked_genes = sorted(genes_pvalues, key=lambda x: x[1])  # ascending by pvalue
    # Placeholder: You can integrate a real GSEA tool here
    return {"gsea_results": f"GSEA ran on {len(ranked_genes)} genes"}


def pathway_enrichment_view(request):
    if request.method == 'POST':
        form = GeneListUploadForm(request.POST, request.FILES)
        if form.is_valid():
            file = form.cleaned_data['gene_list']
            genes_pvalues = parse_gene_pvalue(file)
            if not genes_pvalues:
                messages.error(request, "Invalid or empty gene list file.")
                return render(request, 'pathway_enrichment/upload.html', {'form': form})

            gene_symbols = [g[0] for g in genes_pvalues]

            enrichr_results = run_enrichr(gene_symbols)
            gsea_results = run_gsea(genes_pvalues)

            context = {
                'enrichr_results': enrichr_results,
                'gsea_results': gsea_results,
            }
            return render(request, 'pathway_enrichment/results.html', context)
    else:
        form = GeneListUploadForm()

    return render(request, 'pathway_enrichment/upload.html', {'form': form})
