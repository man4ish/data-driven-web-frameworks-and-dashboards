import pandas as pd
from django.shortcuts import render
from django.http import HttpResponse
from .forms import GeneExpressionForm
import matplotlib.pyplot as plt
import seaborn as sns
from io import BytesIO
import base64

# View to handle the file upload and process the gene expression data
def gene_expression_upload(request):
    if request.method == 'POST' and request.FILES['file']:
        file = request.FILES['file']
        data = pd.read_csv(file)  # assuming CSV for simplicity
        plot_url = generate_expression_plot(data)
        return render(request, 'webdemo/gene_expression_dashboard.html', {'plot_url': plot_url})
    
    form = GeneExpressionForm()
    return render(request, 'webdemo/gene_expression_upload.html', {'form': form})

# Function to generate plot from gene expression data
def generate_expression_plot(data):
    # Assuming the first column is Gene names, and rest are expression levels for samples
    data.set_index(data.columns[0], inplace=True)
    
    # Plot heatmap for gene expression
    plt.figure(figsize=(10, 6))
    sns.heatmap(data.T, annot=True, cmap="YlGnBu")
    plt.title("Gene Expression Heatmap")

    # Convert plot to PNG and then to base64 to embed in HTML
    img = BytesIO()
    plt.savefig(img, format='png')
    img.seek(0)
    plot_url = base64.b64encode(img.getvalue()).decode('utf-8')

    return plot_url
