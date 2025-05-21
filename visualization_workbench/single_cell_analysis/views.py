import os
import matplotlib
matplotlib.use('Agg')  # Use non-GUI backend for macOS and server compatibility

import scanpy as sc
import matplotlib.pyplot as plt
from django.shortcuts import render, redirect
from django.conf import settings
from django.contrib import messages
from django.http import HttpResponse
import pandas as pd
from .forms import UploadH5ADForm, SelectColorForm
from .forms import GeneExpressionForm

# Define upload directory for h5ad files
UPLOAD_DIR = os.path.join(settings.BASE_DIR, 'data', 'single_cell')
os.makedirs(UPLOAD_DIR, exist_ok=True)

# Define directory for UMAP plot
FIG_DIR = os.path.join(settings.BASE_DIR, 'static', 'figures')
os.makedirs(FIG_DIR, exist_ok=True)

def upload_file(request):
    if request.method == 'POST':
        form = UploadH5ADForm(request.POST, request.FILES)
        if form.is_valid():
            file = request.FILES['file']
            file_path = os.path.join(UPLOAD_DIR, file.name)
            
            with open(file_path, 'wb+') as dest:
                for chunk in file.chunks():
                    dest.write(chunk)

            request.session['h5ad_file'] = file_path
            return redirect('umap_plot')
    else:
        form = UploadH5ADForm()
    
    return render(request, 'single_cell_analysis/upload.html', {'form': form})



def umap_plot(request):
    file_path = request.session.get('h5ad_file')
    if not file_path or not os.path.exists(file_path):
        messages.error(request, "No .h5ad file found. Please upload a valid file.")
        return redirect('upload_h5ad')

    try:
        adata = sc.read_h5ad(file_path)
    except Exception as e:
        return render(request, 'single_cell_analysis/umap.html', {
            'error': f"Failed to read .h5ad file: {str(e)}"
        })

    if adata.shape[0] == 0:
        return render(request, 'single_cell_analysis/umap.html', {
            'error': 'The uploaded dataset contains no cells.'
        })

    # Minimal preprocessing
    if 'X_pca' not in adata.obsm:
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        sc.pp.pca(adata)
    if 'neighbors' not in adata.uns:
        sc.pp.neighbors(adata)
    if 'X_umap' not in adata.obsm:
        sc.tl.umap(adata)

    # Determine color options
    color_options = list(adata.obs.columns)
    default_color = color_options[0] if color_options else None

    selected_color = request.GET.get('color_by', default_color)

    # Save UMAP plot
    plot_filename = 'umap_plot.png'
    full_plot_path = os.path.join(FIG_DIR, plot_filename)

    try:
        sc.pl.umap(adata, color=selected_color, show=False, save=None)
        plt.savefig(full_plot_path)
        plt.close()
    except Exception as e:
        return render(request, 'single_cell_analysis/umap.html', {
            'error': f"Failed to generate UMAP plot: {str(e)}"
        })

    # Create form with dropdown
    form = SelectColorForm(color_choices=color_options, initial={'color_by': selected_color})

    return render(request, 'single_cell_analysis/umap.html', {
        'form': form,
        'plot_path': f'/static/figures/{plot_filename}',
        'filename': os.path.basename(file_path),
    })



def download_metadata(request):
    file_path = request.session.get('h5ad_file')
    if not file_path or not os.path.exists(file_path):
        return HttpResponse("No file found.", status=404)

    adata = sc.read_h5ad(file_path)
    df = adata.obs
    response = HttpResponse(content_type='text/csv')
    response['Content-Disposition'] = 'attachment; filename=metadata.csv'
    df.to_csv(path_or_buf=response)
    return response

def gene_umap(request):
    file_path = request.session.get('h5ad_file')
    if not file_path or not os.path.exists(file_path):
        return redirect('upload_h5ad')

    adata = sc.read_h5ad(file_path)

    gene_name = None
    img_path = None
    message = ""
    
    if request.method == 'POST':
        form = GeneExpressionForm(request.POST)
        if form.is_valid():
            gene_name = form.cleaned_data['gene']
            if gene_name in adata.var_names:
                fig_path = os.path.join(settings.MEDIA_ROOT, 'gene_umap.png')
                sc.pl.umap(adata, color=gene_name, show=False, save=False)
                plt.savefig(fig_path)
                img_path = os.path.join(settings.MEDIA_URL, 'gene_umap.png')
            else:
                message = f"Gene '{gene_name}' not found in dataset."
    else:
        form = GeneExpressionForm()

    return render(request, 'single_cell_analysis/gene_umap.html', {
        'form': form,
        'img_path': img_path,
        'gene_name': gene_name,
        'message': message
    })
