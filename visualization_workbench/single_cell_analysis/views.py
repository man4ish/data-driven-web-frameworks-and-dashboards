import os
import matplotlib
matplotlib.use('Agg')  # Use non-GUI backend for macOS and server compatibility

import scanpy as sc
import matplotlib.pyplot as plt
from django.shortcuts import render, redirect
from django.conf import settings
from django.contrib import messages
from .forms import UploadH5ADForm
from django.http import HttpResponse
import pandas as pd
from .forms import SelectColorForm

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

    # Minimal preprocessing
    if 'X_pca' not in adata.obsm:
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        sc.pp.pca(adata)
    if 'neighbors' not in adata.uns:
        sc.pp.neighbors(adata)
    if 'X_umap' not in adata.obsm:
        sc.tl.umap(adata)

    # Get all color options
    color_options = list(adata.obs.columns)
    color = request.GET.get('color_by') or (color_options[0] if color_options else None)

    # Plot
    plot_filename = 'umap_plot.png'
    full_plot_path = os.path.join(FIG_DIR, plot_filename)

    try:
        sc.pl.umap(adata, color=color, show=False, save=None)
        plt.savefig(full_plot_path)
        plt.close()
    except Exception as e:
        return render(request, 'single_cell_analysis/umap.html', {
            'error': f"Failed to generate UMAP plot: {str(e)}"
        })

    # Add dropdown form
    form = SelectColorForm(initial={'color_by': color}, color_choices=color_options)

    return render(request, 'single_cell_analysis/umap.html', {
        'plot_path': f'/static/figures/{plot_filename}',
        'file_name': os.path.basename(file_path),
        'color': color,
        'form': form
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