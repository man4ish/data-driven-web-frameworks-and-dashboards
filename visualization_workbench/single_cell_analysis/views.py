import os
import scanpy as sc
import matplotlib.pyplot as plt
from django.shortcuts import render, redirect
from django.conf import settings
from .forms import UploadH5ADForm

UPLOAD_DIR = os.path.join(settings.BASE_DIR, 'data', 'single_cell')
os.makedirs(UPLOAD_DIR, exist_ok=True)

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
        return redirect('upload_h5ad')

    adata = sc.read_h5ad(file_path)
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)

    sc.pl.umap(adata, color='louvain' if 'louvain' in adata.obs else adata.obs.columns[0],
               save='_plot.png', show=False)

    img_path = os.path.join(settings.BASE_DIR, 'figures', 'umap_plot.png')
    os.makedirs(os.path.dirname(img_path), exist_ok=True)
    os.rename('figures/umap_plot.png', img_path)

    return render(request, 'single_cell_analysis/umap.html', {'plot_path': '/static/figures/umap_plot.png'})
