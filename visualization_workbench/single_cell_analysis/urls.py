from django.urls import path
from . import views
from .views import gene_umap

app_name = "single_cell_analysis"

urlpatterns = [
    path('', views.upload_file, name='upload_h5ad'),  # Default route → upload
    path('umap/', views.umap_plot, name='umap_plot'),
    path('download/', views.download_metadata, name='download_metadata'),
    path('gene_umap/', gene_umap, name='gene_umap'),
]
