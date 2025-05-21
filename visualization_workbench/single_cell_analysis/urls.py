from django.urls import path
from . import views

urlpatterns = [
    path('', views.upload_file, name='upload_h5ad'),  # Default route → upload
    path('umap/', views.umap_plot, name='umap_plot'),
    path('download/', views.download_metadata, name='download_metadata'),
]
