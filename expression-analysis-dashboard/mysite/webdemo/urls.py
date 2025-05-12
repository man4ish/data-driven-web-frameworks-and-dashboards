# webdemo/urls.py

from django.urls import path
from . import views

urlpatterns = [
    path('upload/', views.gene_expression_upload, name='upload'),
    # If you want to make the homepage the upload page
    path('', views.gene_expression_upload, name='home'),  # Map the root URL to the upload view

]
