from django.urls import path
from . import views

urlpatterns = [
    path('', views.home, name='pathway_enrichment_home'),  # Add a home or index view
    path('enrichment/', views.run_enrichr, name='run_enrichment'),
]
