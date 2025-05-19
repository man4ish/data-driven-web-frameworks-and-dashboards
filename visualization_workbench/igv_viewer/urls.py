from django.urls import path
from . import views

urlpatterns = [
    path('', views.igv_home, name='igv_home'),
    path('genomes/', views.genome_list, name='genome_list'),
]
