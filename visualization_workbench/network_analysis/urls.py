# network_analysis/urls.py
from django.urls import path
from . import views

urlpatterns = [
    path('', views.network_home, name='network_home'),
    path('api/network/', views.fetch_network_data, name='fetch_network_data'),
    path('api/modules/', views.fetch_modules, name='fetch_modules'),
]