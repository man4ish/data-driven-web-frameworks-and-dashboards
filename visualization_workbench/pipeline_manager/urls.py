app_name = 'pipeline_manager'  # <-- add this

from django.urls import path
from . import views

urlpatterns = [
    path('', views.index, name='pipeline_home'),
]
