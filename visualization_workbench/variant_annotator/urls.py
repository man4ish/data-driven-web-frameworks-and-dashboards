# urls.py
from django.urls import path
from . import views

app_name = 'variant_annotator'

urlpatterns = [
    path('upload/', views.upload_view, name='upload'),
    path('monitor/<str:run_id>/', views.monitor_view, name='monitor'),
]