"""
URL configuration for visualization_workbench project.

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/5.2/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""

from django.contrib import admin
from django.urls import path, include

from django.contrib import admin
from django.urls import path, include

urlpatterns = [
    path('admin/', admin.site.urls),
    path('', include('home.urls')),  # Home page at root
    path('igv/', include('igv_viewer.urls')),
    path('pipeline/', include('pipeline_manager.urls')),
    path('ml_predictor/', include('ml_predictor.urls')),
    path('variant-annotation/', include('variant_annotation.urls')),
    path('pathway-enrichment/', include(('pathway_enrichment.urls', 'pathway_enrichment'), namespace='pathway_enrichment')),

]
