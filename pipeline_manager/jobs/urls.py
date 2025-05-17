from django.urls import path
from . import views

urlpatterns = [
    path('', views.jobs_home, name='jobs_home'),  # This handles /jobs/
    path('submit/', views.submit_job, name='submit_job'),
    path('status/<int:job_id>/', views.job_status, name='job_status'),
]
