import os
from celery import Celery

os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'visualization_workbench.settings')

app = Celery('visualization_workbench')
app.config_from_object('django.conf:settings', namespace='CELERY')
app.autodiscover_tasks()