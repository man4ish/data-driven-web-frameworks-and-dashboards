import os
from celery import Celery

os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'pipeline_manager.settings')

app = Celery('pipeline_manager')
app.config_from_object('django.conf:settings', namespace='CELERY')
app.autodiscover_tasks()
