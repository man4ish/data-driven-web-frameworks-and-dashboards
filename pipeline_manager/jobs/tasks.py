from celery import shared_task
import time
from .models import Job

@shared_task(bind=True)
def run_dummy_pipeline(self, job_id):
    job = Job.objects.get(id=job_id)
    job.status = 'RUNNING'
    job.save()

    # Simulate a pipeline running
    for i in range(5):
        job.log += f"Step {i+1}/5: running...\n"
        job.save()
        time.sleep(2)  # simulate work

    job.status = 'SUCCESS'
    job.log += "Pipeline finished successfully.\n"
    job.save()
