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

@shared_task
def run_rnaseq_nextflow_pipeline(job_id):
    job = Job.objects.get(id=job_id)
    job.status = 'RUNNING'
    job.save()

    import subprocess
    import os

    reads = job.reads_path
    outdir = job.out_dir
    log_path = os.path.join(outdir, f"job_{job.id}.log")
    job.log_file = log_path
    job.save()

    cmd = [
        "nextflow", "run", "/path/to/your/rnaseq_pipeline.nf",
        "--reads", reads,
        "--outdir", outdir
    ]

    try:
        with open(log_path, "w") as logf:
            subprocess.run(cmd, stdout=logf, stderr=logf, check=True)
        job.status = "SUCCESS"
        job.log += "Nextflow pipeline completed successfully.\n"
    except subprocess.CalledProcessError as e:
        job.status = "FAILED"
        job.log += f"Pipeline failed: {str(e)}\n"
    finally:
        job.save()

@shared_task
def run_wgs_pipeline(job_id):
    # Example: Run Whole Genome Sequencing (WGS) pipeline logic here
    # You can fetch the Job from DB by job_id, process input file, start pipeline, etc.
    print(f"Running WGS pipeline for job {job_id}")
    # Add your WGS pipeline execution code here
    # e.g., call subprocess to execute the pipeline script
    return f"WGS pipeline completed for job {job_id}"

@shared_task
def run_wes_pipeline(job_id):
    # Example: Run Whole Exome Sequencing (WES) pipeline logic here
    print(f"Running WES pipeline for job {job_id}")
    # Add your WES pipeline execution code here
    return f"WES pipeline completed for job {job_id}"

@shared_task
def run_methylation_pipeline(job_id):
    # Example: Run Methylation analysis pipeline logic here
    print(f"Running Methylation pipeline for job {job_id}")
    # Add your methylation pipeline execution code here
    return f"Methylation pipeline completed for job {job_id}"        