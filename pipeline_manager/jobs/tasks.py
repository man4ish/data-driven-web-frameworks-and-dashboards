from celery import shared_task
import time
from .models import Job
import subprocess
import os

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
    import os
    import subprocess

    job = Job.objects.get(id=job_id)
    job.status = 'RUNNING_RNASEQ'
    job.save()

    reads = job.input_file.path
    outdir = job.out_dir.strip() if job.out_dir else ''

    if not outdir:
        base_output_dir = "/Users/manishkumar/Desktop/llm-bio-webapps/pipeline_manager/media/job_outputs"
        outdir = os.path.join(base_output_dir, f"job_{job.id}")
        job.out_dir = outdir
        job.save()

    job.log += f"Output directory before creation: '{outdir}'\n"
    job.save()

    if outdir == '':
        job.log += "ERROR: output directory path is empty after processing.\n"
        job.status = "FAILED"
        job.save()
        return  # or raise an exception to halt

    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)

    job.log += f"Output directory exists or created: '{outdir}'\n"
    job.save()
    
    # Absolute path to nextflow executable - update if needed
    nextflow_path = "/usr/local/bin/nextflow"  # Update if needed
    rnaseq_pipeline_path = "/Users/manishkumar/Desktop/llm-bio-webapps/pipelines/rnaseq_pipeline.nf"  # YOUR ACTUAL PIPELINE PATH

    cmd = [
        nextflow_path, "run", rnaseq_pipeline_path,
        "--reads", reads,
        "--outdir", outdir
    ]

    job.log += f"Running command: {' '.join(cmd)}\n"
    job.save()

    try:
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
        job.log += f"STDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}\n"
        job.status = "SUCCESS"
    except subprocess.CalledProcessError as e:
        job.log += f"Pipeline failed:\n{e.stderr}\n"
        job.status = "FAILED"
    except Exception as e:
        job.log += f"Unexpected error: {str(e)}\n"
        job.status = "FAILED"
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