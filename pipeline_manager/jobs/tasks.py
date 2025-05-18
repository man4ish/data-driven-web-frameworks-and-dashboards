from celery import shared_task
import time
from .models import Job
import subprocess
import os

# Define paths to your pipeline scripts
PIPELINE_SCRIPTS = {
    "rnaseq": "/Users/manishkumar/Desktop/llm-bio-webapps/pipelines/rnaseq_pipeline.nf",
    "wgs": "/Users/manishkumar/Desktop/llm-bio-webapps/pipelines/wgs_pipeline.nf",
    "wes": "/Users/manishkumar/Desktop/llm-bio-webapps/pipelines/wes_pipeline.nf",
    "methylation": "/Users/manishkumar/Desktop/llm-bio-webapps/pipelines/methylation_pipeline.nf"
}

NEXTFLOW_PATH = "/usr/local/bin/nextflow"
BASE_OUTPUT_DIR = "/Users/manishkumar/Desktop/llm-bio-webapps/pipeline_manager/media/job_outputs"

# Dummy task for testing
@shared_task(bind=True)
def run_dummy_pipeline(self, job_id):
    job = Job.objects.get(id=job_id)
    job.status = 'RUNNING'
    job.save()

    for i in range(5):
        job.log += f"Step {i+1}/5: running...\n"
        job.save()
        time.sleep(2)

    job.status = 'SUCCESS'
    job.log += "Pipeline finished successfully.\n"
    job.save()

# Create output directory
def setup_output_directory(job, pipeline_name):
    outdir = job.out_dir.strip() if job.out_dir else ''
    if not outdir:
        outdir = os.path.join(BASE_OUTPUT_DIR, f"{pipeline_name}_job_{job.id}")
        job.out_dir = outdir
        job.save()

    job.log += f"Output directory before creation: '{outdir}'\n"

    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)
        job.log += f"Output directory created: '{outdir}'\n"
    else:
        job.log += f"Output directory already exists: '{outdir}'\n"

    job.save()
    return outdir

# Execute pipeline command using subprocess
def run_pipeline_command(job, cmd):
    job.log += f"Running command: {' '.join(cmd)}\n"
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

# Main pipeline execution function
def run_nextflow_pipeline(job_id, pipeline_key):
    job = Job.objects.get(id=job_id)
    job.status = f"RUNNING_{pipeline_key.upper()}"
    job.log += "Starting run_nextflow_pipeline()\n"
    job.save()

    pipeline_script = PIPELINE_SCRIPTS.get(pipeline_key)
    if not pipeline_script:
        job.log += f"ERROR: Pipeline key not found: '{pipeline_key}'\n"
        job.status = "FAILED"
        job.save()
        return
    if not os.path.exists(pipeline_script):
        job.log += f"ERROR: Pipeline script does not exist: '{pipeline_script}'\n"
        job.status = "FAILED"
        job.save()
        return

    # Ensure input file exists
    if not job.input_file or not os.path.exists(job.input_file.path):
        job.log += f"ERROR: Input file not found: {job.input_file.path if job.input_file else 'None'}\n"
        job.status = "FAILED"
        job.save()
        return

    reads = job.input_file.path
    job.log += f"Input file path: {reads}\n"

    # Create output dir
    outdir = setup_output_directory(job, pipeline_key)
    job.log += "Output directory set up complete.\n"

    # Prepare and run nextflow command
    cmd = [NEXTFLOW_PATH, "run", pipeline_script, "--reads", reads, "--outdir", outdir]
    job.log += "Prepared pipeline command. Launching...\n"
    job.save()

    run_pipeline_command(job, cmd)

# Pipeline-specific Celery tasks
@shared_task
def run_rnaseq_nextflow_pipeline(job_id):
    run_nextflow_pipeline(job_id, "rnaseq")

@shared_task
def run_wgs_pipeline(job_id):
    run_nextflow_pipeline(job_id, "wgs")

@shared_task
def run_wes_pipeline(job_id):
    run_nextflow_pipeline(job_id, "wes")

@shared_task
def run_methylation_pipeline(job_id):
    run_nextflow_pipeline(job_id, "methylation")
