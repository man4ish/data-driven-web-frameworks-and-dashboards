# Nextflow Pipeline Manager with Django and Celery

This project provides a Django-based web interface to manage and run bioinformatics pipelines using Nextflow. The pipelines supported include RNA-seq, WGS, WES, and methylation analyses.

---

## Features

- Upload input files (e.g., FASTQ, VCF)
- Choose a pipeline type (rnaseq, wgs, wes, methylation)
- Run pipelines asynchronously using Celery workers
- Track job status and logs via the Django admin or user interface
- Store output files and logs per job

---

## Requirements

- Python 3.8+
- Django 3.0+
- Celery 5+
- Nextflow installed and accessible via command line
- RabbitMQ or Redis for Celery broker
- Unix-like OS recommended (macOS/Linux)

---

## Installation

1. Clone the repository:

   ```bash
   git clone https://github.com/yourusername/nextflow-pipeline-manager.git
   cd nextflow-pipeline-manager

- Install Python dependencies:

```bash
pip install -r requirements.txt
```
- Ensure Nextflow is installed and executable:

```bash
nextflow -v
```

If Nextflow is not installed, install it from https://nextflow.io/

- Setup Celery broker (RabbitMQ or Redis) and start the broker service.

- Apply Django migrations and create a superuser:

``` bash
python manage.py migrate
python manage.py createsuperuser
```

## Configuration
Pipeline Scripts
Make sure your pipeline scripts are located under:
```
/Users/manishkumar/Desktop/llm-bio-webapps/pipelines/
```
Example scripts:

```
rnaseq_pipeline.nf

wgs_pipeline.nf

wes_pipeline.nf

methylation_pipeline.nf
```

Update the PIPELINE_SCRIPTS dictionary in your Django app's tasks or settings file accordingly:

```
PIPELINE_SCRIPTS = {
    "rnaseq": "/Users/manishkumar/Desktop/llm-bio-webapps/pipelines/rnaseq_pipeline.nf",
    "wgs": "/Users/manishkumar/Desktop/llm-bio-webapps/pipelines/wgs_pipeline.nf",
    "wes": "/Users/manishkumar/Desktop/llm-bio-webapps/pipelines/wes_pipeline.nf",
    "methylation": "/Users/manishkumar/Desktop/llm-bio-webapps/pipelines/methylation_pipeline.nf"
}
```

## Nextflow Path
Update the Nextflow executable path if different:

```
NEXTFLOW_PATH = "/usr/local/bin/nextflow"
```

### Usage
Start Django server:

```bash
python manage.py runserver
```
Start Celery worker in a separate terminal:

```bash
celery -A your_project_name worker --loglevel=info
```

- Upload files and create jobs via the Django interface.

- Trigger pipelines using the appropriate Celery tasks (run_rnaseq_nextflow_pipeline, run_wgs_pipeline, etc.).

- Monitor job logs and statuses via the Django admin or API.

## Debugging & Troubleshooting
### Common Issues
1. Pipeline script not found
```
ERROR: Pipeline script does not exist: '/Users/manishkumar/Desktop/llm-bio-webapps/pipelines/wes_pipeline.nf'
```

- Make sure the .nf script file exists in the specified directory.

- Verify the path in the PIPELINE_SCRIPTS dictionary matches the actual file location.

- Use absolute paths to avoid confusion.

- Check permissions to ensure the Django process can read the file.

2. Job status stuck at PENDING or no logs generated

- Confirm that Celery workers are running and connected to the broker.

- Check Celery logs for task receipt and execution errors.

- Verify that the Celery task updates the Job model's status and log fields properly.

- Ensure the pipeline command runs without immediate failure — try running the Nextflow command manually with the same parameters to debug.

- Make sure output directories are writable and correctly configured.

3. Output directory not created or empty
- Confirm BASE_OUTPUT_DIR is set to a writable directory.

- Check if the output folder exists on disk after job runs.

- Validate that Nextflow pipeline actually produces outputs at the specified --outdir.

- Review pipeline scripts to confirm output paths.

4. Nextflow command fails silently or returns error
- Enable verbose logging in Nextflow using -trace nextflow or -with-report options.

- Capture full stdout and stderr in your Celery task logs.

- Run Nextflow command manually in terminal to reproduce and debug.

- Check Nextflow version compatibility.

## Viewing Job Logs in Django Shell
You can inspect logs directly from Django shell:

```
from jobs.models import Job
job = Job.objects.get(id=JOB_ID)  # Replace JOB_ID with your job id
print(job.log)
```

This log contains real-time debug information including:

- Commands run

- Output directory creation

- Pipeline stdout/stderr

- Errors and status updates