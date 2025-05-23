# airflow_api.py
import requests
import json

def trigger_annotation_dag(tool, file_paths):
    dag_id = "variant_annotation_dag"
    run_id = str(uuid.uuid4())
    data = {
        "conf": {
            "tool": tool,
            "vcf_files": file_paths
        },
        "run_id": run_id
    }
    response = requests.post(f"http://localhost:8080/api/v1/dags/{dag_id}/dagRuns", json=data,
                             auth=("airflow", "airflow"))
    response.raise_for_status()
    return run_id

def get_dag_status(run_id):
    response = requests.get(f"http://localhost:8080/api/v1/dagRuns/{run_id}", auth=("airflow", "airflow"))
    if response.status_code == 200:
        return response.json().get('state', 'unknown')
    return 'unknown'

def get_log_content(run_id):
    # Simplified for demonstration; ideally fetch actual log content
    return f"Logs for run ID: {run_id}"

