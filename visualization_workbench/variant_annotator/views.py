# views.py
from django.shortcuts import render, redirect
from .forms import UploadForm
from .airflow_api import trigger_annotation_dag, get_dag_status, get_log_content
import os
import uuid
from django.conf import settings

def upload_view(request):
    if request.method == 'POST':
        form = UploadForm(request.POST, request.FILES)
        if form.is_valid():
            files = request.FILES.getlist('vcf_files')
            tool = form.cleaned_data['tool']
            upload_dir = os.path.join(settings.MEDIA_ROOT, str(uuid.uuid4()))
            os.makedirs(upload_dir, exist_ok=True)
            paths = []
            for file in files:
                path = os.path.join(upload_dir, file.name)
                with open(path, 'wb+') as destination:
                    for chunk in file.chunks():
                        destination.write(chunk)
                paths.append(path)

            run_id = trigger_annotation_dag(tool, paths)
            return redirect('variant_annotator:monitor', run_id=run_id)
    else:
        form = UploadForm()
    return render(request, 'variant_annotator/upload.html', {'form': form})

def monitor_view(request, run_id):
    status = get_dag_status(run_id)
    log = get_log_content(run_id)
    return render(request, 'variant_annotator/monitor.html', {
        'run_id': run_id,
        'status': status,
        'log': log
    })