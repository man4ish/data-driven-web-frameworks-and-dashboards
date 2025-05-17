from django.shortcuts import render, redirect
from django.contrib.auth.decorators import login_required
from .forms import JobForm
from .models import Job
from .tasks import run_dummy_pipeline, run_rnaseq_nextflow_pipeline, run_wgs_pipeline, run_wes_pipeline, run_methylation_pipeline

@login_required
def submit_job(request):
    if request.method == 'POST':
        form = JobForm(request.POST, request.FILES)
        if form.is_valid():
            job = form.save(commit=False)
            job.user = request.user
            job.save()
            
            # Launch pipeline based on selected type
            pipeline_type = form.cleaned_data.get('pipeline_type')
            if pipeline_type == 'rnaseq':
                run_rnaseq_nextflow_pipeline.delay(job.id)
            elif pipeline_type == 'wgs':
                run_wgs_pipeline.delay(job.id)
            elif pipeline_type == 'wes':
                run_wes_pipeline.delay(job.id)
            elif pipeline_type == 'methylation':
                run_methylation_pipeline.delay(job.id)
            else:
                run_dummy_pipeline.delay(job.id)

            return redirect('job_status', job_id=job.id)
    else:
        form = JobForm()
    return render(request, 'jobs/submit_job.html', {'form': form})

@login_required
def job_status(request, job_id):
    job = Job.objects.get(id=job_id, user=request.user)
    return render(request, 'jobs/job_status.html', {'job': job})

@login_required
def jobs_home(request):
    jobs = Job.objects.filter(user=request.user).order_by('-created_at')
    return render(request, 'jobs/home.html', {'jobs': jobs})
