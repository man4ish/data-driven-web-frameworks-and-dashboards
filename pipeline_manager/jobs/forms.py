from django import forms
from .models import Job

class JobForm(forms.ModelForm):
    PIPELINE_CHOICES = [
        ('rnaseq', 'RNA-Seq'),
        ('wgs', 'Whole Genome Sequencing'),
        ('wes', 'Whole Exome Sequencing'),
        ('methylation', 'Methylation'),
        # add more pipeline types here as needed
    ]
    
    pipeline_type = forms.ChoiceField(choices=PIPELINE_CHOICES, required=True, label="Pipeline Type")
    
    class Meta:
        model = Job
        fields = ['input_file', 'pipeline_type', 'reference_genome']
