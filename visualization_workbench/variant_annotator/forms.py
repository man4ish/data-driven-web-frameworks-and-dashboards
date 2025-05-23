# forms.py
from django import forms

TOOL_CHOICES = [
    ("snpeff", "SnpEff"),
    ("vep", "VEP"),
    ("annovar", "ANNOVAR"),
]

class UploadForm(forms.Form):
    vcf_files = forms.FileField(widget=forms.ClearableFileInput(attrs={'multiple': True}))
    tool = forms.ChoiceField(choices=TOOL_CHOICES)