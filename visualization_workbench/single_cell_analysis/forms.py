from django import forms

class UploadH5ADForm(forms.Form):
    file = forms.FileField(label='Upload .h5ad file')
