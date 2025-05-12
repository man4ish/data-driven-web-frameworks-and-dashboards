from django import forms

class GeneExpressionForm(forms.Form):
    file = forms.FileField(label='Upload Gene Expression File (CSV)')
