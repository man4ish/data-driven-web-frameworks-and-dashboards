from django import forms

class SelectColorForm(forms.Form):
    color_by = forms.ChoiceField(label="Color by", choices=[], required=False)

    def __init__(self, *args, **kwargs):
        color_choices = kwargs.pop('color_choices', [])
        super().__init__(*args, **kwargs)
        self.fields['color_by'].choices = [(c, c) for c in color_choices]
    
class UploadH5ADForm(forms.Form):
    file = forms.FileField(label='Upload .h5ad file')