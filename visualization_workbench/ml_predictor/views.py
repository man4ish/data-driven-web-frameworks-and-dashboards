from django.shortcuts import render
from .forms import UploadDataForm
import pandas as pd
from .ml_model.predictor import make_prediction

def upload_view(request):
    if request.method == 'POST':
        form = UploadDataForm(request.POST, request.FILES)
        if form.is_valid():
            file = request.FILES['file']
            df = pd.read_csv(file)
            predictions = make_prediction(df)
            return render(request, 'ml_predictor/results.html', {
                'predictions': predictions.to_html()
            })
    else:
        form = UploadDataForm()
    return render(request, 'ml_predictor/upload.html', {'form': form})