from flask import request, render_template, redirect, url_for
import os
from werkzeug.utils import secure_filename

UPLOAD_FOLDER = 'uploads'

def home():
    if request.method == 'POST':
        uploaded_file = request.files.get('file')
        selected_genes = request.form.get('genes')  # Optional: From form input

        if not uploaded_file:
            return "No file uploaded.", 400

        filename = secure_filename(uploaded_file.filename)
        filepath = os.path.join(UPLOAD_FOLDER, filename)
        os.makedirs(UPLOAD_FOLDER, exist_ok=True)
        uploaded_file.save(filepath)
        print(filepath)
        # Redirect to /visualize with filepath and genes
        return redirect(url_for('visualize', filepath=filepath, genes=selected_genes))
        
    return render_template('home.html')
