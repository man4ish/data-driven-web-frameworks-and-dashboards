# routes/download.py
from flask import send_file
import os

def download():
    # Construct an absolute path to the results/output.csv file
    base_dir = os.path.dirname(os.path.abspath(__file__))  # Current file's directory
    file_path = os.path.join(base_dir, '..', 'results', 'output.csv')  # Navigate up and into 'results/'

    # Normalize the path
    file_path = os.path.abspath(file_path)

    # Check if file exists before sending
    if os.path.exists(file_path):
        return send_file(file_path, as_attachment=True)
    else:
        return "<p><strong>File not found:</strong> Make sure the output.csv has been generated in the 'results' folder.</p>", 404
