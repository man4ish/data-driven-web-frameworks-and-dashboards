from flask import Flask, render_template, request
import subprocess
import os
import uuid

app = Flask(__name__)
UPLOAD_FOLDER = 'uploads'
os.makedirs(UPLOAD_FOLDER, exist_ok=True)

@app.route('/', methods=['GET', 'POST'])
def index():
    result = ""
    if request.method == 'POST':
        sequence = request.form.get('sequence')
        file = request.files.get('file')

        # Create a unique filename
        input_file = os.path.join(UPLOAD_FOLDER, f"input_{uuid.uuid4().hex}.fasta")
        with open(input_file, 'w') as f:
            if file and file.filename:
                f.write(file.read().decode('utf-8'))
            elif sequence:
                f.write(sequence)
            else:
                return "No sequence input."

        # Output file
        output_file = input_file.replace(".fasta", "_out.txt")

        # Run BLAST (example with blastn, adjust as needed)
        blast_command = f"blastn -query {input_file} -db mydb -out {output_file} -outfmt 7"
 
        # blast_command = f"blastn -query {input_file} -db nt -out {output_file} -outfmt 7"
        try:
            subprocess.run(blast_command, shell=True, check=True)
            with open(output_file, 'r') as f:
                result = f.read()
        except subprocess.CalledProcessError as e:
            result = f"Error running BLAST: {e}"

    return render_template('index.html', result=result)

if __name__ == '__main__':
    app.run(debug=True)

