from flask import Flask, request, render_template, redirect, url_for, flash
import os
import subprocess

app = Flask(__name__)
app.secret_key = 'your-secret-key'  # needed for flashing messages

UPLOAD_FOLDER = 'data'
ANNOVAR_PATH = '/Users/manishkumar/Desktop/llm-bio-webapps/variant-annotation-app/bin/annovar'
HUMANDB_PATH = os.path.join(ANNOVAR_PATH, 'humandb')

if not os.path.exists(UPLOAD_FOLDER):
    os.makedirs(UPLOAD_FOLDER)

def run_annovar(vcf_file_path, output_prefix):
    cmd = [
        os.path.join(ANNOVAR_PATH, 'table_annovar.pl'),
        vcf_file_path,
        HUMANDB_PATH,
        '-buildver', 'hg19',
        '-out', output_prefix,
        '-remove',
        '-protocol', 'refGene,clinvar_20210501',
        '-operation', 'g,f',
        '-nastring', '.',
        '-vcfinput'
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise Exception(f"Annotation failed: {result.stderr}")

def parse_annovar_output(output_prefix):
    """
    Parse annovar output file (output_prefix.hg19_multianno.txt) into a list of dict variants.
    For now, just read lines and split by tabs. Adjust columns as needed.
    """
    variants = []
    filepath = f"{output_prefix}.hg19_multianno.txt"
    if not os.path.exists(filepath):
        return variants
    
    with open(filepath, 'r') as f:
        headers = f.readline().strip().split('\t')
        for line in f:
            fields = line.strip().split('\t')
            variant = dict(zip(headers, fields))
            variants.append(variant)
    return variants

@app.route('/', methods=['GET'])
def index():
    return render_template('index.html')

@app.route('/annotate', methods=['POST'])
def annotate():
    if 'vcf_file' not in request.files:
        flash('No file part')
        return redirect(url_for('index'))
    
    file = request.files['vcf_file']
    if file.filename == '':
        flash('No selected file')
        return redirect(url_for('index'))
    
    filepath = os.path.join(UPLOAD_FOLDER, file.filename)
    file.save(filepath)

    output_prefix = os.path.join(UPLOAD_FOLDER, 'annotated_output')

    try:
        run_annovar(filepath, output_prefix)
        variants = parse_annovar_output(output_prefix)
    except Exception as e:
        flash(str(e))
        return redirect(url_for('index'))
    
    return render_template('result.html', variants=variants)

if __name__ == '__main__':
    app.run(debug=True)
