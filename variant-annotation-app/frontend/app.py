from flask import Flask, request, render_template, redirect, flash
import os
import subprocess
from collections import Counter

app = Flask(__name__)
app.secret_key = 'your-secret-key'

UPLOAD_FOLDER = 'data'
ANNOVAR_PATH = '/Users/manishkumar/Desktop/llm-bio-webapps/variant-annotation-app/bin/annovar'
HUMANDB_PATH = os.path.join(ANNOVAR_PATH, 'humandb')

if not os.path.exists(UPLOAD_FOLDER):
    os.makedirs(UPLOAD_FOLDER)


def summarize_variants_by_chromosome(variants):
    chroms = [v.get('Chr', 'Unknown') for v in variants]
    counts = Counter(chroms)
    return dict(counts)


class BaseAnnotator:
    def __init__(self, vcf_file_path, output_prefix):
        self.vcf_file_path = vcf_file_path
        self.output_prefix = output_prefix

    def run(self):
        raise NotImplementedError

    def parse_output(self):
        raise NotImplementedError


class AnnovarAnnotator(BaseAnnotator):
    def run(self):
        cmd = [
            os.path.join(ANNOVAR_PATH, 'table_annovar.pl'),
            self.vcf_file_path,
            HUMANDB_PATH,
            '-buildver', 'hg19',
            '-out', self.output_prefix,
            '-remove',
            '-protocol', 'refGene,clinvar_20210501',
            '-operation', 'g,f',
            '-nastring', '.',
            '-vcfinput'
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise Exception(f"ANNOVAR annotation failed: {result.stderr}")

    def parse_output(self):
        variants = []
        filepath = f"{self.output_prefix}.hg19_multianno.txt"
        if not os.path.exists(filepath):
            return variants

        with open(filepath, 'r') as f:
            headers = f.readline().strip().split('\t')
            for line in f:
                fields = line.strip().split('\t')
                variant = dict(zip(headers, fields))
                variants.append(variant)
        return variants


class VepAnnotator(BaseAnnotator):
    def run(self):
        # TODO: Implement VEP run command
        pass

    def parse_output(self):
        # TODO: Implement VEP output parsing
        return []


class SnpEffAnnotator(BaseAnnotator):
    def run(self):
        # TODO: Implement SnpEff run command
        pass

    def parse_output(self):
        # TODO: Implement SnpEff output parsing
        return []


class BcftoolsAnnotator(BaseAnnotator):
    def run(self):
        # TODO: Implement BCFtools run command
        pass

    def parse_output(self):
        # TODO: Implement BCFtools output parsing
        return []


@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        if 'vcf_file' not in request.files:
            flash('No file part')
            return redirect('/')

        file = request.files['vcf_file']
        if file.filename == '':
            flash('No selected file')
            return redirect('/')

        annotator_choice = request.form.get('annotator', 'annovar').lower()
        filepath = os.path.join(UPLOAD_FOLDER, file.filename)
        file.save(filepath)

        output_prefix = os.path.join(UPLOAD_FOLDER, 'annotated_output')

        annotator_map = {
            'annovar': AnnovarAnnotator,
            'vep': VepAnnotator,
            'snpeff': SnpEffAnnotator,
            'bcftools': BcftoolsAnnotator
        }

        annotator_class = annotator_map.get(annotator_choice)
        if not annotator_class:
            flash(f"Unknown annotator selected: {annotator_choice}")
            return redirect('/')

        annotator = annotator_class(filepath, output_prefix)

        try:
            annotator.run()
            variants = annotator.parse_output()
        except Exception as e:
            flash(str(e))
            return redirect('/')

        variant_counts = summarize_variants_by_chromosome(variants)
        return render_template('result.html', variants=variants, variant_counts=variant_counts)

    return render_template('index.html')
    

if __name__ == '__main__':
    app.run(debug=True)
