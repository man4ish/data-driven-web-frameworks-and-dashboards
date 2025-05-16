import os 
import subprocess
from typing import List
import pandas as pd

def run_annovar(vcf_path: str, output_dir: str, humandb_path: str) -> str:
    base = os.path.basename(vcf_path).replace(".vcf", "")
    out_prefix = os.path.join(output_dir, base)
    import subprocess

    cmd = [
        "/Users/manishkumar/Desktop/llm-bio-webapps/variant-annotation-app/bin/annovar/table_annovar.pl",
        "data/input.vcf",
        "/Users/manishkumar/Desktop/llm-bio-webapps/variant-annotation-app/bin/annovar/humandb",
        "-buildver", "hg19",
        "-out", "data/input",
        "-remove",
        "-protocol", "refGene,clinvar_20210501",
        "-operation", "g,f",
        "-nastring", ".",
        "-vcfinput"
    ]

    try:
        subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        print("Command failed.")
        print("STDOUT:\n", e.stdout.decode())
        print("STDERR:\n", e.stderr.decode())
    return out_prefix + ".hg19_multianno.txt"