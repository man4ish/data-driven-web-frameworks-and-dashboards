import os
from Bio import SeqIO

def parse_file_and_summarize(file_path: str, prompt: str) -> str:
    if file_path.endswith(".fasta") or file_path.endswith(".fa"):
        records = list(SeqIO.parse(file_path, "fasta"))
        count = len(records)
        lengths = [len(r.seq) for r in records]
        return f"FASTA file with {count} sequences. Avg length: {sum(lengths)//count}. Prompt: {prompt}"

    elif file_path.endswith(".vcf"):
        with open(file_path) as f:
            lines = [l for l in f if not l.startswith("#")]
        return f"VCF file with {len(lines)} variants. Prompt: {prompt}"

    elif file_path.endswith(".csv") or file_path.endswith(".tsv"):
        import pandas as pd
        df = pd.read_csv(file_path, sep=None, engine='python')
        return f"Expression file with {df.shape[0]} genes and {df.shape[1]} samples. Prompt: {prompt}"

    else:
        return "Unsupported file format"

