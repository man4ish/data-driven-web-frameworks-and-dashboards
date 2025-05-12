# Local BLAST Web Server using Flask

This project is a simple Flask web application that allows users to input a DNA sequence (via text or file upload), run a BLAST search in the background using a local BLAST database, and view the results directly on the webpage.

---

## Features

- Paste or upload a FASTA sequence
- Runs `blastn` using NCBI BLAST+ tools
- Displays results on the same page
- Works completely offline with local databases
- Easy to set up and run on macOS or Linux

---

## Requirements

- Python 3.x
- Flask
- NCBI BLAST+ tools
- macOS (tested) or Linux

---

## Setup Instructions

### 1. Install Dependencies

```bash
pip install flask
brew install blast
```

### 2. Clone or Create the Project Directory

```bash
git clone https://github.com/man4ish/blast_app.git
cd blast_app
```

### 3. Set Up a Local BLAST Database

Create a simple FASTA file:

```fasta
>seq1
ATGCGTACGTAGCTAGCTAGCTAG
>seq2
TTATCGATCGATCGATCGATCGA
```

Save it as `db_sequences.fasta`.

Then, run:

```bash
makeblastdb -in db_sequences.fasta -dbtype nucl -out mydb
```

### 4. Run the Flask App

```bash
python app.py
```

Open your browser and go to: `http://127.0.0.1:5000`

---

## File Structure

```
blast_app/
├── app.py                  # Main Flask backend
├── db_sequences.fasta      # Input FASTA for BLAST DB
├── mydb.*                  # BLAST database files (created via makeblastdb)
├── templates/
│   └── index.html          # Web UI template
└── uploads/                # Temporary sequence input files (auto-created)
```

---

## Example Sequence

Paste this into the textbox or save as a `.fasta` file and upload:

```fasta
>mytest
ATGCGTACGTAGC
```

---

## Notes

- The BLAST database must be in the same folder or you need to adjust the path in `app.py`.
- This app uses `blastn` and expects nucleotide sequences.
- You can change to `blastp` for protein sequences by updating the database and command.

---

## References

- NCBI BLAST+ Command Line Tools: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
- Flask Documentation: https://flask.palletsprojects.com/


