# Visualization Workbench

A modular Django-based bioinformatics web application that integrates tools for genomic data visualization, machine learning prediction, variant annotation, and pathway enrichment analysis.

---

## Project Structure

```
visualization_workbench/
├── manage.py                # Django management script
├── db.sqlite3               # SQLite DB for development
├── start_app.sh             # Shell script to run the app
├── data/                    # Input/output data directory
├── home/                    # Landing page app
├── igv_viewer/              # Genomic viewer integration (IGV)
├── pipeline_manager/        # Workflow management tools
├── ml_predictor/            # Machine learning prediction tools
├── variant_annotation/      # Variant annotation functionality
├── pathway_enrichment/      # Pathway enrichment & GSEA
└── visualization_workbench/ # Main Django project config
```

---

## Getting Started

### 1. Clone the Repository

```bash
git clone https://github.com/your-username/visualization_workbench.git
cd visualization_workbench
```

### 2. Set Up Virtual Environment

```bash
python3 -m venv env
source env/bin/activate  # or `source env/bin/activate.fish` if using fish shell
```

### 3. Install Dependencies

```bash
pip install -r requirements.txt
```

> 💡 If you don’t have `requirements.txt`, install common dependencies manually:

```bash
pip install django gseapy pandas
```

### 4. Run Migrations

```bash
python manage.py migrate
```

### 5. Start the Server

```bash
./start_app.sh
# or
python manage.py runserver
```

Visit [http://127.0.0.1:8000](http://127.0.0.1:8000)

---

## Modules Overview

| Module               | Description                                                  |
|----------------------|--------------------------------------------------------------|
| `home/`              | Landing page for the application                             |
| `igv_viewer/`        | Visualize genomic tracks via IGV.js                          |
| `pipeline_manager/`  | Launch and monitor pipelines or workflows                    |
| `ml_predictor/`      | Apply ML models to biological datasets                       |
| `variant_annotation/`| Annotate VCF or variant files using bioinformatics tools     |
| `pathway_enrichment/`| Run Enrichr or GSEA on gene lists and visualize results      |

---

## Pathway Enrichment Input Formats

- **Enrichr**: `.tsv` file with **one gene per line**
- **GSEA**: `.tsv` file with **two columns**: `gene` and `score`

---

## Utilities

- `start_app.sh`: Custom shell script to launch the Django app
- `data/`: Temporary directory for uploaded and processed files

