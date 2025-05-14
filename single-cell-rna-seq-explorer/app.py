from flask import Flask
from routes.home import home
from routes.visualization import visualize
from routes.download import download
from routes.differential_analysis import differential_analysis  # Import the new differential analysis module
from routes.gsea import gsea
from routes.go_enrichment import go_enrichment
from routes.gene_coexpression import gene_coexpression  # Import the new gene co-expression module
from routes.pathway_enrichment import pathway_enrichment  # Import the new pathway enrichment module
import os

app = Flask(__name__)

# Set the folder for uploads
app.config['UPLOAD_FOLDER'] = 'uploads'  # Path to the folder where files will be uploaded
app.config['STATIC_FOLDER'] = 'static'   # Path to the folder where static files (like plots) are stored

# Set a maximum file size (optional, but useful for large files)
app.config['MAX_CONTENT_LENGTH'] = 16 * 1024 * 1024  # Max file size is 16MB

# Ensure the upload and static directories exist
if not os.path.exists(app.config['UPLOAD_FOLDER']):
    os.makedirs(app.config['UPLOAD_FOLDER'])

if not os.path.exists(os.path.join(app.config['STATIC_FOLDER'], 'plots')):
    os.makedirs(os.path.join(app.config['STATIC_FOLDER'], 'plots'))

# Homepage for file upload & gene input
app.add_url_rule('/', 'home', home, methods=['GET', 'POST'])

# Visualization (GET for query param-based, POST for form submission)
app.add_url_rule('/visualize', 'visualize', visualize, methods=['GET', 'POST'])

# Download generated result
app.add_url_rule('/download', 'download', download, methods=['GET'])

# Add GSEA route
app.add_url_rule('/gsea', 'gsea', gsea, methods=['POST', 'GET']) 

# Differential analysis route
app.add_url_rule('/differential', 'differential', differential_analysis, methods=['GET', 'POST'])

# Gene Ontology (GO) enrichment analysis route
app.add_url_rule('/go_enrichment', 'go_enrichment', go_enrichment, methods=['POST', 'GET'])

# Add Gene Co-expression Network route
app.add_url_rule('/gene_coexpression', 'gene_coexpression', gene_coexpression, methods=['POST', 'GET'])  # New route

# Add Pathway Enrichment route
app.add_url_rule('/pathway_enrichment', 'pathway_enrichment', pathway_enrichment, methods=['POST', 'GET'])  # New route


if __name__ == '__main__':
    app.run(host='127.0.0.1', port=5000, debug=True)