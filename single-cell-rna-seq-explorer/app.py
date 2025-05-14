from flask import Flask
from routes.home import home
from routes.visualization import visualize
from routes.download import download
from routes.differential_analysis import differential_analysis  # Import the new differential analysis module
from routes.gsea import gsea
from routes.go_enrichment import go_enrichment
from routes.gene_coexpression import gene_coexpression  # Import the new gene co-expression module

app = Flask(__name__)

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

if __name__ == '__main__':
    app.run(host='127.0.0.1', port=5000, debug=True)