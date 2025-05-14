import matplotlib
matplotlib.use('Agg')  # Use the 'Agg' backend for non-interactive plotting

from flask import request, render_template, jsonify
import pandas as pd
import networkx as nx  # Assuming you are using networkx for gene co-expression network
import matplotlib.pyplot as plt
import io
import base64

# This function will process the file and generate the co-expression network
def process_coexpression_data(file):
    # Read data from the uploaded file (assuming CSV format for simplicity)
    data = pd.read_csv(file)
    
    # Remove the gene names column before calculating correlation
    data = data.drop(columns=['Gene'])
    
    # Create a correlation matrix of gene expression values
    correlation_matrix = data.corr()

    # Create a graph from the correlation matrix (network of genes)
    G = nx.Graph()

    # Add nodes and edges with weights based on correlation
    for i in range(len(correlation_matrix.columns)):
        for j in range(i+1, len(correlation_matrix.columns)):
            weight = correlation_matrix.iloc[i, j]
            if weight > 0.8:  # Only include edges with high correlation (can be adjusted)
                G.add_edge(correlation_matrix.columns[i], correlation_matrix.columns[j], weight=weight)
    
    return G

# This function will generate the plot of the co-expression network
def generate_network_plot(G):
    # Create a plot of the co-expression network
    plt.figure(figsize=(10, 8))
    pos = nx.spring_layout(G)  # Layout of the network (can be changed)
    nx.draw_networkx_nodes(G, pos)
    nx.draw_networkx_edges(G, pos)
    nx.draw_networkx_labels(G, pos)
    
    # Save the plot to a bytes buffer and encode it for embedding in HTML
    img_buf = io.BytesIO()
    plt.savefig(img_buf, format='png')
    img_buf.seek(0)
    img_base64 = base64.b64encode(img_buf.read()).decode('utf-8')
    plt.close()
    
    return img_base64

# Define the route for gene co-expression analysis
def gene_coexpression():
    if request.method == 'POST':
        # Check if the file is included in the request
        if 'file' not in request.files:
            return jsonify({'error': 'No file part'}), 400
        file = request.files['file']
        
        if file.filename == '':
            return jsonify({'error': 'No selected file'}), 400
        
        if file:
            # Process the file and generate the gene co-expression network
            G = process_coexpression_data(file)
            
            # Generate the network plot (network image)
            network_image = generate_network_plot(G)
            
            # Render the result in an HTML page with the image
            return render_template('gene_coexpression_result.html', network_image=network_image)
    
    return render_template('gene_coexpression.html')
