import pandas as pd
import numpy as np

# Number of cells and genes
num_cells = 50
num_genes = 10

# Generate random data for gene expression
data = np.random.rand(num_cells, num_genes) * 10  # Random values between 0 and 10

# Generate sample gene names and cell names
genes = [f'Gene{i+1}' for i in range(num_genes)]
cells = [f'Cell{i+1}' for i in range(num_cells)]

# Create a DataFrame
df = pd.DataFrame(data, columns=genes, index=cells)

# Save as CSV
df.to_csv("synthetic_scRNA_seq_data.csv")

print("Synthetic data saved as 'synthetic_scRNA_seq_data.csv'.")
