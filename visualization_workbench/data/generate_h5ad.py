import scanpy as sc
import pandas as pd
import numpy as np

# Create a small dummy dataset
adata = sc.AnnData(
    X=np.random.poisson(1.0, (100, 50)),  # 100 cells, 50 genes
    obs=pd.DataFrame(index=[f"cell_{i}" for i in range(100)]),
    var=pd.DataFrame(index=[f"gene_{i}" for i in range(50)])
)

# Save to file
adata.write("test_data.h5ad")
