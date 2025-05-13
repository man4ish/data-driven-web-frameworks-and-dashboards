import pandas as pd
import scanpy as sc

def load_and_preprocess(filepath):
    df = pd.read_csv(filepath, index_col=0)
    df = df.apply(pd.to_numeric, errors='coerce').dropna()

    if df.empty:
        return {
            'error': f"""
                <h3>No valid numeric data.</h3>
                <pre>{pd.read_csv(filepath).head().to_string()}</pre>
            """
        }

    adata = sc.AnnData(df)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    return {'df': df, 'adata': adata}
