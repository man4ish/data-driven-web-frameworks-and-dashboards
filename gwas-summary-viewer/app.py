from flask import Flask, request, render_template, redirect, url_for
import os
import pandas as pd
import numpy as np
import plotly.express as px
import uuid

app = Flask(__name__)
UPLOAD_FOLDER = 'uploads'
os.makedirs(UPLOAD_FOLDER, exist_ok=True)

@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        file = request.files['file']
        if file:
            filepath = os.path.join(UPLOAD_FOLDER, f"{uuid.uuid4()}_{file.filename}")
            file.save(filepath)
            return redirect(url_for('visualize', filepath=filepath))
    return render_template('index.html')

@app.route('/visualize')
def visualize():
    filepath = request.args.get('filepath')
    df = pd.read_csv(filepath, sep=None, engine='python')  # handles CSV/TSV

    # Basic column name standardization
    df.columns = [col.upper() for col in df.columns]
    
    # Required columns: SNP, CHR, BP, P
    if not {'CHR', 'BP', 'P'}.issubset(df.columns):
        return "Uploaded file must contain CHR, BP, and P columns."

    # Manhattan Plot
    # Apply the transformation using numpy's log10
    df['-LOG10P'] = df['P'].apply(lambda p: pd.NA if p <= 0 else -np.log10(p))
    manhattan_fig = px.scatter(df, x='BP', y='-LOG10P', color='CHR', title='Manhattan Plot')
    manhattan_plot = manhattan_fig.to_html(full_html=False)

    # QQ Plot
    df_sorted = df.sort_values('P')
    df_sorted = df_sorted[df_sorted['P'] > 0]
    df_sorted['EXP_P'] = [i / len(df_sorted) for i in range(1, len(df_sorted)+1)]
    df_sorted['EXP_LOG10P'] = -df_sorted['EXP_P'].apply(lambda p: np.log10(p))
    df_sorted['OBS_LOG10P'] = -df_sorted['P'].apply(lambda p: np.log10(p))
    qq_fig = px.scatter(df_sorted, x='EXP_LOG10P', y='OBS_LOG10P', title='QQ Plot')
    qq_plot = qq_fig.to_html(full_html=False)

    return render_template('visualize.html', manhattan_plot=manhattan_plot, qq_plot=qq_plot)

if __name__ == '__main__':
    app.run(debug=True)