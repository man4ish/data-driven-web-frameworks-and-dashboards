# Single-Cell RNA-seq Explorer

This web application allows users to upload a CSV file containing single-cell RNA sequencing data and visualize various aspects of the data using plots and heatmaps. Users can also specify genes of interest to generate gene expression heatmaps.

## Features

- **Upload CSV File**: Allows users to upload single-cell RNA-seq data in CSV format.
- **Visualization**: Generates PCA, UMAP, KMeans clustering, and autoencoder UMAP plots.
- **Gene Expression Heatmap**: Optionally generate a heatmap for selected genes.
- **Download Results**: Users can download the generated CSV file containing processed results.

## Technologies Used

- **Flask**: Web framework for Python.
- **Pandas**: Data processing and manipulation.
- **Matplotlib & Seaborn**: Data visualization.
- **Bootstrap**: Front-end styling.

## Prerequisites

Before you begin, ensure you have the following installed:

- Python 3.7 or higher
- Flask
- Pandas
- Matplotlib
- Seaborn

## Installation

1. Clone the repository to your local machine:

    ```bash
    git clone https://github.com/your-username/single-cell-rna-seq-explorer.git
    cd single-cell-rna-seq-explorer
    ```

2. Create and activate a virtual environment (optional but recommended):

    ```bash
    python3 -m venv venv
    source venv/bin/activate  # On Windows, use venv\Scripts\activate
    ```

3. Install the required dependencies:

    ```bash
    pip install -r requirements.txt
    ```

## Usage

1. Start the Flask server:

    ```bash
    python app.py
    ```

2. Visit `http://127.0.0.1:5000` in your web browser.

3. On the home page, upload a CSV file containing the RNA-seq data and specify any genes (comma-separated) you wish to visualize.

4. Click **Visualize** to generate the visualizations.

5. After visualizations are generated, you can download the processed results by clicking the **Download Results** button.

## File Structure
```
/single-cell-rna-seq-explorer
├── app.py                # Flask application entry point
├── /routes               # Contains route handlers for the app
│   ├── home.py           # Home page route
│   ├── visualization.py  # Visualization route
│   └── download.py       # Download route
├── /templates            # HTML templates
│   ├── visualize.html    # HTML page for visualizations
│   ├── home.html         # HTML page for the home page
│   └── error.html        # HTML page for error handling
├── /static               # Static files (e.g., CSS, JS)
├── /results              # Folder where processed results will be saved
├── /utils                # Utility functions for processing and plotting
│   ├── preprocessing.py  # Data loading and preprocessing functions
│   └── plotting.py       # Plotting functions
└── requirements.txt      # Python dependencies

```