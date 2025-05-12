# GWAS Summary Viewer

This web application visualizes GWAS (Genome-Wide Association Study) summary statistics, specifically using Manhattan and QQ plots. Users can upload a `.csv` or `.tsv` file containing their GWAS data, and the app will generate interactive plots for further analysis.

## Features

- **Upload GWAS Results**: Users can upload `.csv` or `.tsv` files containing GWAS data.
- **Manhattan Plot**: A plot visualizing the association between SNPs (Single Nucleotide Polymorphisms) and their respective p-values.
- **QQ Plot**: A quantile-quantile plot to compare observed versus expected p-values.
- **Dynamic Visualizations**: Interactive plots built using Plotly for detailed data exploration.

## Required Columns in the Data File

The uploaded file must contain the following columns:

- **CHR**: Chromosome number.
- **BP**: Base pair position.
- **P**: P-value for the SNP.

## Getting Started

### Prerequisites

Make sure you have Python 3.7 or higher installed. The following packages are required:

- Flask
- Pandas
- Plotly

You can install them using pip or conda:
```
pip install -r requirements.txt
```


### Installation

1. Clone this repository:
   ```bash
   git clone https://github.com/yourusername/gwas-summary-viewer.git
Navigate to the project directory:

```bash
cd gwas-summary-viewer
```

Install the required dependencies:

```bash
pip install -r requirements.txt
```

### Running the App
To run the application locally, use the following command:
```
python app.py
```
The app will be accessible at http://127.0.0.1:5000/.

### Uploading Data
- Visit the main page of the app.

- Upload your .csv or .tsv file with the required columns (CHR, BP, P).

The app will process the file and generate the visualizations.

### File Format
- Supported file types: .csv or .tsv.

- Required columns: The file must contain the columns CHR, BP, and P.

### Example Data
Here’s an example of a valid GWAS data file:

```
CHR,BP,P
1,1000,0.001
1,2000,0.02
2,1500,0.03
3,10000,0.005
```