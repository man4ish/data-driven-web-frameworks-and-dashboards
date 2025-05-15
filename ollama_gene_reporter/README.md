# Bioinformatics Report Generator with Ollama LLM Integration

This application allows you to upload gene expression CSV files and generate AI-driven summary reports with interactive charts. It integrates multiple Ollama large language models (LLMs) for natural language analysis and supports customizable chart types and model parameters.

---

## Features

- Upload CSV gene expression files for analysis.
- Generate summary reports with AI models: LLaMA 3, Mistral, Gemma, Deepseek-r1.
- Choose different chart types: box plot, scatter plot, violin plot, heatmap.
- Configure model temperature for response creativity.
- Interactive Plotly charts embedded in the report.
- Easy-to-use web interface with FastAPI and Jinja2 templates.

---

## Requirements

- Python 3.8+
- Ollama CLI installed ([https://ollama.com](https://ollama.com))
- Docker (optional, if using Docker for Ollama models)
- Required Python packages:
  ```bash
  pip install fastapi uvicorn pandas requests plotly jinja2
Setup and Running
1. Install Ollama CLI
Download and install Ollama CLI from https://ollama.com following their instructions.

2. Pull Ollama Models
Use the ollama pull command to download the desired models locally:

```bash
ollama pull llama3
ollama pull mistral
ollama pull gemma
ollama pull deepseek-r1
```

Verify the models are downloaded:

```bash
ollama ls
```

Example output:
```
NAME                  ID              SIZE      MODIFIED
llama3:latest         365c0bd3c000    4.7 GB    4 hours ago
mistral:latest        abcdef123456    3.5 GB    2 days ago
gemma:latest          7890abcd1234    5.1 GB    1 day ago
deepseek-r1:latest    0a8c26691023    4.7 GB    3 months ago
```
3. Start Ollama Server
Start the Ollama API server (make sure port 11434 is free):
```
ollama serve
```
4. Run the FastAPI Application
From the project directory, run:

```bash
uvicorn main:app --reload
```
The app will be accessible at: http://localhost:8000

## Usage
- Open the web page in your browser.

- Upload your gene expression CSV file.

- Enter a natural language prompt for analysis.

- Select the chart type.

- Select the Ollama LLM model.

- Adjust the temperature (0.0 for deterministic, up to 1.0 for creative).

- Click Generate Report.

- View AI-generated summaries and interactive charts.

## API
The app exposes two endpoints:

- GET / - The main upload page.

- POST /analyze - Handles file upload, parameters, calls Ollama API, and returns results.

## Troubleshooting
Port 11434 in use: If you get address already in use error, find and kill the process:

```bash
lsof -i :11434
kill <PID>
```

- Model not found: Make sure you pulled the model with ollama pull <model_name>.

- Ollama server not running: Run ollama serve before starting the FastAPI app.

## Example CSV Format
The CSV should have a Gene column and sample expression columns, e.g.:

```
Gene	Sample1	Sample2	Sample3
GeneA	5.2	3.4	6.1
GeneB	2.8	4.0	3.7
```