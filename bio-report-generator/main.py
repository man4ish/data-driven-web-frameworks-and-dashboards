from fastapi import FastAPI, Request, UploadFile, Form
from fastapi.responses import HTMLResponse
from fastapi.templating import Jinja2Templates
import pandas as pd
import requests
import io
import plotly.express as px
import plotly
import json

app = FastAPI()
templates = Jinja2Templates(directory="templates")

def ask_ollama(prompt: str) -> str:
    response = requests.post("http://localhost:11434/api/generate", json={
        "model": "llama3",
        "prompt": prompt,
        "stream": False
    })
    response.raise_for_status()
    return response.json()["response"]

def generate_summary_from_csv(file: UploadFile, user_prompt: str) -> (str, str):
    # Read CSV content
    contents = file.file.read()
    df = pd.read_csv(io.BytesIO(contents))

    # Basic summary
    gene_means = df.mean(numeric_only=True).round(2).to_dict()
    stats_summary = "\n".join([f"{k}: {v}" for k, v in gene_means.items()])

    # LLM prompt
    llm_prompt = f"""This is gene expression data:\n{df.head().to_csv(index=False)}\n
    Here are some summary statistics:\n{stats_summary}\n
    {user_prompt}"""

    summary = ask_ollama(llm_prompt)

    # Create a chart
    melted = df.melt(id_vars=["Gene"], var_name="Sample", value_name="Expression")
    fig = px.box(melted, x="Sample", y="Expression", color="Sample", title="Expression Distribution per Sample")
    chart_html = plotly.io.to_html(fig, include_plotlyjs='cdn')

    return summary, chart_html

@app.get("/", response_class=HTMLResponse)
async def index(request: Request):
    return templates.TemplateResponse("index.html", {"request": request})

@app.post("/analyze", response_class=HTMLResponse)
async def analyze(request: Request, file: UploadFile, prompt: str = Form(...)):
    summary, chart_html = generate_summary_from_csv(file, prompt)
    return templates.TemplateResponse("result.html", {
        "request": request,
        "summary": summary,
        "chart_html": chart_html
    })
