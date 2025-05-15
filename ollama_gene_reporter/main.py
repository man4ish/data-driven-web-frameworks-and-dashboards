from fastapi import FastAPI, Request, UploadFile, Form
from fastapi.responses import HTMLResponse
from fastapi.templating import Jinja2Templates
import pandas as pd
import requests
import io
import plotly.express as px
import plotly.graph_objects as go
import plotly
import json

app = FastAPI()
templates = Jinja2Templates(directory="templates")


def ask_ollama(prompt: str, model: str, temperature: float) -> str:
    print(model)
    response = requests.post("http://localhost:11434/api/generate", json={
        "model": model,
        "prompt": prompt,
        "temperature": temperature,
        "stream": False
    })
    
    response.raise_for_status()
    return response.json()["response"]


def generate_summary_from_csv(file: UploadFile, user_prompt: str, chart_type: str, model: str, temperature: float) -> (str, str):
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

    summary = ask_ollama(llm_prompt, model=model, temperature=temperature)

    # Create chart
    chart_html = ""
    if "Gene" not in df.columns:
        df.insert(0, "Gene", [f"Gene_{i}" for i in range(len(df))])  # fallback

    melted = df.melt(id_vars=["Gene"], var_name="Sample", value_name="Expression")

    if chart_type == "box":
        fig = px.box(melted, x="Sample", y="Expression", color="Sample", title="Expression Distribution per Sample")
    elif chart_type == "scatter":
        fig = px.scatter(melted, x="Sample", y="Expression", color="Sample", title="Expression Scatter Plot")
    elif chart_type == "violin":
        fig = px.violin(melted, x="Sample", y="Expression", color="Sample", box=True, points="all", title="Violin Plot")
    elif chart_type == "heatmap":
        pivot_df = df.set_index("Gene")
        fig = px.imshow(pivot_df.T, labels=dict(x="Gene", y="Sample", color="Expression"),
                        title="Gene Expression Heatmap")
    else:
        fig = px.box(melted, x="Sample", y="Expression", color="Sample", title="Default Box Plot")

    chart_html = plotly.io.to_html(fig, include_plotlyjs='cdn')

    return summary, chart_html


@app.get("/", response_class=HTMLResponse)
async def index(request: Request):
    return templates.TemplateResponse("index.html", {"request": request})


@app.post("/analyze", response_class=HTMLResponse)
async def analyze(
    request: Request,
    file: UploadFile,
    prompt: str = Form(...),
    chart_type: str = Form("box"),
    model: str = Form("llama3"),
    temperature: float = Form(0.7)
):
    summary, chart_html = generate_summary_from_csv(file, prompt, chart_type, model, temperature)
    return templates.TemplateResponse("result.html", {
        "request": request,
        "summary": summary,
        "chart_html": chart_html
    })
