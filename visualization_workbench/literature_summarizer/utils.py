from Bio import Entrez
import requests

Entrez.email = "your_email@example.com"  # Required

def fetch_pubmed_abstracts(query, max_results=5):
    handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
    record = Entrez.read(handle)
    ids = record["IdList"]
    
    handle = Entrez.efetch(db="pubmed", id=ids, rettype="abstract", retmode="text")
    abstracts = handle.read()
    return abstracts

def summarize_with_llm(text, model="your-model", temperature=0.7):
    response = requests.post(
        "http://localhost:11434/v1/generate",
        json={"prompt": f"Summarize this biomedical text:\n\n{text}", "model": model, "temperature": temperature},
    )
    response.raise_for_status()
    return response.json()["response"]
