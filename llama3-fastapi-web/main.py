from fastapi import FastAPI, Request, Form
from fastapi.responses import HTMLResponse
from fastapi.templating import Jinja2Templates
import httpx
import logging

app = FastAPI()
templates = Jinja2Templates(directory="templates")

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

async def ask_ollama(question: str) -> str:
    url = "http://localhost:11434/api/generate"
    payload = {
        "model": "llama3",
        "prompt": question,
        "stream": False
    }
    logger.info(f"Sending question to Ollama: {question}")
    async with httpx.AsyncClient() as client:
        response = await client.post(url, json=payload, timeout=10)
        response.raise_for_status()
        answer = response.json().get("response", "")
        logger.info(f"Received answer from Ollama: {answer}")
        return answer

@app.get("/", response_class=HTMLResponse)
async def home(request: Request):
    return templates.TemplateResponse("index.html", {"request": request})

@app.post("/ask", response_class=HTMLResponse)
async def ask_question(request: Request, question: str = Form(...)):
    question = question.strip()
    if not question:
        answer = "Please enter a valid question."
    else:
        try:
            answer = await ask_ollama(question)
        except httpx.RequestError as e:
            logger.error(f"Request error while contacting Ollama: {e}")
            answer = "Sorry, I couldn't reach the Ollama server. Please try again later."
        except Exception as e:
            logger.error(f"Unexpected error: {e}")
            answer = "An unexpected error occurred. Please try again."
    return templates.TemplateResponse("index.html", {
        "request": request,
        "question": question,
        "answer": answer
    })
