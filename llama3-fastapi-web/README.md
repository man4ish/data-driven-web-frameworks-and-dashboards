# LLaMA3 FastAPI Web Interface

A simple and elegant web interface built with **FastAPI** and **Jinja2** to interact with the [Ollama](https://ollama.com) LLaMA3 model running locally.

## Features

- Ask questions via a clean web form
- Sends requests asynchronously to the Ollama API
- Handles errors gracefully with user-friendly messages
- Logs requests and responses for easy debugging
- Simple, minimalistic UI powered by Jinja2 templates

## Requirements

- Python 3.8+
- [FastAPI](https://fastapi.tiangolo.com/)
- [httpx](https://www.python-httpx.org/) (for async HTTP requests)
- [Uvicorn](https://www.uvicorn.org/) (ASGI server)
- Ollama app running locally with the `llama3` model accessible at `http://localhost:11434/api/generate`

## Installation

### Clone this repository:

   ```bash
   git clone https://github.com/yourusername/llama3-fastapi-app.git
   cd llama3-fastapi-app
   ```

Create and activate a virtual environment (optional but recommended):

```bash
python3 -m venv venv
source venv/bin/activate
```

### Install dependencies:

```bash
pip install fastapi uvicorn httpx jinja2
```

Make sure the Ollama app is installed and running locally with the llama3 model:

```bash
ollama run llama3
```

### Usage
Run the FastAPI application using Uvicorn:

```bash
uvicorn main:app --reload
```

Open your browser and navigate to:

```
http://localhost:8000/
```

Type your question and submit the form to get an answer from the LLaMA3 model.

