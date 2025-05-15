#!/bin/bash

# Make sure Ollama server is running:
echo "Make sure you run 'ollama serve' in another terminal window before starting this app."

# Run FastAPI server
uvicorn main:app --reload --host 127.0.0.1 --port 8000
