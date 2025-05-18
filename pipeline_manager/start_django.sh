#!/bin/bash

PORT=8000

echo "Checking if port $PORT is in use..."

PID=$(lsof -ti tcp:$PORT)

if [ -n "$PID" ]; then
    echo "Port $PORT is in use by process ID(s): $PID"
    echo "Killing process(es)..."
    kill -9 $PID
    echo "Port $PORT is now free."
else
    echo "Port $PORT is already free."
fi

echo "Checking for model changes..."
CHANGES=$(python manage.py makemigrations --check --dry-run 2>&1)

if echo "$CHANGES" | grep -q "No changes detected"; then
    echo "No model changes detected."
else
    echo "Making migrations..."
    python manage.py makemigrations

    echo "Applying migrations..."
    python manage.py migrate
fi

echo "Starting Django development server on port $PORT..."
python manage.py runserver $PORT

