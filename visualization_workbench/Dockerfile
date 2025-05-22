# Use official Python slim image
FROM python:3.12-slim

# Set working directory inside container
WORKDIR /app

# Copy requirements.txt and install dependencies
COPY requirements.txt .

RUN pip install --no-cache-dir -r requirements.txt

# Copy all project files into the container
COPY . .

# Expose port 8000 for Django
EXPOSE 8000

# Command to run Django server
CMD ["python", "manage.py", "runserver", "0.0.0.0:8000"]
