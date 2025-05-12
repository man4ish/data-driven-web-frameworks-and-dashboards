#!/bin/bash
set -e

IMAGE_NAME="dash-app"

# Check if image exists
if [[ "$(docker images -q $IMAGE_NAME 2> /dev/null)" == "" ]]; then
  echo "🔧 Building Docker image..."
  docker build -t $IMAGE_NAME .
else
  echo "✅ Docker image '$IMAGE_NAME' already exists. Skipping build."
fi

# Remove old container if running
docker rm -f dash-app-container 2>/dev/null || true

# Run the container
docker run -d --name dash-app-container -p 8050:8050 $IMAGE_NAME
