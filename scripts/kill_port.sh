#!/bin/bash

# Check if port number is provided
if [ -z "$1" ]; then
  echo "Usage: $0 <port_number>"
  exit 1
fi

PORT=$1

# Find and kill process using the specified port
PIDS=$(lsof -i :$PORT | awk 'NR != 1 {print $2}' | sort -u)

if [ -z "$PIDS" ]; then
  echo "No process found using port $PORT."
  exit 0
fi

echo "Killing processes on port $PORT: $PIDS"
echo "$PIDS" | xargs kill -9

echo "Done."

