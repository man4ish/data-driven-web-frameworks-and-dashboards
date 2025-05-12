# Dash Bioinformatics Dashboard

This repository contains a Dash web application for visualizing bioinformatics data interactively. The app allows users to explore genomic data, create dynamic plots, and generate insights through a user-friendly web interface. The app is built with Python and Dockerized for easy deployment.

## Project Structure
```
dash-bio-dashboard/
├── dash_app.py # Main Dash app Python script
├── Dockerfile # Docker configuration file for app deployment
├── requirements.txt # List of dependencies required to run the app
└── README.md # Documentation for the repository
```


## Requirements

Before running the app, ensure you have the following dependencies installed:

- Python 3.12 or higher
- Docker (for containerization)

You can install the required Python packages using the `requirements.txt` file:

```bash
pip install -r requirements.txt
```

## Running the App Locally
### Clone the repository:
```bash
git clone https://github.com/your-username/dash-bio-dashboard.git
cd dash-bio-dashboard
```

## Install the required dependencies:

```bash
pip install -r requirements.txt
```
Run the Dash app:
```bash
python dash_app.py
```
The app will be accessible at http://localhost:8050 in your web browser.

### Dockerizing the App
This app is also Dockerized for easy deployment. To build and run the app in a Docker container:

Build the Docker image:

```bash
docker build -t dash-app .
```

### Run the Docker container:

```bash
docker run -p 8050:8050 dash-app
```

The app will be accessible at http://localhost:8050 in your web browser.

## Usage
- Interact with various dynamic plots in the dashboard.

- Explore bioinformatics data visually.

- Use filters and selection options to modify the content displayed on the plots.

