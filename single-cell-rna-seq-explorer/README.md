# Single-Cell RNA-seq Explorer

This web application allows users to upload a CSV file containing single-cell RNA sequencing data and visualize various aspects of the data using plots and heatmaps. Users can also specify genes of interest to generate gene expression heatmaps.

## Features

- **Upload CSV File**: Allows users to upload single-cell RNA-seq data in CSV format.
- **Visualization**: Generates PCA, UMAP, KMeans clustering, and autoencoder UMAP plots.
- **Gene Expression Heatmap**: Optionally generate a heatmap for selected genes.
- **Download Results**: Users can download the generated CSV file containing processed results.

## Technologies Used

- **Flask**: Web framework for Python.
- **Pandas**: Data processing and manipulation.
- **Matplotlib & Seaborn**: Data visualization.
- **Bootstrap**: Front-end styling.

## Prerequisites

Before you begin, ensure you have the following installed:

- Python 3.7 or higher
- Flask
- Pandas
- Matplotlib
- Seaborn

## Installation

1. Clone the repository to your local machine:

    ```bash
    git clone https://github.com/your-username/single-cell-rna-seq-explorer.git
    cd single-cell-rna-seq-explorer
    ```

2. Create and activate a virtual environment (optional but recommended):

    ```bash
    python3 -m venv venv
    source venv/bin/activate  # On Windows, use venv\Scripts\activate
    ```

3. Install the required dependencies:

    ```bash
    pip install -r requirements.txt
    ```

## Usage

1. Start the Flask server:

    ```bash
    python app.py
    ```

2. Visit `http://127.0.0.1:5000` in your web browser.

3. On the home page, upload a CSV file containing the RNA-seq data and specify any genes (comma-separated) you wish to visualize.

4. Click **Visualize** to generate the visualizations.

5. After visualizations are generated, you can download the processed results by clicking the **Download Results** button.

## File Structure
```
/single-cell-rna-seq-explorer
├── app.py                # Flask application entry point
├── /routes               # Contains route handlers for the app
│   ├── home.py           # Home page route
│   ├── visualization.py  # Visualization route
│   └── download.py       # Download route
├── /templates            # HTML templates
│   ├── visualize.html    # HTML page for visualizations
│   ├── home.html         # HTML page for the home page
│   └── error.html        # HTML page for error handling
├── /static               # Static files (e.g., CSS, JS)
├── /results              # Folder where processed results will be saved
├── /utils                # Utility functions for processing and plotting
│   ├── preprocessing.py  # Data loading and preprocessing functions
│   └── plotting.py       # Plotting functions
└── requirements.txt      # Python dependencies

```

## 🐳 Dockerizing and Deploying to AWS

### 1️⃣ Create a `.dockerignore` File

A `.dockerignore` file prevents unnecessary files from being added to your Docker image (similar to `.gitignore`).

**Example:**

```
__pycache__/
.git/
.env
*.pyc
*.pyo
```

---

### 2️⃣ Build the Docker Image

From your project directory, build the Docker image:

```bash
docker build -t flask-app .
```

---

### 3️⃣ Run the Docker Container Locally

After the image is built, run the application in a container:

```bash
docker run -p 80:80 flask-app
```

Visit [http://localhost](http://localhost) in your browser to verify it's running.

---

### 4️⃣ Push the Docker Image to AWS

Once confirmed locally, push your image to AWS.

#### ✅ Push to Amazon ECR (Elastic Container Registry)

**1. Create an ECR Repository:**

Go to the [ECR Console](https://console.aws.amazon.com/ecr) and create a new repository.

**2. Authenticate Docker to ECR:**

```bash
aws ecr get-login-password --region <aws-region> | docker login --username AWS --password-stdin <aws-account-id>.dkr.ecr.<aws-region>.amazonaws.com
```

**3. Tag the Docker Image:**

```bash
docker tag flask-app:latest <aws-account-id>.dkr.ecr.<aws-region>.amazonaws.com/my-repository:latest
```

**4. Push the Docker Image:**

```bash
docker push <aws-account-id>.dkr.ecr.<aws-region>.amazonaws.com/my-repository:latest
```

---

### 5️⃣ Deploy to AWS ECS (Elastic Container Service)

**1. Create an ECS Cluster:**

Use either the **Fargate** or **EC2** launch type.

**2. Create a Task Definition:**

Define your container settings (CPU, memory, ports, etc.).

**3. Run the Task:**

Launch the task in your ECS cluster.

**4. (Optional) Set Up Load Balancer:**

Set up an **Application Load Balancer (ALB)** to distribute traffic and manage availability.

---

### 6️⃣ (Optional) Set Up Auto-Scaling

ECS supports auto-scaling to manage traffic spikes and optimize resource use.

---

### 🔁 Alternatives to ECS

- **Elastic Beanstalk:** Simple deployment with managed infrastructure, load balancing, and scaling.
- **AWS Lambda + API Gateway:** Ideal for lightweight or serverless applications.
