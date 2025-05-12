# Machine Learning Inference Web Platform

This is a simple web application that allows users to input a text message and predict whether the message is "ham" or "spam" using a machine learning model. The app is built with Flask and uses a machine learning model trained on the SMS Spam Collection dataset.

## Technologies Used

- **Flask**: A lightweight WSGI web application framework in Python.
- **scikit-learn**: A machine learning library for Python, used for model training and prediction.
- **pandas**: A data manipulation and analysis library for Python.
- **Docker**: A platform to develop, ship, and run applications inside containers.

## Project Structure

- **app.py**: Main Flask application that handles the web server and routing.
- **train_model.py**: Script to train the machine learning model using the SMS Spam dataset.
- **model.pkl**: Pickled machine learning model used for predictions.
- **templates/**: Contains HTML templates for the web pages.
  - **index.html**: Main page where users can input a message.
- **static/**: Contains static files like CSS and images.
- **data/**: Contains the dataset used to train the model (SMS Spam Collection Dataset).

## Setup

### Clone the repository:
   ```bash
   git clone https://github.com/your-username/ml-inference-web-platform.git
   cd ml-inference-web-platform
   ```

### Create a conda environment:

```bash
conda create -n ml-inference-web-platform python=3.12
conda activate ml-inference-web-platform
```

### Install the required dependencies:

```bash
pip install -r requirements.txt
```
### Train the model (only required once):

```bash
python train_model.py
```

### Run the Flask app:

```bash
python app.py
```

### Open the web application:

Go to http://127.0.0.1:5000 in your web browser.

Enter a message to predict whether it's "ham" or "spam".

### Docker Setup (Optional)
You can also run the application inside a Docker container. Follow the steps below:

### Build the Docker image:

```bash
docker build -t ml-inference-web-platform .
```

### Run the Docker container:

```bash
docker run -p 5000:5000 ml-inference-web-platform
```

Access the web application at http://127.0.0.1:5000.

## Model Training
The model is trained using the SMS Spam Collection Dataset. The training script (train_model.py) does the following:

- Loads and preprocesses the dataset.

- Converts text data into numerical features using TF-IDF.

- Trains a machine learning model (e.g., Logistic Regression, Random Forest, etc.).

- Saves the trained model as a pickle file (model.pkl), which is used in the Flask application for predictions.