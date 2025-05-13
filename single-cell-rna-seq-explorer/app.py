from flask import Flask
from routes.home import home
from routes.visualization import visualize
from routes.download import download  # Import download route

app = Flask(__name__)

# Homepage for file upload & gene input
app.add_url_rule('/', 'home', home, methods=['GET', 'POST'])

# Visualization (GET for query param-based, POST for form submission)
app.add_url_rule('/visualize', 'visualize', visualize, methods=['GET', 'POST'])

# Download generated result
app.add_url_rule('/download', 'download', download, methods=['GET'])

if __name__ == '__main__':
    app.run(host='127.0.0.1', port=5000, debug=True)
