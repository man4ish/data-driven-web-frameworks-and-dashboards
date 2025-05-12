import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report
import pickle

# Read CSV with correct encoding and handle extra commas
df = pd.read_csv("data/spam.csv", encoding="latin1", usecols=[0, 1], names=["label", "text"], header=0)

# Optional: Check the first few rows
print(df.head())

# Preprocess the data
X = df['text']
y = df['label']

# Split the data
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Vectorizing the text data
vectorizer = TfidfVectorizer(stop_words='english')
X_train_vectorized = vectorizer.fit_transform(X_train)
X_test_vectorized = vectorizer.transform(X_test)

# Train a Random Forest Classifier
model = RandomForestClassifier()
model.fit(X_train_vectorized, y_train)

# Evaluate the model
y_pred = model.predict(X_test_vectorized)
print(classification_report(y_test, y_pred))

# Save the model and vectorizer to disk
with open('model.pkl', 'wb') as model_file:
    pickle.dump(model, model_file)
    
with open('vectorizer.pkl', 'wb') as vectorizer_file:
    pickle.dump(vectorizer, vectorizer_file)

