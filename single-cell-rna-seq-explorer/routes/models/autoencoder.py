# routes/models/autoencoder.py
from keras.models import Sequential
from keras.layers import Dense

def build_autoencoder(X):
    """
    Build a simple autoencoder for dimensionality reduction.
    """
    model = Sequential()
    model.add(Dense(128, activation='relu', input_dim=X.shape[1]))
    model.add(Dense(64, activation='relu'))
    model.add(Dense(32, activation='relu'))
    model.add(Dense(2, activation='linear'))  # 2 for 2D visualization
    model.add(Dense(32, activation='relu'))
    model.add(Dense(64, activation='relu'))
    model.add(Dense(128, activation='relu'))
    model.add(Dense(X.shape[1], activation='sigmoid'))

    model.compile(optimizer='adam', loss='mse')
    model.fit(X, X, epochs=50, batch_size=32, verbose=0)
    return model
