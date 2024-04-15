import pandas as pd
from sklearn.dummy import DummyClassifier
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score
import tensorflow as tf
import os
import numpy as np
import datetime

SLUG='disorder'
# Define PROJECT_HOME as two levels up from the current script
PROJECT_HOME = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW_DIR = f'{PROJECT_HOME}/data/raw'
INTERIM_DIR = f'{PROJECT_HOME}/data/interim'
PROCESSED_DIR = f'{PROJECT_HOME}/data/processed'
TABLE_DIR = f'{PROJECT_HOME}/tables/{SLUG}'
MODEL_DIR = f'{PROJECT_HOME}/models/{SLUG}'

data_file = os.path.join(TABLE_DIR, 'disorder_711.tsv')
df = pd.read_csv(data_file, sep='\t')

# Filter only numeric columns (excluding 'structured' for features)
X = df.select_dtypes(include=[np.number]).drop(columns='structured')
y = df['structured']
X.replace([np.inf, -np.inf, np.nan], 0.0, inplace=True)
y.replace([np.inf, -np.inf, np.nan], 0.0, inplace=True)
# Normalize the features
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# Split the data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X_scaled, y, test_size=0.2, random_state=42)

assert(not np.any(np.isnan(X_train)))
assert(not np.any(np.isinf(X_train)))
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Conv1D, MaxPooling1D, Flatten, Dense, Dropout

# Building the model
model = Sequential([
    Conv1D(32, 3, activation='relu', input_shape=(X_train.shape[1], 1)),
    MaxPooling1D(2),
    Conv1D(64, 3, activation='relu'),
    MaxPooling1D(2),
    Conv1D(64, 3, activation='relu'),
    MaxPooling1D(2),
    Flatten(),
    Dense(64, activation='relu'),
    Dropout(0.5),
    Dense(1, activation='sigmoid')
])
optimizer = tf.keras.optimizers.Adam(learning_rate=0.0001)
model.compile(optimizer=optimizer, loss='binary_crossentropy', metrics=['accuracy'])
model.summary()


# Reshape data for 1D ConvNet
X_train_reshaped = X_train.reshape((X_train.shape[0], X_train.shape[1], 1))
X_test_reshaped = X_test.reshape((X_test.shape[0], X_test.shape[1], 1))

dummy_clf = DummyClassifier(strategy='most_frequent') 
# Fit the Dummy Classifier on the training data
dummy_clf.fit(X_train, y_train)

# Make predictions on the test set
y_pred = dummy_clf.predict(X_test)

# Calculate the accuracy of the Dummy Classifier
dummy_accuracy = accuracy_score(y_test, y_pred)
print(f"Dummy Classifier Accuracy: {dummy_accuracy:.2f}")


# Train the model
history = model.fit(X_train_reshaped, y_train, epochs=10, validation_data=(X_test_reshaped, y_test))

now = datetime.datetime.now()

# Format date and time
# Format as YYYY-MM-DD
date_str = now.strftime("%Y-%m-%d")
# Calculate minutes since midnight
minutes_since_midnight = now.hour * 60 + now.minute

# Create a formatted string that includes date and minutes since midnight
formatted_time = "{}_{}min".format(date_str, minutes_since_midnight)
if not os.path.exists(MODEL_DIR):
    os.makedirs(MODEL_DIR)
outfile = f'{MODEL_DIR}/disorder_convnet_{formatted_time}.h5'
print("saving model to file", outfile)
model.save(outfile)
# Evaluate the model
model.evaluate(X_test_reshaped, y_test)

# Make predictions (optional)
predictions = model.predict(X_test_reshaped)
predicted_classes = (predictions > 0.5).astype(int)

