import pandas as pd
import pickle
import copy
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras.callbacks import EarlyStopping
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Conv1D, MaxPooling1D, Flatten, Dropout
from tensorflow.keras.layers import Input,BatchNormalization, ReLU, Add
from tensorflow.keras.models import Model
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import os
from os.path import basename
import numpy as np
from tensorflow.keras.preprocessing.sequence import pad_sequences
from tensorflow.keras.preprocessing.sequence import pad_sequences
import glob
from train_disorder_from_confidence import prepare_data_sliding_window

# Define PROJECT_HOME as two levels up from the current script
PROJECT_HOME = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW_DIR = f'{PROJECT_HOME}/data/raw'
INTERIM_DIR = f'{PROJECT_HOME}/data/interim'
PROCESSED_DIR = f'{PROJECT_HOME}/data/processed'
TABLE_DIR = f'{PROJECT_HOME}/tables/disorder'
MODEL_DIR = f'{PROJECT_HOME}/models/disorder'


def prepare_data_for_per_position_classification_untiled(df, sequence_id_col, label_col, fixed_length = 100):
    grouped = df.groupby(sequence_id_col)

    sequences = []
    labels = []

    for _, group in grouped:
        features = group.select_dtypes(include=[np.number]).drop([label_col], axis=1).values
        label = group[label_col].values  # Assuming this is already prepared per position

        # Pad or truncate each sequence to the fixed_length
        padded_features = pad_sequences([features], maxlen=fixed_length, dtype='float32', padding='post', truncating='post')[0]
        padded_label = pad_sequences([label], maxlen=fixed_length, value=-1, padding='post', truncating='post')[0]  # Use a specific value for padding

        sequences.append(padded_features)
        labels.append(padded_label)

    return np.array(sequences), np.array(labels)

def prepare_data_for_per_position_classification(df, sequence_id_col, label_col, fixed_length=100):
    grouped = df.groupby(sequence_id_col)

    sequences = []
    labels = []

    for _, group in grouped:
        total_length = len(group)  # Total length of the protein sequence
        features = group.select_dtypes(include=[np.number]).drop([label_col], axis=1).values
        label = group[label_col].values  # Assuming this is already prepared per position
        
        num_tiles = len(features) // fixed_length
        for i in range(num_tiles):
            start_index = i * fixed_length
            end_index = start_index + fixed_length

            tile_features = features[start_index:end_index]
            tile_label = label[start_index:end_index]

            # Calculate and append the position encoding
            position_encoding = end_index / total_length
            position_feature = np.full((fixed_length, 1), position_encoding, dtype='float32')
            tile_features = np.hstack((tile_features, position_feature))

            sequences.append(tile_features)
            labels.append(tile_label)

        # Handle the remainder if necessary
        remainder = len(features) % fixed_length
        if remainder > 0:
            start_index = num_tiles * fixed_length
            tile_features = features[start_index:]
            tile_label = label[start_index:]

            # Append position encoding
            position_encoding = len(features) / total_length
            position_feature = np.full((remainder, 1), position_encoding, dtype='float32')
            tile_features = np.hstack((tile_features, position_feature))

            # Pad the remainder to fixed_length
            tile_features = np.pad(tile_features, ((0, fixed_length - remainder), (0, 0)), mode='constant', constant_values=0)
            tile_label = np.pad(tile_label, (0, fixed_length - remainder), mode='constant', constant_values=-1)

            sequences.append(tile_features)
            labels.append(tile_label)

    return np.array(sequences), np.array(labels)

# def prepare_data_sliding_window(df: pd.DataFrame, window_size, random_shuffle=False):
#     # Ensure only numeric columns are considered for features and target (except 'Accession' which is categorical)
#     numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
#     df_numeric = df[numeric_cols]

#     # We still need 'Accession' for grouping but ensure it's not in the numeric processing
#     accession_groups = df.groupby('Accession')

#     X_transformed_list = []
#     y_transformed_list = []

#     for accession, group in accession_groups:
#         # Make sure to use only numeric columns for X and y
#         X = group[numeric_cols[:-1]].to_numpy()  # all numeric columns except the last one are features
#         y = group[numeric_cols[-1]].to_numpy()  # the last numeric column is the target
        
#         num_samples, num_features = X.shape
#         half_window = window_size // 2

#         if num_samples >= window_size:
#             X_group_transformed = np.empty((num_samples - 2 * half_window, num_features * window_size))
#             y_group_transformed = y[half_window:-half_window]

#             for i in range(half_window, num_samples - half_window):
#                 window_features = X[i - half_window:i + half_window + 1].flatten()
#                 X_group_transformed[i - half_window] = window_features

#             X_transformed_list.append(X_group_transformed)
#             y_transformed_list.append(y_group_transformed)

#     X_transformed = np.concatenate(X_transformed_list, axis=0)
#     y_transformed = np.concatenate(y_transformed_list, axis=0)

#     if random_shuffle:
#         # Shuffle the dataset to randomize rows
#         indices = np.random.permutation(X_transformed.shape[0])
#         X_transformed = X_transformed[indices]
#         y_transformed = y_transformed[indices]

#     return X_transformed, y_transformed


def read_csv_if_string(df, sep=None):
    if isinstance(df, str):
        if sep is None:
            if df.lower().endswith(".tsv"):
                sep = '\t'
            elif df.lower().endswith(" "):
                sep = ' '
            else:
                sep = ',' # regular CSV file
        print("Reading data table", df)
        df = pd.read_csv(df, sep=sep)
    return df


def run_disorder_prediction(df_test, model_file, scaler_file=None,
    window_size=17,
    model_name='dense',
    shuffle_columns=False, # can be a list of column names
    drop_columns=['Unnamed: 0'], # , 'ResNo'],
    ):
    '''
    We are not relying on `test_train_split` because we
    want to avoid the same accession codes ending up in
    both test and train data. Instead the split is managed
    by the script `data_disprotstructure.py`

    Args:
    df_test: dataframe or filename of test data
    model_file(str): filename of model to load
    model_name(str): 'dense', 'conv1d', 'resnet'
    '''
    df_test = read_csv_if_string(df_test)
    df_test = df_test.drop(columns=drop_columns, errors='ignore')
    print("Column names:", list(df_test.columns))
    if model_name in ['conv1d', 'resnet']:
        assert isinstance(df_test,pd.DataFrame)
        X_test, y_test=prepare_data_for_per_position_classification(df_test, "Accession", "structured")
    elif model_name == 'dense':
        # Remove text columns
        # df_train_numeric = df_train.select_dtypes(include=[np.number])
        # Separate features and target
        # X = df_numeric.iloc[:, :-1]  # all columns except the last one
        prep_result = prepare_data_sliding_window(df_test, window_size=window_size)
        X_test=prep_result['X']
        y_test = prep_result['y']
        feature_names = prep_result['columns']
        # y = df_numeric.iloc[:, -1]   # the last column
    else:
        raise ValueError("Error: supported model modes are resnet, dense, conv1d")
    if shuffle_columns is not None:
        if isinstance(shuffle_columns, bool) and shuffle_columns:
            for i in range(X_test.shape[1]):
                # randomly shuffle column i
                np.random.shuffle(X_test[:, i])
    # Split data into training and test sets
    # X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
    print("Reading model from file", model_file)
    model = tf.keras.models.load_model(model_file)
    if scaler_file is None:
        scaler_file = model_file.replace(".keras", '_scaler.pkl')
    if not os.path.exists(scaler_file):
        raise FileNotFoundError("Cannot find file for scaler:", scaler_file)
    with open(scaler_file, 'rb') as f:
        scaler = pickle.load(f)
    # Normalize the features
    if model_name == 'dense':
        X_test_scaled = scaler.transform(X_test)
    else:
        X_test_reshaped = X_test.reshape(-1, X_test.shape[-1])
        X_test_scaled = scaler.transform(X_test_reshaped)
        X_test_scaled = X_test_scaled.reshape(X_test.shape)
    print("Test feature data shape:", X_test_scaled.shape)
    print("Test label data shape:", y_test.shape)

    model.summary()
    # Train the model
    # model.fit(X_train_scaled, y_train, epochs=epochs, validation_split=0.1)

    # Evaluate the model on test data
    print("Model evaluation results:")
    model.evaluate(X_test_scaled, y_test)
    predictions = model.predict(X_test_scaled)
    if shuffle_columns:
        for i in range(X_test.shape[1]): # number of cases:
            X_test_scaled_cp = copy.deepcopy(X_test_scaled)
            np.random.shuffle(X_test_scaled_cp[:, i])
            print("Evaluating performance for shuffling of", i,feature_names[i])
            model.evaluate(X_test_scaled_cp, y_test)

def run_main_test_all(
    input_dir=os.path.join(PROCESSED_DIR,'disorder'),
    model_name='dense',
    model_filenames = None,
    model_dir = MODEL_DIR, 
    window_size=15,
    shuffle_columns=False,
    regions=['flexible linker/spacer',
    'molecular adaptor activity', 'molecular function regulator',
    'phosphorylation display site',
    'lipid binding',
    'ion binding',
    'pre-molten globule',
    'nucleic acid binding','disorder', 'protein binding', 'disorder to order']):
    print("Reading test data from directory", input_dir)
    print("Reading models from directory", model_dir)
    for region in regions:
        region_name = region.replace("/","_").replace(" ", "_")
        test_filenames = glob.glob(f"{input_dir}/disorder_{region_name}_*_test.tsv")
        print(f"Test file names for region {region}:", test_filenames)
        if model_filenames is None:
            model_filenames = glob.glob(f"{model_dir}/disorder_{region_name}_*.keras")
        for test_file in test_filenames:
            for model_file in model_filenames:
                print(f"Starting testing model {basename(model_file)} using data from {basename(test_file)} window-size: w{window_size}")
                if not os.path.exists(model_file):
                    print("Warning: could not find model file at", model_file)
                    continue
                try:
                    run_disorder_prediction(df_test=test_file,model_name=model_name, model_file=model_file,
                    window_size=window_size,
                    shuffle_columns=shuffle_columns)
                except Exception as e:
                    print(f"Encountered exception during evaluation: {e}")
                    raise e


if __name__ == "__main__":
    run_main_test_all(model_name='dense',
    model_filenames=['/Users/eckart/spyder/discomplex/models/disorder/disorder_disorder_5584_708_dense_w5.keras'],
    window_size=5,
    shuffle_columns=False)
