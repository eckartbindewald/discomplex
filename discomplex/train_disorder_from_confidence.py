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
import sys
import seaborn as sns
import matplotlib.pyplot as plt
import xgboost as xgb
from sklearn.metrics import accuracy_score, roc_auc_score

# Define PROJECT_HOME as two levels up from the current script
PROJECT_HOME = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW_DIR = f'{PROJECT_HOME}/data/raw'
INTERIM_DIR = f'{PROJECT_HOME}/data/interim'
PROCESSED_DIR = f'{PROJECT_HOME}/data/processed'
TABLE_DIR = f'{PROJECT_HOME}/tables/disorder'
MODEL_DIR = f'{PROJECT_HOME}/models/disorder'

# def prepare_data_for_per_position_classification(df, sequence_id_col, label_col):
#     grouped = df.groupby(sequence_id_col)
#     max_length = grouped.size().max()  # Find the maximum sequence length for padding

#     sequences = []
#     labels = []

#     for _, group in grouped:
#         features = group.drop([label_col, sequence_id_col,"AACode"], axis=1).values
#         label = group[label_col].values  # Assuming this is already prepared per position

#         # Padding each sequence to the max_length
#         padded_features = pad_sequences([features], maxlen=max_length, dtype='float32', padding='post')[0]
#         padded_label = pad_sequences([label], maxlen=max_length, value=-1, padding='post')[0]  # Use a specific value for padding

#         sequences.append(padded_features)
#         labels.append(padded_label)

#     return np.array(sequences), np.array(labels)


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

def prepare_data_sliding_window(df: pd.DataFrame, window_size, random_shuffle=False):
    # Ensure only numeric columns are considered for features and target (except 'Accession' which is categorical)
    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
    X_colnames = numeric_cols[:-1] # everything except last column (target)
    # df_numeric = df[numeric_cols]
    if window_size == 1: # no groups, each row is by itself:
        X = df[numeric_cols[:-1]].to_numpy()  # all numeric columns except the last one are features
        y = df[numeric_cols[-1]].to_numpy() 
        assert not random_shuffle, 'random_shuffle not supported for window-size 1'
        return {'X':X,'y':y, 'columns':X_colnames}

    # We still need 'Accession' for grouping but ensure it's not in the numeric processing
    accession_groups = df.groupby('Accession')

    X_transformed_list = []
    y_transformed_list = []
    half_window = window_size // 2

    for accession, group in accession_groups:
        # print("group for accession", accession, ":")
        # print(group.to_string())
        # assert False
        # Make sure to use only numeric columns for X and y
        X = group[numeric_cols[:-1]].to_numpy()  # all numeric columns except the last one are features
        y = group[numeric_cols[-1]].to_numpy()  # the last numeric column is the target
        
        num_samples, num_features = X.shape
        
        if num_samples >= window_size:
            X_group_transformed = np.empty((num_samples - 2 * half_window, num_features * window_size))
            y_group_transformed = y[half_window:-half_window]

            for i in range(half_window, num_samples - half_window):
                window_features = X[i - half_window:i + half_window + 1].flatten()
                X_group_transformed[i - half_window] = window_features

            X_transformed_list.append(X_group_transformed)
            y_transformed_list.append(y_group_transformed)

    X_transformed = np.concatenate(X_transformed_list, axis=0)
    y_transformed = np.concatenate(y_transformed_list, axis=0)

    if random_shuffle:
        # Shuffle the dataset to randomize rows
        indices = np.random.permutation(X_transformed.shape[0])
        X_transformed = X_transformed[indices]
        y_transformed = y_transformed[indices]
    colnames_all = []
    for i in range(-half_window, half_window+1):
        istr = str(i)
        if i > 0:
            istr = '+' + istr
        elif i == 0:
            istr = ''
        cols = copy.deepcopy(X_colnames)
        for j in range(len(cols)):
            cols[j] = cols[j] + istr
        colnames_all += cols
    return {'X':X_transformed, 'y':y_transformed, 'columns':colnames_all}


def _test_prepare_data_sliding_window(df, window_sizes=[3,1]):
    print("input dataframe:")
    print(df)
    for window_size in window_sizes:
        trafo_result = prepare_data_sliding_window(df, window_size=window_size)
        X_transformed = trafo_result['X']
        y_transformed = trafo_result['y']
        colnames = trafo_result['columns']
        num_samples = len(df) - 2 * (window_size // 2) * len(df['Accession'].unique())
        num_features = (df.shape[1] - 2) * window_size  # Minus 'Accession' and 'Target', times the window size # NO
        
        print(f"Testing window_size={window_size}:")
        print("X_transformed shape:", X_transformed.shape)
        print("y_transformed shape:", y_transformed.shape)
        print("X transformed:")
        print(colnames)
        print(X_transformed)
        print(y_transformed)
        # Check if the number of samples is correct
        expected_samples = df.groupby('Accession').apply(lambda g: max(0, len(g) - (window_size - 1))).sum()
        assert X_transformed.shape[0] == expected_samples, "Number of samples does not match expected value."
        
        # Check if the number of features is correct
        assert X_transformed.shape[1] == num_features, "Number of features does not match expected value."

        # Check the first and last element in y_transformed if possible
        if expected_samples > 0:
            assert y_transformed[0] == df[df['Accession'] == 'A']['Target'].iloc[window_size // 2], "First element of y_transformed is incorrect."
            if expected_samples > 1:
                assert y_transformed[-1] == df[df['Accession'] == 'B']['Target'].iloc[-1 - (window_size // 2)], "Last element of y_transformed is incorrect."

def test_prepare_data_sliding_window(window_sizes=[3,1]):
    # Run the test with a variety of window sizes
    # Sample DataFrame creation
    data = {
        'Accession': ['A', 'A', 'A', 'B', 'B', 'B', 'B', 'A', 'A', 'A'],
        'Feature1': [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7,-0.3,-0.2,-0.1],
        'Feature2': [1, 2, 3, 4, 5, 6, 7,-1,-2,-3],
        'Target': [0, 1, 0, 1, 1, 0, 1, 0, 1, 0]
    }
    df = pd.DataFrame(data)
    # test_window_sizes = [1, 2, 3, 4]  # Testing window size of 1 to 4
    _test_prepare_data_sliding_window(df, window_sizes)


def gen_conv1d_model_per_position(X):
    # Extract the sequence length and number of features directly from X
    sequence_length = X.shape[1]  # Assuming X is shaped as [samples, sequence_length, num_features]
    num_features = X.shape[2]
    
    # Create the model
    model = Sequential([
        Conv1D(filters=64, kernel_size=7, padding='same', activation='relu', input_shape=(sequence_length, num_features)),
        Conv1D(filters=32, kernel_size=5, padding='same', activation='relu'),
        # Dropout(0.5),
        Conv1D(filters=16, kernel_size=3, padding='same', activation='relu'),
        Conv1D(filters=1, kernel_size=1, padding='same', activation='sigmoid')  # Output layer: 1 filter per position
    ])
    
    return model


def residual_block(x, filters, kernel_size, downsample=False):
    stride = 2 if downsample else 1
    y = Conv1D(filters, kernel_size=kernel_size, padding='same', strides=stride)(x)
    y = BatchNormalization()(y)
    y = ReLU()(y)
    y = Conv1D(filters, kernel_size=kernel_size, padding='same')(y)
    y = BatchNormalization()(y)

    if downsample:
        x = Conv1D(filters, kernel_size=1, strides=2, padding='same')(x)

    out = Add()([x, y])
    out = ReLU()(out)
    return out


def build_resnet(input_shape, num_classes, num_blocks):
    inputs = Input(shape=input_shape)
    x = Conv1D(64, kernel_size=7, padding='same', strides=1)(inputs)
    x = BatchNormalization()(x)
    x = ReLU()(x)

    # Adding Residual Blocks
    for _ in range(num_blocks):
        x = residual_block(x, 64, 5)
        x = residual_block(x, 64, 5, downsample=False)  # Avoid downsampling if it changes the output sequence length

    # Final Conv1D to adjust the number of output channels to match the number of classes (for binary: 1)
    x = Conv1D(1, kernel_size=3, padding='same')(x)
    x = BatchNormalization()(x)
    x = ReLU()(x)

    # Using a sigmoid activation function to get outputs in the range of 0 to 1 for binary classification at each position
    outputs = Conv1D(1, kernel_size=1, activation='sigmoid', padding='same')(x)  # This ensures output shape matches the target shape (None, sequence_length, 1)

    model = Model(inputs=inputs, outputs=outputs)
    return model


def gen_conv1d_model(input_size:int):
    '''
    Defines Conv1D architecture model
    '''
    model = Sequential([
        Conv1D(filters=64, kernel_size=3, activation='relu', input_shape=(input_size, 1)),
        MaxPooling1D(pool_size=2),
        Conv1D(filters=128, kernel_size=3, activation='relu'),
        MaxPooling1D(pool_size=2),
        Flatten(),
        Dense(64, activation='relu'),
        Dropout(0.5),
        Dense(1, activation='sigmoid')  # Change to 'softmax' and adjust the last layer units if multi-class
    ])
    return model

def gen_dense_model(input_size:int, dropout=0.2):
    '''
    Defines a dense (i.e. fully connected) model with three hidden layers and dropout for regularization.
    
    Args:
    - input_size (int): Number of features in the input dataset.
    
    Returns:
    - keras.Model: A compiled dense neural network model.
    '''
    model = keras.Sequential([
        keras.layers.Dense(input_size//2, activation='relu', input_shape=(input_size,)),
        keras.layers.Dropout(dropout),  #  regularization
        keras.layers.Dense(input_size//4, activation='relu'),
        # keras.layers.Dropout(dropout), 
        # keras.layers.Dense(150, activation='relu'),
        # keras.layers.Dropout(0.2), 
        # keras.layers.Dense(input_size//4, activation='relu'),
        keras.layers.Dense(1, activation='sigmoid')  # Use 'sigmoid' for binary classification
    ])
    return model

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

def show_correlation(df, columns=None, plot_font_size=4):
    # Calculate correlation matrix
    if isinstance(df, np.ndarray):
        df = pd.DataFrame(df)
        if columns is not None:
            df.columns=columns
    correlation_matrix = df.select_dtypes(include=[np.number]).corr()
    plt.figure(figsize=(8, 6))
    ax=sns.heatmap(correlation_matrix, annot=True, cmap='coolwarm', fmt=".2f", annot_kws={'size': plot_font_size})
    ax.tick_params(axis='both', which='major', labelsize=4)  # Reducing the font size of the tick labels
    # Set the ticks to be at every data point
    # ax.set_xticks(np.arange(df.shape[1]) + 0.5, minor=False)
    # ax.set_yticks(np.arange(df.shape[0]) + 0.5, minor=False)

    # Set tick labels to be the default or any other labels you want
    # ax.set_xticklabels(np.arange(1, df.shape[1] + 1))
    # ax.set_yticklabels(np.arange(1, df.shape[1] + 1))
    plt.title('Correlation Matrix of Numeric Features')
    plt.show()

def main(df_train, df_test, outbase='.',
    mode='train',
    model_name='dense',
    epochs=10, # 30,
    window_size=17,
    # prep_data=False,
    drop_columns=['Unnamed: 0'], # 'ResNo'],
    early_stopping = EarlyStopping( # Setup early stopping
        monitor='val_loss',  # Monitor validation loss
        min_delta=0.01,      # Minimum change to qualify as an improvement
        patience=5,         # Stop after this many epochs without improvement
        verbose=1,           # Print messages when stopping
        mode='min',          # we are minimizing
        restore_best_weights=True  # Restore model weights from the epoch with the best value of the monitored quantity
    ),
    use_smote=True
    
    ):
    '''
    We are not relying on `test_train_split` because we
    want to avoid the same accession codes ending up in
    both test and train data. Instead the split is managed
    by the script `data_disprotstructure.py`

    Args:
    df_train: DataFrame or filenamewith one residue per row, last column is target variable
    df_test: dataframe or filename of test data

    model(str): 'dense', 'conv1d', 'resnet'
    '''
    df_train = read_csv_if_string(df_train)
    df_test = read_csv_if_string(df_test)
    df_train = df_train.drop(columns=drop_columns, errors='ignore')
    df_test = df_test.drop(columns=drop_columns, errors='ignore')
    print("Column names:", df_train.columns)
    if mode == 'corr':
        show_correlation(df_train)

    if model_name in ['conv1d', 'resnet']:
        assert isinstance(df_train,pd.DataFrame)
        X_train,y_train=prepare_data_for_per_position_classification(df_train, "Accession", "structured")
        X_test, y_test=prepare_data_for_per_position_classification(df_test, "Accession", "structured")
    elif model_name == 'dense':
        # Remove text columns
        # df_train_numeric = df_train.select_dtypes(include=[np.number])
        # Separate features and target
        # X = df_numeric.iloc[:, :-1]  # all columns except the last one
        train_data_result = prepare_data_sliding_window(df_train, window_size=window_size)
        X_train = train_data_result['X']
        y_train = train_data_result['y']
        X_columns = train_data_result['columns']
        test_data_result = prepare_data_sliding_window(df_test, window_size=window_size)
        X_test = test_data_result['X']
        y_test = test_data_result['y']  
        # print(X_train[:10,:])
        # assert False
        # y = df_numeric.iloc[:, -1]   # the last column
    else:
        raise ValueError("Error: supported model modes are resnet, dense, conv1d")
    if mode == 'corr':
        show_correlation(np.concatenate([X_train, y_train.reshape(-1, 1)], axis=1), X_columns + ['structured'])
    # print("X_train rows:")
    # print(X_columns)
    # print(X_train[:5,:])
    # print(y_train[:5])

    if use_smote:
        from imblearn.over_sampling import SMOTE
        smote = SMOTE()
        X_train, y_train = smote.fit_resample(X_train, y_train)

    # Normalize the features
    scaler = StandardScaler()
    if model_name == 'dense':
        X_train_scaled = scaler.fit_transform(X_train)
        X_test_scaled = scaler.transform(X_test)
    else:
        X_train_reshaped = X_train.reshape(-1, X_train.shape[-1])  # Reshape to scale
        X_train_scaled = scaler.fit_transform(X_train_reshaped)
        X_train_scaled = X_train_scaled.reshape(X_train.shape)  # Reshape back to original shape
        X_test_reshaped = X_test.reshape(-1, X_test.shape[-1])
        X_test_scaled = scaler.transform(X_test_reshaped)
        X_test_scaled = X_test_scaled.reshape(X_test.shape)
    input_size=X_train_scaled.shape[1]
    print("Training feature data shape:", X_train_scaled.shape, input_size)
    print("Training label data shape:", y_train.shape)
    if model_name == 'dense':
        model = gen_dense_model(input_size=input_size)
    elif model_name == 'conv1d':
        model = gen_conv1d_model_per_position(X_train_scaled)#input_size=input_size)
    elif model_name == 'resnet':
        input_shape = X_train_scaled.shape[1:]
        model = build_resnet(input_shape, num_classes=1, num_blocks=2)
    else:
        raise ValueError("Supported model types are dense, conv1d, resnet")

    # Compile the model with accuracy and precision as additional metrics
    model.compile(optimizer='adam', 
                loss='binary_crossentropy',  # Use 'binary_crossentropy' for binary classification
                metrics=['accuracy','auc', tf.keras.metrics.Precision()])
    print(f"{model_name} model summary:")
    model.summary()
    # Train the model
    # model.fit(X_train_scaled, y_train, epochs=epochs, validation_split=0.1)

    history = model.fit(
        X_train_scaled, y_train,
        validation_split=0.1,
        # validation_data=(X_val, y_val),
        epochs=epochs,
        callbacks=[early_stopping]
    )

    # Evaluate the model on test data
    model.evaluate(X_test_scaled, y_test)
    model_outfile = os.path.join(MODEL_DIR, outbase + '.keras')
    scaler_outfile = model_outfile.replace(".keras", '_scaler.pkl')
    history_outfile = model_outfile.replace(".keras","_history.pkl")
    if not os.path.exists(MODEL_DIR):
        print("Creating model output directory", MODEL_DIR)
        os.makedirs(MODEL_DIR)
    model.save(model_outfile)
    # Saving the scaler object to a pickle file
    with open(scaler_outfile, 'wb') as f:
        print("Saving scaler to", scaler_outfile)
        pickle.dump(scaler, f)
    # Save the history to a pickle file
    with open(history_outfile, 'wb') as file_pi:
        print("Saving pickle file to", history_outfile)
        pickle.dump(history.history, file_pi)

    # XGBOOST
    # fit baseline model:
    # Initialize and train XGBoost model
    xgb_model = xgb.XGBClassifier(objective='binary:logistic', 
        eval_metric='logloss')
    xgb_model.fit(X_train_scaled, y_train)
    # Make predictions with XGBoost
    y_pred_xgb = xgb_model.predict(X_test_scaled)
    y_pred_proba_xgb = xgb_model.predict_proba(X_test_scaled)[:, 1]

    # Evaluate XGBoost model
    accuracy_xgb = accuracy_score(y_test, y_pred_xgb)
    auc_xgb = roc_auc_score(y_test, y_pred_proba_xgb)
    print(f"XGBoost Accuracy: {accuracy_xgb}, AUC: {auc_xgb}")
    xgb_model_outfile = os.path.join(MODEL_DIR, outbase + '_xgb_model.json')
    print("Saving XGBoost model to", xgb_model_outfile)
    xgb_model.save_model(xgb_model_outfile)


def run_main_model75(data_file=os.path.join(TABLE_DIR, "disorder_75.tsv")):
    main(df=data_file,outbase='disorder_disorder_75')

def run_main_model839(data_file=os.path.join(TABLE_DIR, "disorder_839.tsv")):
    main(df=data_file,outbase='disorder_disorder_839')

def run_main_model_protein_binding(data_file=os.path.join(TABLE_DIR, "disorder_protein binding_541.tsv")):
    main(df=data_file, outbase='disorder_protein_binding_541')

def run_main_model_disorder_to_order(data_file=os.path.join(TABLE_DIR, "disorder_disorder to order_396.tsv")):
    main(df=data_file, outbase='disorder_disorder to order_396')

def run_main_train_all(
    input_dir=os.path.join(PROCESSED_DIR,'disorder'),
    output_dir=f"{MODEL_DIR}",
    model_name='dense',
    window_size=17,
    mode='train',
    regions=['flexible linker/spacer',
    'molecular adaptor activity', 'molecular function regulator',
    'phosphorylation display site',
    'lipid binding',
    'ion binding',
    'pre-molten globule',
    'nucleic acid binding','disorder', 'protein binding', 'disorder to order']):
    
    for region in regions:
        region_name = region.replace("/","_").replace(" ", "_")
        filenames = glob.glob(f"{input_dir}/disorder_{region_name}_*_train.tsv")
        print(f"Train file names for region {region}:", filenames)
        for train_file in filenames:
            test_file = train_file.replace("_train", "_test")
            if not os.path.exists(test_file):
                print("Warning: could not find test file:",
                    test_file)
            # data_file=os.path.join(PROCESSED_DIR,'disorder', f"disorder_disorder to order_396.tsv")
            outbase = os.path.join(output_dir, f'{basename(train_file.replace("_train.tsv", ""))}_{model_name}_w{window_size}')
            print(f"Starting training model from {train_file} and output basename {outbase}")
            main(df_train=train_file, df_test=test_file, outbase=outbase, model_name=model_name,
            window_size=window_size, mode=mode)

def run_all_tests():
    test_prepare_data_sliding_window()

if __name__ == "__main__":
    mode = 'train'
    if len(sys.argv) > 1:
        if sys.argv[1] == 'test':
            run_all_tests()
            sys.exit(0)
        mode = sys.argv[1]
    run_main_train_all(model_name='dense', window_size=5,   
        mode=mode) # regions=['disorder'],
    # run_main_model_disorder_to_order()
    # run_main_model_protein_binding()
    # run_main_model839()