import pandas as pd
import numpy as np
import os
import argparse
from sklearn.model_selection import KFold, LeaveOneOut
from sklearn.metrics import confusion_matrix, matthews_corrcoef, roc_curve, auc, precision_recall_curve
from sklearn.svm import LinearSVC
from sklearn.preprocessing import MinMaxScaler

# Set up command line arguments
parser = argparse.ArgumentParser(description='Process species and tissue type.')
parser.add_argument('species', type=str, help='The species to process (e.g., mouse, human)')
parser.add_argument('tissue', type=str, help='The tissue type to process (e.g., heart, liver)')
args = parser.parse_args()

# Assign arguments to variables
species = args.species
tissue = args.tissue

# File paths
X_path = f'../../HinSAGE/{species}/lncRNA_embeddings_{tissue}.csv'
esslnc_path = f'../../data/benchmark/{species}/ess_lpi.csv'
nonesslnc_path = f'../../data/benchmark/{species}/noness_lpi.csv'

# Load data
X = pd.read_csv(X_path, index_col=0, header=None)
esslnc = pd.read_csv(esslnc_path)
nonesslnc = pd.read_csv(nonesslnc_path)

# Get the lncRNA IDs for positive and negative samples
esslnc_id = set(esslnc['lncRNA_ID'])
nonesslnc_id = set(nonesslnc['lncRNA_ID'])

# Prepare data arrays based on the IDs
X_positive = X[X.index.isin(esslnc_id)].values
X_negative = X[X.index.isin(nonesslnc_id)].values
ids_positive = X[X.index.isin(esslnc_id)].index
ids_negative = X[X.index.isin(nonesslnc_id)].index

# Combine datasets
X_all = np.vstack((X_positive, X_negative))
ids_all = np.hstack((ids_positive, ids_negative))

# Ensure the output directory for metrics exists
save_path = f'svm_{species}_{tissue}_shuffled_performance.csv'

# Create the CSV file and write the header only once
if not os.path.exists(save_path):
    header = ['Shuffle_File', 'sen', 'spe', 'ppv', 'f1', 'acc', 'mcc', 'roc_auc', 'pr_auc']
    metrics_df = pd.DataFrame(columns=header)
    metrics_df.to_csv(save_path, mode='w', header=True, index=False)

# Loop through shuffled datasets (0 to 1000)
for x in range(0, 1001):
    shuffle_file = f"shuffled_{x}.csv"
    
    # Load the shuffled dataset and check if it has been processed before
    shuffled_data = pd.read_csv(os.path.join(f'./{species}_shuffled', shuffle_file))

    # Get shuffled labels
    shuffled_y_all = shuffled_data['Label'].values

    # Initialize KFold and LeaveOneOut
    if species == 'mouse':
        cv = LeaveOneOut()
        C = 10
    else:
        cv = KFold(n_splits=10, shuffle=True, random_state=42)
        C = 10

    # Initialize lists to store true labels and decision scores for performance evaluation
    all_true_labels = []
    all_decision_scores = []
    # Prepare DataFrame to save experimental records
    experiment_records = pd.DataFrame()

    # KFold cross-validation
    for fold, (train_index, test_index) in enumerate(cv.split(X_all)):
        X_train, X_test = X_all[train_index], X_all[test_index]
        y_train, y_test = shuffled_y_all[train_index], shuffled_y_all[test_index]
        ids_train, ids_test = ids_all[train_index], ids_all[test_index]

        # Apply MinMaxScaler
        scaler = MinMaxScaler()
        X_train_scaled = scaler.fit_transform(X_train)
        X_test_scaled = scaler.transform(X_test)

        # Train LinearSVC model
        svm = LinearSVC(C=C, dual=False)
        svm.fit(X_train_scaled, y_train)

        # Get decision scores and predictions
        decision_scores = svm.decision_function(X_test_scaled)
        predictions = (decision_scores >= 0).astype(int)

        # Store true labels and decision scores
        all_true_labels.extend(y_test)
        all_decision_scores.extend(decision_scores)

        fold_data = {
                'Fold': fold + 1,
                'Train_IDs': [list(ids_train)],
                'Train_Labels': [list(y_train)],
                'Test_IDs': [list(ids_test)],
                'Test_Labels': [list(y_test)],
                'Predictions': [list(predictions)],
                'Decision_Scores': [list(decision_scores)]
            }
        fold_df = pd.DataFrame(fold_data)
        experiment_records = pd.concat([experiment_records, fold_df], ignore_index=True)

    # Convert lists to arrays for performance evaluation
    all_true_labels = np.array(all_true_labels)
    all_decision_scores = np.array(all_decision_scores)

    # Compute confusion matrix using threshold at 0
    tn, fp, fn, tp = confusion_matrix(all_true_labels, (all_decision_scores >= 0).astype(int)).ravel()

    # Compute performance metrics
    sensitivity = tp / (tp + fn)
    specificity = tn / (tn + fp)
    ppv = tp / (tp + fp) if (tp + fp) > 0 else 0
    accuracy = (tp + tn) / (tp + tn + fp + fn)
    f1_score = 2 * (ppv * sensitivity) / (ppv + sensitivity) if (ppv + sensitivity) > 0 else 0
    mcc = matthews_corrcoef(all_true_labels, (all_decision_scores >= 0).astype(int))

    # Compute ROC curve and PR curve data
    fpr, tpr, _ = roc_curve(all_true_labels, all_decision_scores)
    roc_auc = auc(fpr, tpr)

    precision, recall, _ = precision_recall_curve(all_true_labels, all_decision_scores)
    pr_auc = auc(recall, precision)

    # Store metrics for the current shuffle iteration
    metrics_row = {
        'Shuffle_File': shuffle_file,
        'Sensitivity(Recall)': sensitivity,
        'Specificity': specificity,
        'PPV(Precision)': ppv,
        'F1_Score': f1_score,
        'Accuracy': accuracy,
        'MCC': mcc,
        'roc_auc': roc_auc,
        'pr_auc': pr_auc
    }

    # Append the current metrics row to the CSV file immediately
    metrics_df = pd.DataFrame([metrics_row])  # Create a DataFrame for the current row
    metrics_df.to_csv(save_path, mode='a', header=False, index=False)

    # Save experiment records for the current shuffle iteration
    experiment_records.to_csv(f'./experiment_details/{species}/svm/{tissue}/{tissue}_experiment_details_{x}.csv', index=False)

    print(f"Processed and saved results for {shuffle_file}")

