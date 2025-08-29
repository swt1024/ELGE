import pandas as pd
import numpy as np
from sklearn.metrics import (
    confusion_matrix, f1_score, matthews_corrcoef, accuracy_score,
)

# === Step 1: read data ===

gic_file = '../../results/mouse/mouse_GIC_score.csv'
essential_file = '../../data/benchMarking/mouse/ess_lpi.csv'
nonessential_file = '../../data/benchMarking/mouse/noness_lpi.csv'

gic_df = pd.read_csv(gic_file)
ess_ids = pd.read_csv(essential_file, header=None)[0].tolist()
noness_ids = pd.read_csv(nonessential_file, header=None)[0].tolist()

# retain only essential and non-essential lncRNAs
gic_df = gic_df[gic_df['lncRNA_ID'].isin(ess_ids + noness_ids)]

# label the data
# 1: essential, 0: non-essential
gic_df['Label'] = gic_df['lncRNA_ID'].apply(lambda x: 1 if x in ess_ids else (0 if x in noness_ids else np.nan))
gic_df.dropna(inplace=True)

y_true = gic_df['Label'].astype(int).values
y_scores = gic_df['Score'].values

# === Step 2: calculate metrics ===
def compute_metrics(y_true, y_pred):
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
    sen = tp / (tp + fn) if (tp + fn) > 0 else 0
    spe = tn / (tn + fp) if (tn + fp) > 0 else 0
    ppv = tp / (tp + fp) if (tp + fp) > 0 else 0
    acc = accuracy_score(y_true, y_pred)
    f1 = f1_score(y_true, y_pred)
    mcc = matthews_corrcoef(y_true, y_pred)
    return {
        'Sensitivity': sen,
        'Specificity': spe,
        'PPV': ppv,
        'Accuracy': acc,
        'F1': f1,
        'MCC': mcc
    }

# === Step 3: Use 0.5 as threshold ===
threshold = 0.5
y_pred = (y_scores >= threshold).astype(int)
metrics = compute_metrics(y_true, y_pred)

# === Step 4: print metrics ===
print(f"Threshold: {threshold}")
print(metrics)
