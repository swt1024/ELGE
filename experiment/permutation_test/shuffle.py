import pandas as pd
import numpy as np
import random
import os

# File paths
esslnc_path = '../../data/benchmark/mouse/ess_lpi.csv'
nonesslnc_path = '../../data/benchmark/mouse/noness_lpi.csv'

# Load data
esslnc = pd.read_csv(esslnc_path)
nonesslnc = pd.read_csv(nonesslnc_path)

esslnc_id = set(esslnc['lncRNA_ID'])
nonesslnc_id = set(nonesslnc['lncRNA_ID'])

# Prepare data arrays
ids_positive = list(esslnc_id)  # Convert set to list
ids_negative = list(nonesslnc_id)  # Convert set to list

# Combine datasets
y_all = np.hstack((np.ones(len(ids_positive)), np.zeros(len(ids_negative))))
ids_all = np.hstack((ids_positive, ids_negative))

# Create directory to store shuffled datasets
output_dir = './mouse_shuffled'
os.makedirs(output_dir, exist_ok=True)

# Save the original data as shuffle_0.csv before shuffling
original_data = pd.DataFrame({
    'ID': ids_all,
    'Label': y_all
})
original_file_path = os.path.join(output_dir, 'shuffled_0.csv')
original_data.to_csv(original_file_path, index=False)
print(f"Saved original data to {original_file_path}")

# Shuffle and save 1000 times
for shuffle_num in range(1000):
    # Shuffle the labels (y_all) for the given iteration
    shuffled_y_all = random.sample(list(y_all), len(y_all))  # Convert y_all to list if needed

    # Store the shuffled data into a new CSV file for each iteration
    shuffled_file_path = os.path.join(output_dir, f'shuffled_{shuffle_num+1}.csv')
    shuffled_data = pd.DataFrame({
        'ID': ids_all,
        'Label': shuffled_y_all
    })
    shuffled_data.to_csv(shuffled_file_path, index=False)

    print(f"Saved shuffled data {shuffle_num + 1} to {shuffled_file_path}")
