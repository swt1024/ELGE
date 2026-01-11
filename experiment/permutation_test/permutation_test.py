import pandas as pd
import numpy as np

model_names = ['svm', 'mlp']

# Function to compute the p-value based on the observed performance and random performances
def compute_p_value(obs_performance, random_performances):
	p_values = {}  # Dictionary to store the p-values for each metric
	for idx, metric in enumerate(['sen','spe','ppv','f1','acc','mcc','roc_auc','pr_auc']):
		metric_values = random_performances[:, idx]
		p_value = np.sum(metric_values >= obs_performance[idx]) / len(metric_values)
		p_values[metric] = p_value
	return p_values

# List to collect all results


for species in ['human','mouse']:
	results = []
	for model_name in model_names:
		if species == 'human':
			tissues = ['heart','lung','stomach']
		else :
			tissues = ['heart','lung','brain']
		for tissue in tissues:
			# ------------------------ 1. Read CSV File ------------------------
			file_path = f'{model_name}_{species}_{tissue}_shuffled_performance.csv'
			df = pd.read_csv(file_path)

			# ------------------------ 2. Extract Data ------------------------
			obs_performance = df[df['Shuffle_File'] == 'shuffled_0.csv'][['sen','spe','ppv','f1','acc','mcc','roc_auc','pr_auc']].values.flatten()
			random_performances = df[df['Shuffle_File'] != 'shuffled_0.csv'][['sen','spe','ppv','f1','acc','mcc','roc_auc','pr_auc']].values

			# ------------------------ 3. Compute P-Value ------------------------
			p_values = compute_p_value(obs_performance, random_performances)

			# ------------------------ 4. Store Results ------------------------
			result_row = {
				'model': model_name,
				'tissue': tissue
			}
			result_row.update(p_values)  # Add all metrics p-values
			results.append(result_row)

	# Convert results list to DataFrame
	pval_df = pd.DataFrame(results)

	# Save to CSV
	pval_df.to_csv(f'{species}_p_values.csv', index=False)

	print(f"P-values saved to {species}_p_values.csv")
