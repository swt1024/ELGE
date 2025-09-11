import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import rcParams
rcParams['font.family'] = 'Times New Roman'

for species in ['mouse','human']:

	# File paths for Biological Process (BP), Cellular Component (CC), and Molecular Function (MF)
	bp_file = f'GOTERM_BP_{species}.csv' 
	cc_file = f'GOTERM_CC_{species}.csv'  
	mf_file = f'GOTERM_MF_{species}.csv'  

	# Read the three files into separate DataFrames
	bp_df = pd.read_csv(bp_file)
	cc_df = pd.read_csv(cc_file)
	mf_df = pd.read_csv(mf_file)

	# Add a 'Category' column to each DataFrame for identification
	bp_df['Category'] = 'Biological Process'
	cc_df['Category'] = 'Cellular Component'
	mf_df['Category'] = 'Molecular Function'

	# Combine all three dataframes into one
	df = pd.concat([bp_df, cc_df, mf_df], ignore_index=True)

	# Filter out rows with FDR >= 0.05 (non-significant GO terms)
	significant_df = df[df['FDR'] < 0.05]
	significant_df = significant_df[significant_df['Fold Enrichment'] > 1]

	# For each category, get the top 5 GO terms by sorting based on p-value
	top_bp = significant_df[significant_df['Category'] == 'Biological Process'].sort_values('PValue').head(5)
	top_cc = significant_df[significant_df['Category'] == 'Cellular Component'].sort_values('PValue').head(5)
	top_mf = significant_df[significant_df['Category'] == 'Molecular Function'].sort_values('PValue').head(5)

	# Combine the top 5 GO terms from each category
	top_terms_df = pd.concat([top_bp, top_cc, top_mf], ignore_index=True)

	# Set up the color palette for different categories
	palette = {"Biological Process": "#EECA40", 
			"Cellular Component": "#FD763F", 
			"Molecular Function": "#23BAC5"}

	# Create the barplot
	plt.figure(figsize=(10, 5))
	sns.barplot(x='Count', y='Term', hue='Category', data=top_terms_df, palette=palette)

	# Set title and labels
	plt.title(f"GO Terms Enrichment ({species})", fontsize=16)
	plt.xlabel("LncRNA Gene Count", fontsize=14)
	plt.ylabel("GO Term", fontsize=14)

	# Move the legend outside the plot
	plt.legend(
		title='Category', 
		bbox_to_anchor=(1.02, 1), 
		loc='upper left', 
		borderaxespad=0
	)

	# Adjust layout to prevent clipping
	plt.tight_layout(rect=[0, 0, 0.99, 1])

	# Save the plot
	plt.savefig(f'go_terms_enrichment_{species}.svg', dpi=300)
