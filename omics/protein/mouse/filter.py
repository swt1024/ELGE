import os
import pandas as pd

def filter_rpt_files(folder_path, filter_conditions_file, output_file):
    # Read filter conditions
    filters_df = pd.read_csv(filter_conditions_file, dtype=str)

    # Get all expression data files (.rpt) in the folder
    rpt_files = [os.path.join(folder_path, f) for f in os.listdir(folder_path) if f.endswith('.rpt')]

    # Store rows that meet the conditions
    filtered_data = []

    for rpt_file in rpt_files:
        try:
            # Read .rpt file using '|' as the delimiter
            df = pd.read_csv(rpt_file, sep='|', dtype=str)

            # Ensure all required columns are present
            required_columns = ["Anatomical Structure", "Theiler Stage", "Age", "Strain"]
            if all(col in df.columns for col in required_columns):
                # Get filter conditions specific to the current file
                conditions = filters_df[filters_df['File'] == os.path.basename(rpt_file)]
                if conditions.empty:
                    continue

                # Turn conditions into a set of tuples for easy filtering
                conditions_set = set([tuple(x) for x in conditions.iloc[:, 1:].values])
                df_filtered = df[df[required_columns].apply(tuple, axis=1).isin(conditions_set)]

                # Insert a 'File' column to indicate the source file if not empty
                if not df_filtered.empty:
                    df_filtered.insert(0, "File", os.path.basename(rpt_file))
                    filtered_data.append(df_filtered)
            else:
                print(f"File {rpt_file} is missing required columns.")
        except Exception as e:
            print(f"Error processing {rpt_file}: {e}")

    # Concatenate and save
    if filtered_data:
        final_df = pd.concat(filtered_data, ignore_index=True)
        final_df.to_csv(output_file, sep=',', index=False, encoding='utf-8')
        print(f"Filtered data saved to {output_file}")
    else:
        print("No valid data found.")

# Example usage
folder_path = "MGI_exp"  # Folder containing all .rpt files with expression data
filter_conditions_file = "MGI_filter_conditions.csv"  # CSV file with filter conditions
output_file = "filtered_exp.csv"

filter_rpt_files(folder_path, filter_conditions_file, output_file)
