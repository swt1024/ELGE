import pandas as pd
import requests
import os

def download_files(csv_path, download_folder):
    # Ensure the download folder exists
    if not os.path.exists(download_folder):
        os.makedirs(download_folder)

    # Read the CSV file
    data = pd.read_csv(csv_path)
    
    # Extract and deduplicate filenames from the 'File' column
    filenames = data['File'].drop_duplicates()
    
    # URL template
    url_template = "https://www.informatics.jax.org/downloads/reports/gxdrnaseq/{}.gz"

    # Loop through the filenames and download each file
    for filename in filenames:
        url = url_template.format(filename)
        response = requests.get(url)
        
        if response.status_code == 200:
            # Define the path to save the file
            file_path = os.path.join(download_folder, f"{filename}.gz")
            # Write the content to a file
            with open(file_path, 'wb') as f:
                f.write(response.content)
            print(f"Downloaded {filename}")
        else:
            print(f"Failed to download {filename}: Status code {response.status_code}")


file_path = 'MGI_filter_conditions.csv'
download_folder = 'MGI_exp'
# Call the function to download files
download_files(file_path, download_folder)
