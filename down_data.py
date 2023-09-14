import os
import pandas as pd
import requests
import yaml
from urllib.parse import urlparse, parse_qs
# Define the project directory
#project_dir = 'ASTRA'
data_dir = ('ASTRA_Data')


# Create a function to download data
def download_data(url, save_dir):
    response = requests.get(url)
    if response.status_code == 200:
        # Extract the filename from the URL
        

        with open(save_path, 'wb') as f:
            f.write(response.content)
        return file_name  # Return the downloaded file name
    else:
        print(f"Failed to download data from {url}")
        return None
    

# Load metadata
metadata_file = os.path.join(data_dir, 'metadata.csv')
metadata_df = pd.read_csv(metadata_file)

# Iterate through metadata and download data
for index, row in metadata_df.iterrows():
    exp_dir = os.path.join(data_dir, row['Experiment'])
    os.makedirs(exp_dir, exist_ok=True)
    
    # Check if Data_URL is not empty (i.e., the data is available for download)
    if row['Data_URL']:
        data_url = row['Data_URL']
        
        #parsed_url = urlparse(data_url)
        parsed_url=parse_qs(urlparse(data_url).query)
       # file_name = os.path.basename(parsed_url.query)
        file_name=os.path.basename(parsed_url["file"][0])
        save_path = os.path.join(exp_dir, file_name)
        
        # Download data and get the actual downloaded file name
       # downloaded_file_name = download_data(data_url, exp_dir)
        downloaded_file_name = download_data(data_url, save_path)
        if downloaded_file_name:
            # Update the Data_URL in metadata.yaml
            metadata_path = os.path.join(exp_dir, 'metadata.yaml')
            with open(metadata_path, 'r') as f:
                metadata = yaml.safe_load(f)
                metadata['Data_URL'] = downloaded_file_name  # Update with the actual downloaded file name
            with open(metadata_path, 'w') as f:
                yaml.dump(metadata, f)

# Once the data is downloaded and named, you can proceed with your analysis or preprocessing steps.




