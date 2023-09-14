

import os
import pandas as pd
import yaml

data_dir = os.path.join('Data')

# Create data directory structure
experiments = ['experiment1', 'experiment2' ]

urls=['','']

n=range(2)
for i in n:
    exp_dir = os.path.join(data_dir, experiments[i])
    os.makedirs(exp_dir, exist_ok=True)

    # Create metadata file (e.g., YAML)
    #add more relevant information: source, organism, platform,, authors, ....
    metadata = {
        'Experiment': experiments[i],
        'Data_Type': 'SingleCell',
        'Data_URL': urls[i]
    }
    
    metadata_file = os.path.join(exp_dir, 'metadata.yaml')
    with open(metadata_file, 'w') as f:
        yaml.dump(metadata, f)

# Create a sample metadata DataFrame
metadata_df = pd.DataFrame({
    'Experiment': experiments,
    'Data_Type': 'SingleCell',
    'Data_URL': urls
})
    
# Save metadata to CSV
metadata_df.to_csv(os.path.join(data_dir, 'metadata.csv'), index=False)

 
	


