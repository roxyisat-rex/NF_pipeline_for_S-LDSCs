from os import listdir
from os.path import join
import pandas as pd
import argparse
import scipy.stats as stats

parser = argparse.ArgumentParser()
parser.add_argument('res_path')
args = parser.parse_args()
result_path = args.res_path

result_files = [f for f in listdir(result_path) if f.endswith('.results')]

p_table = pd.DataFrame(columns=['Cell_Type', 'Quantile', 'Prop._h2', 'Enrichment', 'Prop._SNPs'])

for i, file in enumerate(result_files):
    # Extract cell type and quantile from the filename
    parts = file.split('.')
    cell_type = parts[0].split('_')[-1]
    quantile = parts[-2] if len(parts) > 2 else 'NA'

    with open(join(result_path, file), 'r') as f:
        header = f.readline().split('\t')
        output = f.readline().split('\t')
        if output[0] != 'L2_0':
            raise Exception('L2_0 is not the second line of results file')

        zvalue = output[header.index('Coefficient_z-score\n')]
        transformedP = (1 - stats.norm.cdf(abs(float(zvalue))))

        # Append a new row to the DataFrame for each quantile
        p_table = p_table.append({'Cell_Type': cell_type, 'Quantile': quantile,
                                  'Prop._h2': float(output[header.index('Prop._h2')]),
                                  'Enrichment': float(output[header.index('Enrichment')]),
                                  'Prop._SNPs': float(output[header.index('Prop._SNPs')]),
                                  'Transformed_P': transformedP}, ignore_index=True)

# Sort the DataFrame by Cell_Type and Quantile in descending order
p_table = p_table.sort_values(by=['Cell_Type', 'Quantile'], ascending=[True, False])

p_table.to_csv('LDSC_results_deciles.csv', index=False)
