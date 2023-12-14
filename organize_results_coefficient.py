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

p_table = pd.DataFrame()

for i, file in enumerate(result_files):
    cell_type = file.split('.')[0].split('_')[-1]
    if cell_type in p_table.index:
        raise Exception('Cell type repeated in file {}'.format(file))
    with open(join(result_path, file), 'r') as f:
        header = f.readline().split('\t')
        output = f.readline().split('\t')
        if output[0] != 'L2_0':
            raise Exception('L2_0 is not the second line of results file')
        zvalue = output[header.index('Coefficient_z-score\n')]
        transfomredP = (1 - stats.norm.cdf(abs(float(zvalue))))
        p_table.loc[cell_type, 'p_value'] = transfomredP

        for col in ['Prop._SNPs', 'Prop._h2', 'Enrichment']:
            p_table.loc[cell_type, col] = float(output[header.index(col)])

p_table.to_csv('LDSC_results.csv')
