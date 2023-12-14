import pandas as pd
from os import listdir, system, getenv
from os.path import join, isdir
import argparse

parser = argparse.ArgumentParser()
# parser.add_argument('bedfile_dir')
parser.add_argument('all_peak_excel_file')
args = parser.parse_args()

threshold = 3

all_type = pd.read_excel(args.all_peak_excel_file)

for c in all_type.columns[3:]:
    all_type.loc[:, ['chr', 'start', 'end']][all_type[c] >= threshold].to_csv(c + '.bed', sep='\t', header=False, index=False)
