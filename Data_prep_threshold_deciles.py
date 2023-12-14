# save it as process_excel_to_bed.py
import pandas as pd
import numpy as np
from os.path import join, isdir
from math import floor
import argparse

def process_excel_to_bed(excel_file, output_directory, threshold):
    # Read the Excel file
    all_type = pd.read_excel(excel_file)

    # Create the output directory if it doesn't exist
    if not isdir(output_directory):
        system(f'mkdir {output_directory}')

    percentile_list = np.arange(0, 1.1, 0.1)  # [0, 0.1, 0.2, ...1]

    # Iterate through each column in the Excel file
    for c in all_type.columns[3:]:
        # Select rows where the value is above the threshold and sort by the column in descending order
        above_threshold = all_type.loc[:, ['chr', 'start', 'end']][all_type[c] >= threshold].sort_values(by=[c], ascending=False)

        nums = above_threshold.shape[0]
        
        # Save subsets based on percentiles
        for i in range(1, len(percentile_list)):
            percentile_file = join(output_directory, f'{c}.{int(100*percentile_list[i]):03d}.bed')
            subset = above_threshold.iloc[floor(nums * percentile_list[i-1]):floor(nums * percentile_list[i])]
            subset.to_csv(percentile_file, sep='\t', header=False, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process Excel file to generate BED files based on a threshold.")
    parser.add_argument("excel_file", help="Path to the input Excel file")
    parser.add_argument("output_directory", help="Path to the output directory")
    parser.add_argument("threshold", type=float, help="Threshold value for filtering")
    
    args = parser.parse_args()

    process_excel_to_bed(args.excel_file, args.output_directory, args.threshold)
