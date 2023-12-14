import pandas as pd
import os
from os import listdir
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument('threshold', type=int)
parser.add_argument('bed_path')
args = parser.parse_args()
bed_path = args.bed_path
threshold = args.threshold


counter = 0
file_by_ct = defaultdict(list)  # file by cell type
bed_files = [f for f in listdir(bed_path) if f.endswith('.bed')]
for file in bed_files:
    cell_type = file.split('.')[0].split('_')[-1]
    file_by_ct[cell_type].append(file)

for ct, files in file_by_ct.items():
    count = 0
    out = []
    if len(files) != 10:
        raise Exception('Expect ten deciles. Now get {}.'.format(len(files)))
    files.sort(key=lambda x: int(x.split('.')[1]))
    print(files)
    for file in files:
        peaks = pd.read_csv(os.path.join(bed_path,file), sep='\t', names=['chr', 'start', 'end', 'overlap'])
        decile = int(file.split('.')[1])
        if decile < 1 or decile > 10:
            raise Exception('decile number is assumed to be integer between 1 and 10. Now get {}'.format(decile))

        peaks['decile'] = decile

        if count + peaks['overlap'].sum() <= threshold:
            out.append(peaks.sort_values('overlap', ascending=False))
            count += peaks['overlap'].sum()
            if count == threshold:
                break
        else:
            # partial selection from within this decile
            peaks.sort_values('overlap', inplace=True, ascending=False)
            idx = peaks['overlap'].cumsum().searchsorted(threshold - count, side='right')
            if peaks['overlap'].iloc[:idx].sum() == threshold - count:
                out.append(peaks.iloc[0:idx])
            else:
                out.append(peaks.iloc[0:idx + 1])

            break
    out = pd.concat(out, ignore_index=True, axis=0)
    out.to_csv('overlap_{}_{}.bed'.format(threshold, ct), sep='\t', index=False, header=False)
