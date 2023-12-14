import pandas as pd
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('bim_path')
parser.add_argument('bedfile')
args = parser.parse_args()
bedfile = args.bedfile
bim_path = args.bim_path


bim_stem = '1000G.EUR.hg38'

bim = {}
for i in range(1, 23):
    temp = pd.read_csv(os.path.join(bim_path, '{}.{}.bim'.format(bim_stem, i)), sep='\t',
                       header=None, usecols=[1, 3], names=['rsid', 'snp'])
    temp.sort_values(by='snp', inplace=True)
    bim['chr{}'.format(i)] = temp


res = []  # per decile
peaks = pd.read_csv(bedfile, sep='\t', usecols=[0, 1, 2], names=['chr', 'start', 'end'])
peaks = peaks[(peaks['chr'] != 'chrX') | (peaks['chr'] != 'chrY')]

count = 0
for chr_name, val in bim.items():
    peak = peaks[peaks['chr'] == chr_name].copy()

    if len(peak) == 0:
        continue
    peak['overlap'] = 0
    peak.sort_values(by='start', inplace=True)

    i, j = 0, 0
    while i < len(peak) and j < len(val):
        start, end = peak.iloc[i, 1], peak.iloc[i, 2]
        snp = val.iloc[j, 1]

        if start <= snp <= end:
            j += 1
            peak.iloc[i, 3] += 1
        elif snp < start:
            j += 1
        else:
            i += 1
    res.append(peak[peak['overlap'] != 0])
    count += 1
    print('{}, {}/{}'.format(chr_name, count, len(bim)))
res = pd.concat(res, ignore_index=True)
basename = os.path.basename(bedfile)
res.to_csv('overlap_'+basename, sep='\t', header=False, index=False)


