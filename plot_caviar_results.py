import argparse
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('results', help='TSV of caviar results')
parser.add_argument('--title', default='CAVIAR')
args = parser.parse_args()

data = pd.read_csv(args.results, sep='\t')
print(data)

plt.scatter(data['SNP Caviar'], data['SV Caviar'], alpha=0.8)
print('\n'.join(list(data[np.logical_and(data['SV Caviar'] > 0.8, data['SV Caviar'] > data['SNP Caviar'])]['Gene'])))
plt.annotate('%.1f%%' % (100 * len(data[data['SNP Caviar'] < data['SV Caviar']]) / len(data)), (0.7, 0.9), xycoords='axes fraction', c='r')
plt.annotate('%.1f%%' % (100 * len(data[data['SNP Caviar'] >= data['SV Caviar']]) / len(data)), (0.9, 0.7), xycoords='axes fraction', c='r')

max_val = max(max(data['SNP Caviar']), max(data['SV Caviar']))
plt.plot(np.linspace(0,max_val,1000), np.linspace(0,max_val,1000), linewidth=.5, linestyle='--', c='r')
plt.ylim([0,1])
plt.xlim([0,1])
plt.ylabel('Max CAVIAR Posterior of SV', fontsize=16)
plt.xlabel('Max CAVIAR Posterior of SNP', fontsize=16)
plt.title(args.title, fontsize=18)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.show()
