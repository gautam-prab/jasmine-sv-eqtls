# given a list of genes, compare the beta and R^2 values of a set of SVs and a set of SNPs

import pandas as pd
import numpy as np
import argparse
import matplotlib.pyplot as plt
import os, glob

parser = argparse.ArgumentParser()
parser.add_argument('inputfile', help='Folder of TSVs from compare_sv_snp.py')
parser.add_argument('--idconvert', help='TSV to convert gene IDs if desired')
parser.add_argument('--r2title', help='Title of R^2 Plot')
parser.add_argument('--betatitle', help='Title of Beta Plot')
args = parser.parse_args()

idconvert = None
if args.idconvert:
    idconvert = pd.read_csv(args.idconvert, sep = '\t', index_col=0, squeeze=True)

data = pd.read_csv(args.inputfile, sep = '\t', index_col=0)

print('There are %d genes' % len(data.index))

plt.scatter(data['Max SNP Rsq'], data['Max SV Rsq'])
label = data[np.logical_and(data['Max SV Rsq'] > 0.2, data['Max SNP Rsq'] < data['Max SV Rsq'])]
for idx, row in label.iterrows():
    txt = idx
    if idconvert is not None and txt in idconvert.index:
        txt = idconvert[txt]
    text = plt.annotate(txt, (row['Max SNP Rsq'], row['Max SV Rsq']))
    text.set_fontsize(7)
plt.annotate('%.1f%%' % (100 * len(data[data['Max SNP Rsq'] < data['Max SV Rsq']]) / len(data.index)), (0.7, 0.9), xycoords='axes fraction', c='r')
plt.annotate('%.1f%%' % (100 * len(data[data['Max SNP Rsq'] >= data['Max SV Rsq']]) / len(data.index)), (0.9, 0.7), xycoords='axes fraction', c='r')

max_rsq = max(max(data['Max SNP Rsq']), max(data['Max SV Rsq']))
plt.plot(np.linspace(0,max_rsq,50), np.linspace(0,max_rsq,50), linewidth=.5, linestyle='--', c='r')
plt.ylabel('Max R^2 of SV', fontsize=16)
plt.xlabel('Max R^2 of SNP', fontsize=16)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.title(args.r2title, fontsize=18)
plt.show()

plt.scatter(abs(data[data['Max SV Beta'] > 0]['Max SNP Beta']), abs(data[data['Max SV Beta'] > 0]['Max SV Beta']), c='#1f77b4', alpha=0.5, label='+ effect')
plt.scatter(abs(data[data['Max SV Beta'] <= 0]['Max SNP Beta']), abs(data[data['Max SV Beta'] <= 0]['Max SV Beta']), c='#ff7f0e', alpha=0.5, label='- effect')
plt.legend(loc='right')
label = data[np.logical_and(abs(data['Max SV Beta']) > 1.5, abs(data['Max SNP Beta']) < abs(data['Max SV Beta']))]
for idx, row in label.iterrows():
    txt = idx
    if idconvert is not None and txt in idconvert.index:
        txt = idconvert[txt]
    text = plt.annotate(txt, (abs(row['Max SNP Beta']), abs(row['Max SV Beta'])), fontsize=9)
plt.annotate('%.1f%%' % (100 * len(data[data['Max SNP Beta'] < data['Max SV Beta']]) / len(data.index)), (0.7, 0.9), xycoords='axes fraction', c='r')
plt.annotate('%.1f%%' % (100 * len(data[data['Max SNP Beta'] >= data['Max SV Beta']]) / len(data.index)), (0.9, 0.7), xycoords='axes fraction', c='r')

max_betas = max(max(abs(data['Max SNP Beta'])), max(abs(data['Max SV Beta'])))
plt.plot(np.linspace(0,max_betas + .1,50), np.linspace(0,max_betas + .1, 50), linewidth=.5, linestyle='--', c='r')
plt.ylabel('Max Beta of SV', fontsize=16)
plt.xlabel('Max Beta of SNP', fontsize=16)
plt.xlim(0, max_betas + .1)
plt.ylim(0, max_betas + .1)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.title(args.betatitle, fontsize=18)
plt.show()
