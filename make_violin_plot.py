# make_violin_plot.py
# Make violin plots of SV/gene pairs found by process_pairs.py

import argparse
from cyvcf2 import VCF
import pickle
import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

def build_args():
    parser = argparse.ArgumentParser(description='Make violin plots of SV/gene pairs found by process_pairs.py')
    parser.add_argument('vcffile', default=None, help='VCF file from process_pairs.py')
    parser.add_argument('--exonfile', default='exons.csv', help='Generated by determine_coding_vs_noncoding.py')
    parser.add_argument('--exonic', default=False, action='store_true')
    parser.add_argument('--type_separation', default='sv', help='Separate by SV Type or Gene Type (enter \"sv\" or \"gene\")')
    return parser.parse_args()

def make_violin_plot(args, dict, title):
    fig, ax = plt.subplots()

    color = {'ALU': 'b', 'CNV': 'orange', 'DEL': 'g', 'DEL_ALU': 'r', 'INV': 'purple', 'LINE1': 'brown', 'SVA': 'pink', 'INS': 'red', 'pseudogene': 'orange', 'protein coding': 'blue', 'lncRNA': 'green'}

    labels = sorted(dict.keys())
    ax.get_xaxis().set_tick_params(direction='out')
    ax.xaxis.set_ticks_position('bottom')
    ax.set_xticks(np.arange(1, len(labels) + 1))
    ax.set_xlim(0.25, len(labels) + 0.75)
    if args.type_separation == 'sv':
        ax.set_xticklabels(labels, fontsize=12)
        ax.set_xlabel('SV Type', fontsize=16)
    elif args.type_separation == 'gene':
        ax.set_xlabel('Gene Type', fontsize=16)
        ax.set_xticklabels(labels, fontsize=12)

    ax.set_ylabel('Effect Size (beta)', fontsize=16)
    ax.set_ylim(-3, 3)
    ax.yaxis.set_tick_params(labelsize=12)
    ax.set_title(title, fontsize=16)

    names = list(sorted(dict.keys()))
    parts = ax.violinplot([dict[key] for key in sorted(dict.keys())])
    for partname in ('cbars','cmins','cmaxes'):
        # parts[partname].set_edgecolor('black')
        parts[partname].set_linewidth(1)
        parts[partname].set_edgecolor([color[names[idx]] for idx in range(len(dict.keys()))])
    for idx, pc in enumerate(parts['bodies']):
        pc.set_color(color[names[idx]])

    [bottom, top] = ax.get_ylim()
    ax.set_ylim(min(bottom, -1), max(top, 1))

    i = 1
    max_beta = 0
    for key in sorted(dict.keys()):
        samp_count = len(dict[key])
        ax.text(x=i,y=0.95,s=samp_count,fontsize=12,color='gray', horizontalalignment='center', transform=ax.get_xaxis_transform())
        i+=1

    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.show()

def main():
    args = build_args()

    if args.vcffile:
        qtls = VCF(args.vcffile)

    if args.type_separation == 'gene':
        exon_df = pd.read_csv(args.exonfile, sep=',', index_col=0)

    if args.exonic:
        ex_dict = {}
        nonex_dict = {}
        for record in qtls:
            gene = record.INFO.get('gene')
            if args.type_separation == 'sv':
                type = record.INFO.get('SVTYPE')
            elif args.type_separation == 'gene':
                type = exon_df.loc[gene.split('.')[0], 'ensembl.type_of_gene']
                if not isinstance(type, str): # no gene type found
                    continue
                elif 'IG' in type:
                    continue
                elif 'pseudogene' in type:
                    type = 'pseudogene'
                elif 'protein' in type:
                    type = 'protein coding'
                elif 'lnc' in type:
                    type = 'lncRNA'
            else:
                print('Error')
                break
            ex = record.INFO.get('exonic')

            if ex == 'y':
                if type not in ex_dict:
                    ex_dict[type] = []
                ex_dict[type].append(record.INFO.get('beta'))

            elif ex == 'n':
                if type not in nonex_dict:
                    nonex_dict[type] = []
                nonex_dict[type].append(record.INFO.get('beta'))
        make_violin_plot(args, ex_dict, 'Effect Sizes of Coding SV-Gene Pairs')
        make_violin_plot(args, nonex_dict, 'Effect Sizes of Noncoding SV-Gene Pairs')

    else:
        dict = {}
        for record in qtls:
            gene = record.INFO.get('gene')
            if args.type_separation == 'sv':
                type = record.INFO.get('SVTYPE')
            elif args.type_separation == 'gene':
                type = exon_df.loc[gene.split('.')[0], 'ensembl.type_of_gene']
                if not isinstance(type, str): # no gene type found
                    continue
                elif 'IG' in type:
                    continue
                elif 'pseudogene' in type:
                    type = 'pseudogene'
                elif 'protein' in type:
                    type = 'protein coding'
                elif 'lnc' in type:
                    type = 'lncRNA'
            else:
                print('Error')
                break

            if type not in dict:
                dict[type] = []

            dict[type].append(record.INFO.get('beta'))
        make_violin_plot(args, dict, 'Effect Sizes of SV-Gene Pairs')

if __name__ == '__main__':
    main()
