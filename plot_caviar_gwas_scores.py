import argparse
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('folder', help='Directory containing directory for each gene')
parser.add_argument('--gene_conversion', help='Gene ID conversion info')
parser.add_argument('--sv_csv', help='CSV of SV associations')
parser.add_argument('--snp_csv', help='CSV of SNP associations')
parser.add_argument('--title', help='Title of plot')
args = parser.parse_args()

max_snp_causal = []
max_sv_causal = []
gene_names = []

total = 0
greater = 0
sig_total = 0
sig_greater = 0
if args.gene_conversion and args.sv_csv and args.snp_csv:
    gene_conversion = pd.read_csv(args.gene_conversion, sep = '\t', names = ['Ensembl', 'Name'], index_col = 'Name')
    sv_csv = pd.read_csv(args.sv_csv)
    snp_csv = pd.read_csv(args.snp_csv)
    print('Gene\tSV\t\t\tSNP\t\t\t')
    print('\tName\tp\tCaviar\tName\tp\tCaviar')

for entry in os.scandir(args.folder):
    if entry.is_dir():
        gene_names.append(entry.name)
        caviar = pd.read_csv(entry.path + '/caviar_out_post', sep='\t')
        sv_list = pd.read_csv(entry.path + '/sv_ids.txt', squeeze=True, names=['SV'])
        snp_list = pd.read_csv(entry.path + '/rs_ids.txt', squeeze=True, names=['SNP'])

        if not caviar.loc[caviar['SNP_ID'].isin(sv_list)].empty:
            max_sv_loc = caviar.loc[caviar['SNP_ID'].isin(sv_list), 'Causal_Post._Prob.'].idxmax()
            max_sv = caviar.loc[max_sv_loc, 'Causal_Post._Prob.']
        else:
            max_sv = 0
        if not caviar.loc[caviar['SNP_ID'].isin(snp_list)].empty:
            max_snp_loc = caviar.loc[caviar['SNP_ID'].isin(snp_list), 'Causal_Post._Prob.'].idxmax()
            max_snp = caviar.loc[max_snp_loc, 'Causal_Post._Prob.']
        else:
            max_snp = 0

        max_sv_causal.append(max_sv)
        max_snp_causal.append(max_snp)

        total += 1
        if max_sv > max_snp:
            greater += 1
        if max_sv > 0.2 or max_snp > 0.2:
            sig_total += 1
            if max_sv > max_snp:
                sig_greater += 1

        if max_sv > 0.01 and args.gene_conversion and args.sv_csv and args.snp_csv:
            gene_id = gene_conversion.loc[entry.name, 'Ensembl']
            sv_id = caviar.loc[max_sv_loc, 'SNP_ID']
            snp_id = caviar.loc[max_snp_loc, 'SNP_ID']
            sv_beta = sv_csv.loc[np.logical_and(sv_csv['Gene ID'].str.contains(gene_id), sv_csv['SV ID'] == sv_id), 'p'].iloc[0]
            snp_beta = snp_csv.loc[np.logical_and(snp_csv['Gene ID'].str.contains(gene_id), snp_csv['SV ID'] == snp_id), 'p'].iloc[0]
            sv_info = '{}\t{}\t{}'.format(sv_id, sv_beta, caviar.loc[max_sv_loc, 'Causal_Post._Prob.'])
            snp_info = '{}\t{}\t{}'.format(snp_id, snp_beta, caviar.loc[max_snp_loc, 'Causal_Post._Prob.'])
            print('{}\t{}\t{}'.format(entry.name, sv_info, snp_info))

print('{}/{} above the line'.format(greater, total))
print('{}/{} above 0.2 and above the line'.format(sig_greater, sig_total))

plt.scatter(max_snp_causal, max_sv_causal, alpha=0.8)
for i, txt in enumerate(gene_names):
    if max_sv_causal[i] > 0.2 and max_snp_causal[i] < max_sv_causal[i]:
        if max_sv_causal[i] > 0.8:
            print(txt)
        # text = plt.annotate(txt, (max_snp_causal[i], max_sv_causal[i]))
        # text.set_fontsize(7)
max_sv_arr = np.array(max_sv_causal)
max_snp_arr = np.array(max_snp_causal)
plt.annotate('%.1f%%' % (100 * len(max_sv_arr[max_snp_arr < max_sv_arr]) / len(max_sv_arr)), (0.7, 0.9), xycoords='axes fraction', c='r')
plt.annotate('%.1f%%' % (100 * len(max_sv_arr[max_snp_arr >= max_sv_arr]) / len(max_sv_arr)), (0.9, 0.7), xycoords='axes fraction', c='r')

max_val = max(max(max_snp_causal), max(max_sv_causal))
plt.plot(np.linspace(0,max_val,1000), np.linspace(0,max_val,1000), linewidth=.5, linestyle='--', c='r')
plt.ylim([0,1])
plt.xlim([0,1])
plt.ylabel('Max CAVIAR Posterior of SV', fontsize=16)
plt.xlabel('Max CAVIAR Posterior of SNP', fontsize=16)
plt.title(args.title, fontsize=18)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.show()
