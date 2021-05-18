# make_manhattan_plot_ld.py
# make a manhattan plot with LD

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib as mpl
import argparse

def build_args():
    parser = argparse.ArgumentParser(description='Plot a Manhattan plot')
    parser.add_argument('inputfile', type=str, help='CSV to make Manhattan plot of')
    parser.add_argument('--column', default='p', help='Name of column to make plot of')
    parser.add_argument('--namecolumn', default='SV ID', help='Name of column of names')
    parser.add_argument('--sv', type=str, nargs='*', help='Name of SV(s)')
    parser.add_argument('--rsid', type=str, nargs='*', help='Name of GWAS SNP(s)')
    parser.add_argument('--title', type=str, default='P Values in CSV', help='Title of plot')
    parser.add_argument('--poslog', action='store_true', help='Use if you want positive log instead of negative log')
    parser.add_argument('--ld', type=str, help='LD matrix from plink')
    return parser.parse_args()

args = build_args()
df = pd.read_csv(args.inputfile, sep=',')
sv = df.loc[df[args.namecolumn].isin(args.sv)]
rs = df.loc[df[args.namecolumn].isin(args.rsid)]

other = df.drop(sv.index).drop(rs.index)

rsid0 = args.rsid[0] if isinstance(args.rsid, list) else args.rsid
ld = pd.read_csv(args.ld, sep=' ', header=None)
print(ld)
rs_idx = df[df[args.namecolumn] == rsid0].index[0]
ld_rs = ld.loc[rs_idx]
print(ld_rs)
cmap = cm.get_cmap('plasma')
other_colors = cmap(np.square(ld_rs[other.index]))
sv_colors = cmap(np.square(ld_rs[sv.index]))

vmin = np.min(np.log10(df[args.column]))
vmax = np.max(np.log10(df[args.column]))

if args.poslog:
    out_res = plt.scatter(other.index, np.log10(other[args.column]), c=other_colors, cmap=cm.plasma, vmin=vmin, vmax=vmax, label='Other SNPs')
    plt.scatter(rs.index, np.log10(rs[args.column]), c='g', cmap=cm.plasma, vmin=vmin, vmax=vmax, label='Reported SNP')
    plt.scatter(sv.index, np.log10(sv[args.column]), c=sv_colors, marker='D', cmap=cm.plasma, vmin=vmin, vmax=vmax, label='Structural Variant')
    plt.ylabel('P (log10)', fontsize=16)
else:
    out_res = plt.scatter(other.index, -np.log10(other[args.column]), c=other_colors, cmap=cm.plasma, vmin=vmin, vmax=vmax, label='Other SNPs')
    plt.scatter(rs.index, -np.log10(rs[args.column]), c='g', cmap=cm.plasma, vmin=vmin, vmax=vmax, label='Reported SNP')
    plt.scatter(sv.index, -np.log10(sv[args.column]), c=sv_colors, marker='D', cmap=cm.plasma, vmin=vmin, vmax=vmax, label='Structural Variant')
    plt.ylabel('P (-log10)', fontsize=16)
plt.legend()
plt.xlabel('Variant along Chromosome', fontsize=16)
out_res.set_cmap(cm.plasma)
bar = plt.colorbar(out_res, ax=plt.gca())
bar.set_label('LD with Reported SNP (r^2)', fontsize=16)
out_res.set_clim(0, 1)
plt.title(args.title, fontsize=18)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.show()
