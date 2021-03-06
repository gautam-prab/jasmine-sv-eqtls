# calculate p and beta values for SV/gene pairs

import numpy as np
import pandas as pd
from cyvcf2 import VCF, Writer
from scipy.stats import mannwhitneyu, linregress
import argparse

def build_args():
    parser = argparse.ArgumentParser(description='Calculate p and beta values for each SV/gene pair')
    parser.add_argument('--rnafile', default='Data/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt', help='TSV of RNA seq data for samples of interest')
    parser.add_argument('inputfile', type=str, help='VCF file containing genotypes and paired genes')
    parser.add_argument('--gene', default=None, help='Ensembl name of a single gene to associate all variants to, if inputfile does not have paired genes')
    parser.add_argument('--outcsv', default=None, help='CSV out filename')
    parser.add_argument('--outvcf', default=None, help='VCF out filename')
    return parser.parse_args()

def resolve_alts(gt, alts, alt_dict):
    ret = 0
    if int(gt[0]) != 0: # not reference allele
        alt = alts[int(gt[0]) - 1]
        ret += alt_dict[alt]
    else: # reference allele is copy number 1
        ret += 1
    if int(gt[1]) != 0: # not reference allele
        alt = alts[int(gt[1]) - 1]
        ret += alt_dict[alt]
    else:
        ret += 1
    return ret

def main(args):
    rna_df = pd.read_csv(args.rnafile, sep='\t', index_col=0)
    rna_samples = rna_df.iloc[:,3:].keys().values

    vcf_reader = VCF(args.inputfile, gts012=True)
    samples = np.array(vcf_reader.samples)

    if args.outvcf:
        vcf_reader.add_info_to_header({'ID': 'p', 'Description': 'association p value', 'Type':'Float', 'Number': '1'})
        vcf_reader.add_info_to_header({'ID': 'beta', 'Description': 'association effect size', 'Type':'Float', 'Number': '1'})
        vcf_reader.add_info_to_header({'ID': 'rsq', 'Description': 'r squared of regression', 'Type':'Float', 'Number': '1'})
        vcf_reader.add_info_to_header({'ID': 'stderr', 'Description': 'std. err. of regression', 'Type':'Float', 'Number': '1'})
        w = Writer(args.outvcf, vcf_reader)

    alt_dict = {'<CN0>' : 0, '<CN2>' : 2, '<CN3>' : 3, '<CN4>' : 4, '<CN5>' : 5, '<CN6>' : 6, '<CN7>' : 7, '<CN8>' : 8, '<CN9>' : 9}

    out_dict = {}
    out_dict['SV ID'] = []
    out_dict['Gene ID'] = []
    out_dict['p'] = []
    out_dict['beta'] = []
    out_dict['r^2'] = []
    out_dict['std err'] = []

    for v in vcf_reader:
        if args.gene:
            genes = args.gene
        else:
            genes = v.INFO.get('gene')
        if genes is None:
            continue
        for gene in genes.split(','):
            gene = gene.split('.')[0] # remove version
            line = rna_df.loc[rna_df.index.str.contains(gene, na=False)]
            if line.empty: # gene not found
                continue
            phenotypes = line.iloc[0].iloc[3:]
            phenotypes = (phenotypes - np.mean(phenotypes)) / np.std(phenotypes)

            gts = np.zeros(len(samples))
            phens = np.zeros(len(samples))

            alts = v.ALT
            type = v.INFO.get('SVTYPE')
            if type == 'CNV' and alts[0] in alt_dict:
                copynum = True
                for idx, sample in enumerate(samples):
                    gt = v.genotypes[list(samples).index(sample)]
                    gts[idx] = resolve_alts(gt, alts, alt_dict) # resolve copy numbers
                    phens[idx] = phenotypes[sample]
            else:
                copynum = False
                gts = v.gt_types
                phens = np.array(phenotypes[samples[gts != 3]], dtype=float)
                gts = gts[gts != 3]

            if copynum:
                _,p = mannwhitneyu(phens[gts != 2], phens[gts == 2], alternative='two-sided')
            else:
                _,p = mannwhitneyu(phens[gts != 0], phens[gts == 0], alternative='two-sided')

            beta, _, r_value, _, std_err = linregress(gts, phens)

            out_dict['SV ID'].append(v.ID)
            out_dict['Gene ID'].append(gene)
            out_dict['p'].append(p)
            out_dict['beta'].append(beta)
            out_dict['r^2'].append(r_value ** 2)
            out_dict['std err'].append(std_err)
            if args.outvcf:
                v.INFO['p'] = p
                v.INFO['beta'] = beta
                v.INFO['rsq'] = (r_value ** 2)
                v.INFO['stderr'] = std_err
                w.write_record(v)

    if args.outcsv:
        pd.DataFrame.from_dict(out_dict, orient='columns').to_csv(args.outcsv)

    vcf_reader.close()
    if args.outvcf:
        w.close()

if __name__ == '__main__':
    args = build_args()
    main(args)
