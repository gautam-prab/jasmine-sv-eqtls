import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from cyvcf2 import VCF
from scipy.stats import mannwhitneyu, linregress
import argparse

def build_args():
    parser = argparse.ArgumentParser(description='Plot the genotype-phenotype regression')
    parser.add_argument('inputfile', type=str, help='VCF (with gene annotations) of variants to make boxplots of')
    parser.add_argument('--gene', type=str, help='(optional) gene to pair SVs with')
    parser.add_argument('--gene_names', type=str, help='Gene names TSV from BioMart')
    parser.add_argument('--rnafile', default='Data/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt')
    parser.add_argument('--sv', type=str, help='(optional) only use one SV ID')
    parser.add_argument('--remove_outliers', default=None, type=float, help='Remove extreme expression outliers >n standard deviations')
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
    samples = vcf_reader.samples

    names_df = pd.read_csv(args.gene_names, sep='\t', index_col=0)

    alt_dict = {'<CN0>' : 0, '<CN2>' : 2, '<CN3>' : 3, '<CN4>' : 4, '<CN5>' : 5, '<CN6>' : 6, '<CN7>' : 7, '<CN8>' : 8, '<CN9>' : 9}

    for v in vcf_reader:
        if (args.sv and v.ID != args.sv):
            continue
        if (args.gene):
            gene = args.gene
        else:
            gene = v.INFO.get('gene')
        phenotypes = rna_df.loc[gene].iloc[3:]
        uncorrected_phenotypes = phenotypes
        phenotypes = (phenotypes - np.mean(phenotypes)) / np.std(phenotypes)

        gts = v.gt_types
        phens = np.array(phenotypes[np.asarray(samples)[gts != 3]], dtype=float)
        gts = gts[gts != 3]

        # alts = v.ALT
        # if v.INFO.get('SVTYPE') == 'CNV':
        #     gt = v.genotypes[vcf_samples.index(sample)]
        # for idx, sample in enumerate(samples):
        #     gt = v.genotypes[vcf_samples.index(sample)]
        #
        #     if v.INFO.get('SVTYPE') == 'CNV' and alts[0] in alt_dict:
        #         gts[idx] = resolve_alts(gt, alts, alt_dict) # resolve copy numbers
        #     else:
        #         gts[idx] = gt[0] + gt[1]
        #     phens[idx] = phenotypes[sample]

        if args.remove_outliers is not None:
            gts = gts[phens < args.remove_outliers]
            phens = phens[phens < args.remove_outliers]
        # make a boxplot
        data = []
        copy_nums = np.unique(gts)
        cp_dict = {}
        for idx, num in enumerate(copy_nums):
            cp_dict[num] = idx + 1
            data.append(phens[gts == num])
        plt.boxplot(data, showfliers=False)

        # plot a regression line
        line_x = np.arange(1, len(cp_dict) + 1) # this is where the boxplots will be
        gts_as_x = np.array([cp_dict[gt] for gt in gts.astype(int)], dtype=float) # as x coordinates
        slope, intercept, r_value, p_value, std_err = linregress(gts_as_x.astype(float), phens)

        name = names_df.loc[gene.split('.')[0], 'Gene name']
        print(v.ID)
        print(gene)
        print('Beta is %.2f' % slope)
        print('R^2 is %.2f' % (r_value ** 2))
        line_y = line_x * slope + intercept
        plt.plot(line_x,line_y,linestyle='--')

        # add datapoints with 'jitter'
        gts_jitter = np.random.normal(np.zeros(len(gts)), 0.04)
        plt.scatter(gts_as_x + gts_jitter, phens, c=gts_as_x, cmap='rainbow', marker='.', alpha=0.4)

        # make figure
        plt.xlabel('Genotype/Copy Number', fontsize=18)
        plt.ylabel('Phenotype', fontsize=18)
        plt.xticks(line_x, labels=copy_nums, fontsize=16)
        plt.yticks(fontsize=16)
        plt.title('Effect of %s on %s' % (v.ID, name), fontsize=20)
        plt.show()

if __name__ == '__main__':
    args = build_args()
    main(args)
