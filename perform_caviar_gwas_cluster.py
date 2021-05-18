import subprocess
import argparse
import pandas as pd
import numpy as np
import re
from process_pairs import main

parser = argparse.ArgumentParser(description='Run CAVIAR given SNP file and SV file')
parser.add_argument('--snpfile', type=str, required=True, help='file with other SNPs (.vcf.gz)')
parser.add_argument('--rsfile', type=str, required=True, help='file with SNPs of interest (.vcf.gz)')
parser.add_argument('--svfile', type=str, required=True, help='file with SVs (.vcf.gz)')
parser.add_argument('--rnafile', type=str, required=True, help='file with gene RNA-seq (.txt)')
parser.add_argument('--gene', type=str, required=True, help='ENSEMBL ID to associate')
parser.add_argument('--caviarpath', type=str, default='binaries/CAVIAR', help='path to CAVIAR binary')
parser.add_argument('--outfolder', type=str, required=True, help='folder to put output in')
args = parser.parse_args()

# Merge variant files
with open(args.outfolder + 'all_variants.vcf.gz', 'w') as file:
    subprocess.run('bcftools concat -a -Oz {} {} {} | bcftools sort -Oz'.format(args.snpfile, args.svfile, args.rsfile), shell=True, stdout=file)
subprocess.run(['tabix', '-p', 'vcf', args.outfolder + 'all_variants.vcf.gz'])

# Find top 1000 variants
pairs_parser = argparse.ArgumentParser(description='Calculate p and beta values for each SV/gene pair')
pairs_parser.add_argument('--rnafile', default='Data/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt', help='TSV of RNA seq data for samples of interest')
pairs_parser.add_argument('inputfile', type=str, help='VCF file containing genotypes and paired genes')
pairs_parser.add_argument('--gene', default=None, help='Ensembl name of a single gene to associate all variants to, if inputfile does not have paired genes')
pairs_parser.add_argument('--outcsv', default=None, help='CSV out filename')
pairs_parser.add_argument('--outvcf', default=None, help='VCF out filename')
pair_args = pairs_parser.parse_args([args.outfolder + 'all_variants.vcf.gz', '--gene', args.gene, '--outcsv', args.outfolder + 'all_pairs.csv', '--rnafile', '~/scratch/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt'])
main(pair_args)
all_df = pd.read_csv(args.outfolder + 'all_pairs.csv')
all_df = all_df[all_df['p'] != 0]
sorted_df = all_df.sort_values(by='p')
sorted_df.iloc[:1000]['SV ID'].to_csv(args.outfolder + 'ids.txt', header=False, index=False)
with open(args.outfolder + 'sig_variants.vcf.gz', 'w') as file:
    subprocess.run('bcftools filter -i \'ID=@{}\' -Oz {} | bcftools norm -m - -Oz | bcftools sort -Oz'.format(args.outfolder + 'ids.txt', args.outfolder + 'all_variants.vcf.gz'), shell=True, stdout=file)

# Reprocess pairs to get Z in location-sorted order
pair_args2 = pairs_parser.parse_args([args.outfolder + 'sig_variants.vcf.gz', '--gene', args.gene, '--outcsv', args.outfolder + 'sig_pairs.csv', '--rnafile', '~/scratch/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt'])
main(pair_args2)
sig_df = pd.read_csv(args.outfolder + 'sig_pairs.csv')
pd.Series(np.array(sig_df['beta'] / sig_df['std err']), index=np.array(sig_df['SV ID'])).to_csv(args.outfolder + 'Z.txt', sep = '\t', header = False)

# Calculate LD values
subprocess.run(['gunzip', args.outfolder + 'sig_variants.vcf.gz'])
# edit sig_variants.vcf since vcftools won't like it
subprocess.run('module load vcftools', shell=True)
with open(args.outfolder + 'sig_variants.vcf', 'r') as infile, open(args.outfolder + 'sig_variants_vcftools.vcf', 'w') as outfile:
    for line in infile:
            if '##fileformat' in line:
                outfile.write('##fileformat=VCFv4.1\n')
            else:
                outfile.write(re.sub('".*?"', '""', line))
subprocess.run('ulimit -n 3000 && vcftools --vcf {} --plink --out {}'.format(args.outfolder + 'sig_variants_vcftools.vcf', args.outfolder + 'plink'), shell=True)
subprocess.run('module load plink', shell=True)
subprocess.run(['plink', '--file', args.outfolder + 'plink', '--r', '--matrix', '--out', args.outfolder + 'plink_out'])
subprocess.run(['mv', args.outfolder + 'plink_out.ld', args.outfolder + 'ld.txt'])

print('Running CAVIAR...')
subprocess.run('module load gcc', shell=True)
subprocess.run(['/home-2/gprabhu1@jhu.edu/caviar/CAVIAR-C++/CAVIAR', '-z', args.outfolder + 'Z.txt', '-l', args.outfolder + 'ld.txt', '-o', args.outfolder + 'caviar_out'])
