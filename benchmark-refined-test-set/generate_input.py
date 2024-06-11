import random
import pandas as pd

#read data
df_pos = pd.read_csv('tcr_pmhc_pos.csv')
df_neg = pd.read_csv('tcr_pmhc_neg.csv')
num_pos = df_pos.shape[0]
df_neg_subsample = df_neg.sample(n=num_pos, random_state=1)
df = pd.concat([df_pos, df_neg_subsample], ignore_index=True)

out_dir = './input'
#generate input files for different tools
#tabr-bert
df_input = pd.DataFrame()
df_input[['cdr3', 'allele', 'peptide']] = df[['CDR3B', 'MHC antigen', 'epitope']]
df_input.to_csv(f'{out_dir}/input_tabr.csv', index=False)
#pMTnet
df_input = pd.DataFrame()
df_input[['CDR3', 'Antigen', 'HLA']] = df[['CDR3B', 'epitope', 'MHC antigen']]
df_input['HLA'] = df_input['HLA'].apply(lambda x: x.split('-')[1] if '-' in x else x)
df_input.to_csv(f'{out_dir}/input_pmtnet.csv', index=False)
#Panpep
df_input = pd.DataFrame()
df_input[['Peptide', 'CDR3']] = df[['epitope', 'CDR3B']]
df_input.to_csv(f'{out_dir}/input_panpep.csv', index=False)
#TEIM
df_input = pd.DataFrame()
df_input[['cdr3', 'epitope']] = df[['CDR3B', 'epitope']]
df_input.to_csv(f'{out_dir}/input_teim.csv', index=False)
#DLpTCR
df_input = pd.DataFrame()
df_input['TCR_CDR3_Alpha'] = 'nan'
df_input[['TCR_CDR3_Beta', 'Peptide']] = df[['CDR3B', 'epitope']]
df_input.to_excel(f'{out_dir}/input_dlptcr.xlsx', index=False)
#GLIPH
df_input = pd.DataFrame()
df_input[['CDR3b','epitope']] = df[['CDR3B', 'epitope']].copy()
df_input['source'] = 'unknown'
df_gliph = pd.read_csv('../tcr_gliph_known-sars-2.tsv', sep='\t')
df_gliph_format = df_gliph[['CDR3b', 'epitope']].copy()
df_gliph_format['source'] = 'known'
df_input = pd.concat([df_input, df_gliph_format], ignore_index=True)
df_input.to_csv(f'{out_dir}/input_gliph.tsv', sep='\t', index=False)


