import os
import glob
import pandas as pd

#read vdjdb data
df_vdj = pd.read_csv('../vdjdb_2024-05-26.tsv', sep='\t')
df_vdj = df_vdj[df_vdj['MHC class'] == 'MHCI'].reset_index(drop=True)
df_vdj_b = df_vdj[df_vdj['Gene'] == 'TRB']
#modify format
#CDR3B TRBV TRBJ CDR3A TRAV TRAJ epitope MHC antigen 
data = []
for ndx in df_vdj_b.index:
    complex_id = df_vdj_b.at[ndx, 'complex.id']
    if complex_id == 0:
        data.append([df_vdj_b.at[ndx, 'CDR3'], df_vdj_b.at[ndx, 'V'], df_vdj_b.at[ndx, 'J'], '', '', '', df_vdj_b.at[ndx, 'Epitope'], df_vdj_b.at[ndx, 'MHC A'], df_vdj_b.at[ndx, 'Epitope gene'], df_vdj_b.at[ndx, 'Score']])
    else:
        tcra = df_vdj[(df_vdj['complex.id'] == complex_id) & (df_vdj['Gene'] == 'TRA')]
        tcrb = df_vdj_b.iloc[ndx]
        data.append([tcrb['CDR3'], tcrb['V'], tcrb['J'], tcra['CDR3'].values[0], tcra['V'].values[0], tcra['J'].values[0], tcrb['Epitope'], tcrb['MHC A'], tcrb['Epitope gene'], tcrb['Score']])
df_pos = pd.DataFrame(data, columns=['CDR3B', 'TRBV', 'TRBJ', 'CDR3A', 'TRAV', 'TRAJ', 'epitope', 'MHC antigen', 'Epitope gene', 'vdj-Score'])
df_pos['label'] = 1
#df_pos = df_pos[df_pos['vdj-Score'] > 0].copy()
#read healthy control data
df_sc = pd.read_csv('../tcr_seq_noGD_noMAIT.csv')
df_healthy = df_sc[df_sc['PatientID'].str.contains('HC')]
df_neg = df_healthy[['TCRB_cdr3aa', 'TCRB_vgene', 'TCRB_jgene', 'TCRA_cdr3aa', 'TCRA_vgene', 'TCRA_jgene']].copy()
df_neg.drop_duplicates(subset=['TCRB_cdr3aa'], inplace=True,)
df_neg.columns = ['CDR3B', 'TRBV', 'TRBJ', 'CDR3A', 'TRAV', 'TRAJ']
df_neg = df_neg[~df_neg['CDR3B'].isin(df_pos['CDR3B'])].reset_index(drop=True)
#shuffle and sample negative data
df_epi = pd.DataFrame(df_pos['epitope'].unique(), columns=['epitope'])
df_mhc = pd.DataFrame(df_pos['MHC antigen'].unique(), columns=['MHC antigen'])
df_mhc = df_mhc[df_mhc['MHC antigen'].str.len() <= 11].reset_index(drop=True)
df_neg['epitope'] = df_epi.sample(n=df_neg.shape[0], random_state=1, replace=True).reset_index(drop=True)
df_neg['MHC antigen'] = df_mhc.sample(n=df_neg.shape[0], random_state=1, replace=True).reset_index(drop=True)
df_neg['Epitope gene'] = 'nan'
df_neg['label'] = 0

print(f'# of entries in vdjdb: {df_pos.shape[0]}')
print(f'# of entries in healthy controls: {df_neg.shape[0]}')
#remove TCRs targeting epitopes not in epitope list
#remove TCRs appeared in epitope known tcrs
df_gliph_vdj = pd.read_csv('../tcr_gliph_known-sars-2.tsv', sep='\t')
epitopes_old = list(set(df_gliph_vdj['epitope'].values.tolist()))
df_pos = df_pos[df_pos['epitope'].isin(epitopes_old)].reset_index(drop=True)
print(f'remove cdr3bs targeting epitopes not in epitope list')
print(f'# of pos. entries after filtering: {df_pos.shape[0]}')
df_pos = df_pos[~df_pos['CDR3B'].isin(df_gliph_vdj['CDR3b'])].reset_index(drop=True)
df_neg = df_neg[~df_neg['CDR3B'].isin(df_gliph_vdj['CDR3b'])].reset_index(drop=True)
print(f'remove cdr3bs appeared in epitope known tcrs')
print(f'# of pos. entries after filtering: {df_pos.shape[0]}')
##remove TCRs appeared in training datasets 
##TABR-BERT
df_train_tabr = pd.read_csv('../train_data/train_tabr.csv')
df_pos = df_pos[~df_pos['CDR3B'].isin(df_train_tabr['cdr3'])].reset_index(drop=True)
df_neg = df_neg[~df_neg['CDR3B'].isin(df_train_tabr['cdr3'])].reset_index(drop=True)
print(f'Remove cdr3bs appeared in TABR-BERT training data')
print(f'# of pos. entries after filtering: {df_pos.shape[0]}')
##pMTnet
df_train_pMTnet = pd.read_csv('../train_data/train_pmtnet.csv')
df_pos = df_pos[~df_pos['CDR3B'].isin(df_train_pMTnet['CDR3'])].reset_index(drop=True)
df_neg = df_neg[~df_neg['CDR3B'].isin(df_train_pMTnet['CDR3'])].reset_index(drop=True)
print(f'Remove cdr3bs appeared in pMTnet training data')
print(f'# of pos. entries after filtering: {df_pos.shape[0]}')
##PanPep
df_train_panpep = pd.read_csv('../train_data/train_panpep.csv')
df_pos = df_pos[~df_pos['CDR3B'].isin(df_train_panpep['binding_TCR'])].reset_index(drop=True)
df_neg = df_neg[~df_neg['CDR3B'].isin(df_train_panpep['binding_TCR'])].reset_index(drop=True)
print(f'Remove cdr3bs appeared in PanPep training data')
print(f'# of pos. entries after filtering: {df_pos.shape[0]}')
##TEIM
df_train_teim = pd.read_csv('../train_data/train_teim.csv')
df_pos = df_pos[~df_pos['CDR3B'].isin(df_train_teim['cdr3'])].reset_index(drop=True)
df_neg = df_neg[~df_neg['CDR3B'].isin(df_train_teim['cdr3'])].reset_index(drop=True)
print(f'Remove cdr3bs appeared in TEIM training data')
print(f'# of pos. entries after filtering: {df_pos.shape[0]}')
##DLpTCR
df_train_dlptcr = pd.read_csv('../train_data/train_dlptcr.csv')
df_pos = df_pos[~df_pos['CDR3B'].isin(df_train_dlptcr['CDR3'])].reset_index(drop=True)
df_neg = df_neg[~df_neg['CDR3B'].isin(df_train_dlptcr['CDR3'])].reset_index(drop=True)
print(f'Remove cdr3bs appeared in DLpTCR training data')
print(f'# of pos. entries after filtering: {df_pos.shape[0]}')
#number of entries in healthy controls
print(f'# of neg. entries in healthy controls: {df_neg.shape[0]}')
#remove entries with cdr3b length > 20 and epitope length < 8 or > 12
df_pos = df_pos[(df_pos['CDR3B'].str.len() <= 20) & (df_pos['epitope'].str.len() <= 12) & (df_pos['epitope'].str.len() >= 8)].reset_index(drop=True)
df_neg = df_neg[(df_neg['CDR3B'].str.len() <= 20) & (df_neg['epitope'].str.len() <= 12) & (df_neg['epitope'].str.len() >= 8)].reset_index(drop=True)
df_pos.to_csv('tcr_pmhc_pos.csv', index=False)
df_neg.to_csv('tcr_pmhc_neg.csv', index=False)
