import pandas as pd
import matplotlib.pyplot as plt 
from sklearn import metrics

plt.rcParams['font.family'] = 'Nimbus Roman'
plt.rcParams['font.size'] = 18
fig, ax = plt.subplots(figsize=(6,6))
df_tcr_pos = pd.read_csv('../tcr_pmhc_pos.csv')
tcr_pos = df_tcr_pos['CDR3B'].values
cdr3_colnames = ['cdr3', 'CDR3', 'CDR3', 'cdr3', 'TCRB_CDR3']
score_colnames = ['rank', 'Rank', 'Score', 'binding', 'Probability (predicted as a positive sample)']
tools = ['tabr', 'pmtnet', 'panpep', 'teim', 'dlptcr']
tool_labels = ['TABR-BERT', 'pMTnet', 'PanPep', 'TEIM', 'DLpTCR']
colors = ['blue', 'green', 'red', 'purple', 'orange']
ax.plot(ax.get_xlim(), ax.get_ylim(), ls='--', c='.3')
for ndx, tool in enumerate(tools):
    df_out = pd.read_csv(f'../output/output_{tool}.csv')
    df_out['label'] = df_out[cdr3_colnames[ndx]].apply(lambda x: 1 if x in tcr_pos else 0)
    y = df_out['label'].values
    pred = df_out[score_colnames[ndx]].values
    fpr, tpr, thresholds = metrics.roc_curve(y, pred, pos_label=1)
    auc = metrics.roc_auc_score(y, pred)
    label = tool_labels[ndx] + f' (AUC={auc:.2f})'
    plt.plot(fpr, tpr, label=label, color=colors[ndx])

df_out_gliph = pd.read_csv('./gliph/tpr_fpr_gliph.tsv', sep='\t')
fpr = df_out_gliph['fpr'].values
tpr = df_out_gliph['tpr'].values
plt.scatter(fpr, tpr, label='GLIPH', color='black', marker='*', s=150)
ax.set_xlim([0, 1])
ax.set_ylim([0, 1])
ax.set_xlabel('False Positive Rate')
ax.set_ylabel('True Positive Rate')
plt.legend(loc='lower right', frameon=True, prop={'size': 12})
plt.tight_layout()
plt.savefig('roc_curve.png', dpi=600, bbox_inches='tight')
