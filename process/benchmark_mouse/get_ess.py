import pandas as pd

lncRNA = pd.read_csv("../../data/LPI/mouse/lncRNA.csv")
ess_lnc = pd.read_csv("../../data/benchmark/mouse/mouse_esslnc.csv")

ess_lnc['Noncode_id'] = ess_lnc['Noncode_id'].str.split('.').str[0]

def is_essential(x, symbol):
    if x != '-' and any(ess_lnc[col].isin([x]).any() for col in ['Noncode_id', 'Lncbook_id', 'ensembl_id']):
        return 1
    if symbol != '-' and any(ess_lnc[col].isin([symbol]).any() for col in ['gene_name', 'lib_id']):
        return 1
    return 0

lncRNA['essential'] = lncRNA.apply(lambda row: is_essential(row['gene_id'], row['symbol']), axis=1)

ess_lnc_lpi = lncRNA[lncRNA['essential'] == 1]
ess_lnc_lpi = ess_lnc_lpi[['lncRNA_ID']]

ess_lnc_lpi.to_csv('../../data/benchmark/mouse/ess_lpi.csv', index=False)
