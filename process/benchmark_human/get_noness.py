import pandas as pd 

ess_lnc_lpi = pd.read_csv("../../data/benchmark/human/ess_lpi.csv")
lnc_GIC = pd.read_csv("sorted_GIC_score.csv")

annotated_lnc = pd.read_csv("../../Annotate/human/valid_heart_annotation.csv")

unlabel_lnc = lnc_GIC[~lnc_GIC['lncRNA_ID'].isin(ess_lnc_lpi['lncRNA_ID'])]
unlabel_lnc = unlabel_lnc[unlabel_lnc['lncRNA_ID'].isin(annotated_lnc['lncRNA_ID'])]

ess_counts = ess_lnc_lpi.shape[0]

noness_lnc_lpi = unlabel_lnc[:ess_counts]

noness_lnc_lpi = noness_lnc_lpi[['lncRNA_ID']]

noness_lnc_lpi.to_csv("../../data/benchmark/human/noness_lpi.csv", index=False)