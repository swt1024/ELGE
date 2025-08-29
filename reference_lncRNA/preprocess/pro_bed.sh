# Retain rows categorized as 'gene' and the required columns
awk '$4 ~ /^NONHSAG/ {sub(/\..*/, "", $4); print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}' ../human/bed/NONCODEv6_hg38.lncAndGene.bed > ../human/bed/NONCODEv6_hg38.lncRNAGene.bed
awk '$4 ~ /^NONHSAG/ {sub(/\..*/, "", $4); print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}' ../human/bed/NONCODEv5_hg38.lncAndGene.bed > ../human/bed/NONCODEv5_hg38.lncRNAGene.bed
awk '$4 ~ /^NONMMUG/ {sub(/\..*/, "", $4); print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}' ../mouse/bed/NONCODEv5_mm10.lncAndGene.bed > ../mouse/bed/NONCODEv5_mm10.lncRNAGene.bed
awk '$4 ~ /^NONMMUG/ {sub(/\..*/, "", $4); print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}' ../mouse/bed/NONCODEv6_mm10.lncAndGene.bed > ../mouse/bed/NONCODEv6_mm10.lncRNAGene.bed