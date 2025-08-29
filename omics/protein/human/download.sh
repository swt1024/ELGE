wget -O pLoF.txt.bgz https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz
 
bgzip -d -c pLoF.txt.bgz > pLoF.txt
