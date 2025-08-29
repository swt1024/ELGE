# Create directories for human and mouse genome data
mkdir -p ./human
mkdir -p ./mouse

# Download conservation scores for human genome (hg38) from UCSC
wget -O ./human/phastCons.bw http://hgdownload.cse.ucsc.edu/goldenPath/hg38/phastCons100way/hg38.phastCons100way.bw
wget -O ./human/phyloP.bw http://hgdownload.cse.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw

# Download conservation scores for mouse genome (mm10) from UCSC
wget -O ./mouse/phastCons.bw "http://hgdownload.cse.ucsc.edu/goldenPath/mm10/phastCons60way/mm10.60way.phastCons.bw"
wget -O ./mouse/phyloP.bw "http://hgdownload.cse.ucsc.edu/goldenPath/mm10/phyloP60way/mm10.60way.phyloP60way.bw"

echo "All have been downloaded."
