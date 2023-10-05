gtfFile="/projects/b1042/FouratiLab/genome/GRCh38/Annotation/genes.gtf"
bedFile="/projects/b1042/FouratiLab/genome/GRCh38/Annotation/genes.bed"
grep -P "\ttranscript\t" $gtfFile | \
    cut -f1,4,5,7,9 | \
    sed 's/[[:space:]]/\t/g' | \
    sed 's/[;|"]//g' | \
    awk -F $'\t' 'BEGIN { OFS=FS } { print $1,$2-1,$3,$10,".",$4,$14,$16,$18 }' | \
    sort -k1,1 -k2,2n > $bedFile
