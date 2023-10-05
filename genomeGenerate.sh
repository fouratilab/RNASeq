#!/bin/bash
# @author Slim Fourati (slim.fourati@northwestern.edu)
# @version 1.0

# load apps
module load STAR/2.7.9a

# read input arguments
compress=false
while getopts d:g: option
do
    case "${option}" in
	d) dirData=$OPTARG;;
	g) genome=$OPTARG;;
    esac
done

# set global variables for the script
genomeDir="/projects/b1042/FouratiLab/genome/${genome}"
genomeFasta="$genomeDir/Sequence/genome.fa"
gtfFile="$genomeDir/Annotation/genes.gtf"
maxProc=8

# 1. determine mate length
currentDate=$(date +"%Y-%m-%d %X")
echo -ne "$currentDate: determining mate length..."
if [[ ! -z "$dirData" ]]
then
    file=$(find $dirData -name "*.fq.gz" | head -n 1)
    mateLength=$(zcat $file | \
	head -n 4000 | \
	awk 'NR%2==0 {print length($1)}' | \
	sort -rn | \
	head -n 1)
fi
echo -ne "done\n"

genomeDir="$genomeDir/ggOverhang$(($mateLength -1))"

# 2. Generating genomes for alignment with 'STAR'
currentData=$(date +"%Y-%m-%d %X")
echo -ne "$currentData: generating genomes for alignment..."
if [ -e "$genomeDir/SA" ]
then
    rm -r $genomeDir
fi
mkdir -p $genomeDir
STAR \
    --runMode genomeGenerate \
    --genomeDir $genomeDir \
    --genomeFastaFiles $genomeFasta \
    --genomeSAsparseD 2 \
    --runThreadN $maxProc \
    --sjdbGTFfile $gtfFile \
    --sjdbOverhang $(($mateLength - 1)) &>/dev/null
if [ $? != 0 ]
then
    echo -ne "error\n  unable to index genome"
    exit
fi
currentDir=$(pwd)
if [ -e "$currentDir/Log.out" ]
then
    rm $currentDir/Log.out
fi
# change permission of index directory
chmod g+w $genomeDir
echo -ne "done\n"
