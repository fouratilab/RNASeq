#!/bin/bash
# @author Slim Fourati (slim.fourati@northwestern.edu)
# @version 1.0

# read input arguments
email="slim.fourati@northwestern.edu"
pairEnd=false
isoform=false
genome=GRCh38
acceptedGenome=("GRCh38" "GRCm39")

while getopts :d:e:g:pih option
do
    case "${option}" in
	h) echo "Command: bash rna.preprocess_master.sh -d {fastq/directoryfastq} ..."
	    echo "argument: d=[d]irectory with raw data (required)"
	    echo "          g=reference [g]enome"
	    echo "          p=[p]aired-end sequencing"
	    echo "          i=[i]soform transcript/exon counts"
	    echo "          e=[e]mail address"
	    echo "          h=print [h]elp"
	    exit 1;;
	d) dirFastq=$OPTARG;;
	e) email=$OPTARG;;
	g) genome=$OPTARG
	    if [[ ! "${acceptedGenome[@]}" =~ "$genome" ]]
	    then
		echo "Invalid -g argument: choose between ${acceptedGenome[@]}"
		exit 1
	    fi;;
	p) pairEnd=true;;
	i) isoform=true;;
	\?) echo "Invalid option: -$OPTARG"
	    exit 1;;
	:)
	    echo "Option -$OPTARG requires an argument."
	    exit 1;;
    esac
done

# test that directory is provided
if [ -z ${dirFastq+x} ]
then
    echo "error...option -d required."
    exit 1
fi

# test that directory contains seq files
suffix="fq.gz"
nfiles=$(find $dirFastq -name "*_1.$suffix" | wc -l)
if [ $nfiles -lt 1 ]
then
    echo "error...empty input directory"
    exit 1
fi

# lauch genome indexing
sed -ri "s|^#SBATCH --mail-user=.+$|#SBATCH --mail-user=${email}|g" \
    genomeGenerate.slurm
cmd="sbatch genomeGenerate.slurm -d $dirFastq -g $genome"

# echo $cmd
slurmid=$(eval $cmd | sed -r 's|Submitted batch job ([0-9]*)|\1|g')

# modify preprocessing slurm script
sed -ri "s|^#SBATCH --mail-user=.+$|#SBATCH --mail-user=${email}|g" \
    rna.preprocess.slurm
sed -ri "s|^#SBATCH --array=1-.+$|#SBATCH --array=1-${nfiles}|g" \
    rna.preprocess.slurm
sed -ri "s|^#SBATCH --depend=afterok:.+$|#SBATCH --depend=afterok:${slurmid}|g" \
    rna.preprocess.slurm

# lauch preprocessing slurm script
cmd="sbatch rna.preprocess.slurm -d $dirFastq -g $genome"

if $pairEnd
then
    cmd="$cmd -p"
fi

if $isoform
then
    cmd="$cmd -i"
fi

# echo $cmd
slurmid=$(eval $cmd | sed -r 's|Submitted batch job ([0-9]*)|\1|g')
