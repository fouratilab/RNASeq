#!/bin/bash

# user email address
#SBATCH --mail-user=EMAIL

# mail is sent to you when the job starts and when it terminates or aborts
#SBATCH --mail-type=END,FAIL

# name of job
#SBATCH --job-name=rna.pre_%a

# standard output file
#SBATCH --output=rna.pre_%a.out

# number of nodes and processors, memory required
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32gb
#SBATCH --exclusive

# time requirements
#SBATCH --time=48:00:00

# dependencies
#SBATCH --depend=afterok:$SLURM_JOB_ID

# create array
#SBATCH --array=1-BATCH

# add partition
#SBATCH --account=b1042
#SBATCH --partition=genomics

# initialize directories
pairEnd=false
isoform=false
while getopts d:g:pi option
do
    case "$option" in
        d) dirData=$OPTARG;;
        g) genome=$OPTARG;;
	p) pairEnd=true;;
	i) isoform=true;;
    esac
done

# launch executable script
cmd="bash rna.preprocess_seq.sh -d $dirData"
cmd="$cmd -g $genome"

if $pairEnd
then
    cmd="$cmd -p"
fi

if $isoform
then
    cmd="$cmd -i"
fi

# echo $cmd
eval $cmd
