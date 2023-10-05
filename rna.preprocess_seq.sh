#!/bin/bash
# @author Slim Fourati (slim.fourati@northwestern.edu)
# @version 1.0

# load modules
module load samtools/1.14
module load STAR/2.7.9a
module load htseq/2.0.2
module load gcc/8.4.0
module load python/3.6.10

# read arguments
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

# set global variables for the script
bin="/home/iew5629/.local/bin"
seqDependencies="/projects/b1042/FouratiLab/genome/$genome"
genomeFasta="$seqDependencies/Sequence/genome.fa"
gtfFile="$seqDependencies/Annotation/genes.gtf"
bedFile="$seqDependencies/Annotation/genes.bed"
maxProc=8
stranded="no"

# setting project/sample directories
dirDiagnostic=$(echo $dirData | sed -r 's|/[^/]+$|/diagnostic_plot|g')


# 1. Identify sample of interest
fastq=$(find $dirData -maxdepth 1 -name '*_1.fq.gz')
fastq=( $(echo $fastq | \
          tr ' ' '\n' | \
          sort | \
          uniq) )
fastq=${fastq[$(($SLURM_ARRAY_TASK_ID - 1))]}
sample=$(echo $fastq | sed -r "s|_1.fq.gz||g")

# 2. Determine mate length
flag=true
if $flag
then
    currentDate=$(date +"%Y-%m-%d %X")
    echo -ne "$currentDate: determining mate length..."
    mateLength=$(zcat $fastq | \
        head -n 4000 | \
        awk 'NR%2==0 {print length($1)}' | \
        sort -rn | \
        head -n 1)
    # echo $mateLength
    echo "done"
fi
genomeDir="$seqDependencies/ggOverhang$(($mateLength -1))"

# 4. Read Stats
# Counting total reads and saving totalReads.txt in "diagnostic_plot" folder
flag=true
if $flag
then
    if [ ! -d $dirDiagnostic ]
    then
	mkdir -p $dirDiagnostic
    fi
    currentDate=$(date +"%Y-%m-%d %X")
    echo -ne "$currentDate: calculating total reads..."
    count=$(zcat $fastq | wc -l | awk '{print $1/4}')
    if [ $? != 0 ]
    then
        echo -ne "error\n  error counting total reads $sample\n"
        exit
    fi
    sampleid=$(echo $sample | sed -r 's/.+\///g')
    printf $sampleid'\t'$count'\t'"TotalReads"'\n' >> \
	$dirDiagnostic/ReadStats.txt
    echo "done"
fi

# 7. Alignment with 'STAR'
# Alignment of the reads to the reference genome.
# The genome is loaded into memory for alignment.
# After the end of alignment, the genome is removed from shared memory.
# Lauch STAR with options:
#   genomeDir         :  <path/to/dir/where genome has been generated>
#   genomeLoad        : <mode of shared memory usage for the genome files>
#   readFilesIn       : <mate_1.fq.gz> <mate_2.fq.gz>
#   runThreadN        : <numberOfParallelProcess>
#   outFileNamePrefix : <prefixOutputFiles>
flag=true
if $flag
then
    currentDate=$(date +"%Y-%m-%d %X")
    echo -ne "$currentDate: alignment of reads..."
    if $pairEnd
    then
	STAR --genomeDir $genomeDir \
	    --genomeLoad LoadAndRemove \
	    --readFilesIn ${sample}_1.fq.gz \
	    ${sample}_2.fq.gz \
	    --readFilesCommand zcat \
	    --runThreadN $maxProc \
 	    --outSAMtype BAM Unsorted \
	    --outFileNamePrefix ${sample}_star \
	    --outReadsUnmapped Fastx &>/dev/null
        if [ $? != 0 ]
        then
	    echo -ne "error\n  unable to aligned read in directory $sample"
	    exit 1
        fi
        # delete raw FASTQ
	rm ${sample}_starLog.out
	rm ${sample}_starLog.progress.out
	if ! $isoform
	then
	    rm ${sample}_starSJ.out.tab
	fi
	rm ${sample}_starUnmapped.out.mate1
	rm ${sample}_starUnmapped.out.mate2
	#rmdir ${sample}_star_STARtmp
    else
	STAR --genomeDir $genomeDir \
	    --genomeLoad LoadAndRemove \
	    --readFilesIn ${sample}_1.fq.gz \
	    --readFilesCommand zcat \
	    --runThreadN $maxProc \
       	    --outSAMtype BAM Unsorted \
	    --outFileNamePrefix ${sample}_star \
	    --outReadsUnmapped Fastx &>/dev/null                    
        if [ $? != 0 ]
        then
	    echo -ne "error\n  unable to aligned read in directory $sample"
	    exit
        fi
        # delete raw FASTQ
	rm ${sample}_starLog.out
	rm ${sample}_starLog.progress.out
        rm ${sample}_starUnmapped.out.mate1
	rmdir ${sample}_star_STARtmp
    fi
    echo "done"
fi


# 9. Trimmed and Mapped read stats
# Counting reads surviving after trimming and reads aligned and writing to
# 'Readstats.txt' in "diagnostic plot" folder.
flag=true
if $flag
then
     currentDate=$(date +"%Y-%m-%d %X")
     surviving=$(sed -n '6p' ${sample}_starLog.final.out | \
	 sed -r 's/.*\|\t(.*)/\1/')
     Uniqmapped=$(sed -n '9p' ${sample}_starLog.final.out | \
	 sed -r 's/.*\|\t(.*)/\1/')
     multiMap=$(sed -n '24p' ${sample}_starLog.final.out | \
         sed -r 's/.*\|\t(.*)/\1/')
     # delete unused files
     # rm ${sample}_starLog.final.out
	 
     sampleid=$(echo $sample | sed -r 's/.+\///g')
     printf $sampleid'\t'$surviving'\t'"Surviving"'\n' >> \
         $dirDiagnostic/ReadStats.txt
     printf $sampleid'\t'$Uniqmapped'\t'"UniqMapped"'\n' >> \
	 $dirDiagnostic/ReadStats.txt
     printf $sampleid'\t'$multiMap'\t'"Multimapped"'\n' >> \
	 $dirDiagnostic/ReadStats.txt
     echo "done"
fi


# 10. Order BAM file by position
flag=true
if $flag
then
    currentDate=$(date +"%Y-%m-%d %X")
    echo -ne "$currentDate: sorting bam file by name..."
    samtools sort \
	-@ $maxProc \
	-o $sample.sorted.bam \
	${sample}_starAligned.out.bam &>/dev/null
    # delete unsorted bam
    rm ${sample}_starAligned.out.bam
    echo "done"
fi

# 5. Index BAM files
flag=true
if $flag
then
    currentDate=$(date +"%Y-%m-%d %X")
    echo -ne "$currentDate: indexing bam files..."
    samtools index ${sample}.sorted.bam &>/dev/null
    echo "done"
fi

# 11. RSeQC to infer if dataset is strand-specific
flag=true
if $flag
then
    currentDate=$(date +"%Y-%m-%d %X")
    echo -ne "$currentDate: infering strand-specific experiment..."
    # pip3 install RSeQC --user
    $bin/infer_experiment.py \
	-r $bedFile \
	-i ${sample}.sorted.bam &> $sample.rseqc.infer.out
    # extract from infer_experiment output percentage of reads mapping to
    # reference strand (single-end: "++,--" or pair-end "1++,1--,2+-,2-+")
    percentRefStrand=$(grep -E '1\+\+|\"\+\+' $sample.rseqc.infer.out | \
        sed -r 's/.+: //g')
    # extract from infer_experiment output percentage of reads mapping to
    # reverse strand (single-end: "+-,-+" or pair-end "1+-,1-+,2++,2--")
    percentReverseStrand=$(grep -E '1\+\-|\"\+\-' $sample.rseqc.infer.out | \
        sed -r 's/.+: //g')
    # if one of the percentage above 80% it will be considered stranded 
    if (( $(echo "$percentRefStrand" | awk '{print ($1 >= 0.75)}') ))
    then
        stranded="yes"
    fi
    if (( $(echo "$percentReverseStrand" | awk '{print ($1 >= 0.75)}') ))
    then
        stranded="reverse"
    fi
    echo "done"
fi


# 12. Estimating transcript abundance with "HTSeq"
# Given a file with aligned reads and the list of genomic features, the number
# of reads mapped to each feature is counted, which accounts for the transcript
# abundance estimation. For each position 'i' in the read, a set S(i) is the
# set of all features overlapping position 'i'. The set 'S' which is the
# intersection of all non empty sets 'S(i)' is considered and the counts are
# generated.
# Launch HTSeq with options:
#   m        : <modeOfCounting>
#   stranded : <data is from strand-specific assay or no>
#   i        : <GTF attribute to be used as feature to count>
flag=true
if $flag
then
    currentDate=$(date +"%Y-%m-%d %X")
    echo -ne "$currentDate: estimating transcript abundance..."
    htseq-count \
        --mode=union \
        --stranded=$stranded \
        --idattr=gene_id \
	--format=bam \
        --quiet \
        $sample.sorted.bam \
        $gtfFile \
        > ${sample}_counts_gene
    if $isoform
    then
	htseq-count \
            --mode=union \
            --stranded=$stranded \
            --idattr=transcript_id \
	    --format=bam \
            --quiet \
            $sample.sorted.bam \
            $gtfFile \
            > ${sample}_counts_transcript
	htseq-count \
            --mode=union \
            --stranded=$stranded \
            --idattr=exon_id \
	    --format=bam \
            --quiet \
            $sample.sorted.bam \
            $gtfFile \
            > ${sample}_counts_exon
    fi
    # rm $sample
    echo "done"
fi


# 13. Modifying counts output
# Removing the last 5 lines from the count table
flag=true
if $flag
then
    currentDate=$(date +"%Y-%m-%d %X")
    echo -ne "$currentDate: modifying counts output..."
    tail -n 5 ${sample}_counts_gene > ${sample}.align.stat
    head -n -5 ${sample}_counts_gene > \
	${sample}_genecounts
    # deleting unuse files
    #rm ${sample}_counts_gene
    if $isoform
    then
	head -n -5 ${sample}_counts_transcript > \
	    ${sample}_transcriptcounts
	head -n -5 ${sample}_counts_exon > \
	    ${sample}_exoncounts
	#rm ${sample}_counts_transcript
	#rm ${sample}_counts_exon
    fi
    echo "done"
fi
