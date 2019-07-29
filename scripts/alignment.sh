#!/bin/bash

#######################
# step1 bwa alignment #
#######################

#parameters
projDir=$1
fastqFile=$2
outDir=$3
bwa=$4
samtools=$5
ref=$6
coreN=$7
mismatchN=${8:-3}
gapN=${9:-1}
summaryFile=${10:-"${outDir}/freeHiC.summary"}
shift 10
ligateSite=("$@")

name=${fastqFile##*/}
bin=${projDir}/bin
saveFiles=0 #${10:-"1"}
seqLength=25 #${11:-25} #>=25 enforced by mHiC
number='^[0-9]+$'

## refresh summary file
if [ -e "$summaryFile" ]; then
    rm -rf $summaryFile
    touch $summaryFile
fi

#Make sure the out directory exists
if [ ! -d "$outDir" ]; then
    mkdir -p $outDir
fi

# ## min length for chimeric reads is enforced to be 25bp
# if [ "$seqLength" -lt "25" ]; then
#     seqLength=25
# fi

## compile ligateSite to trim chimeric reads
g++ -std=c++11 -o $bin/ligateSite_trimming_freeHiC $bin/ligateSite_trimming_freeHiC.cpp


echo "(We assume the reads length are the same within the same fastq file.)...... "

#BWA alignment
for i in 1 2
do
    echo "Start alignment of read end $i"
    readL="$(cat ${fastqFile}_$i.fastq | head -2 | awk '{if(NR%4==2) print length($1)}')"

    ## if bwa index not exist
    if [ ! -e $ref.bwt ] || [ ! -e $ref.sa ] || [ ! -e $ref.pac ]; then
    	echo "Step 1.0 Building BWA index"
    	$bwa index $ref
    fi
    
    echo "Step1.1 - BWA alignment"
    $bwa aln -n $mismatchN -o $gapN -t $coreN $ref ${fastqFile}_$i.fastq >$outDir/${name}_$i.sai
    $bwa samse $ref $outDir/${name}_$i.sai ${fastqFile}_$i.fastq >$outDir/${name}_$i.sam
    
    if [[ "$ligateSite" =~ $number ]] && [[ "$ligateSite" -eq "0" ]]; then
	
    	echo "Choose not to rescue chimeric reads or the chimeric reads have been trimmed before alignment!"
    else
    	# step1.2 - filter get aligned sam & unmapped sam
    	echo "Step1.2 - Filter and get unmapped alignment sam file."
    	## Filter out the unmapped for future chimeric reads rescuing.
    	$samtools view -h -f 4  $outDir/${name}_$i.sam >$outDir/${name}_unmapped_$i.sam
    	mv $outDir/${name}_$i.sam $outDir/${name}_${i}_raw.sam
    
    	# step1.3 - trim and filter unmapped
    	echo "Step1.3 - Trim unmapped reads until the restriction enzyme cutting site."
    	$samtools fastq $outDir/${name}_unmapped_$i.sam >$outDir/${name}_unmapped_$i.fastq

    	for cut in ${ligateSite[@]}
    	do
	    
    	    $bin/ligateSite_trimming_freeHiC --fastq $outDir/${name}_unmapped_$i.fastq --cutsite $cut --out $outDir/${name}_unmapped_trim_$i.fastq
    	    rm $outDir/${name}_unmapped_$i.fastq
    	    mv $outDir/${name}_unmapped_trim_$i.fastq $outDir/$name\_unmapped_$i.fastq
    	done
    	awk -v minLen=$seqLength -v maxLen=$readL 'BEGIN {OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= minLen && length(seq) < maxLen) {print header, seq, qheader, qseq}}' < $outDir/${name}_unmapped_$i.fastq >$outDir/${name}_unmapped_trim_filter_$i.fastq

    	# step4 - align trimed read
    	echo "Step1.4 - Rescue chimeric reads by re-aligning trimmed unmapped reads."
    	$bwa aln -n $mismatchN -o $gapN -t $coreN $ref $outDir/${name}_unmapped_trim_filter_$i.fastq >$outDir/${name}_unmapped_trim_filter_$i.sai
    	$bwa samse $ref  $outDir/${name}_unmapped_trim_filter_$i.sai $outDir/${name}_unmapped_trim_filter_$i.fastq >$outDir/${name}_unmapped_trim_filter_$i.sam
    	$samtools view $outDir/${name}_unmapped_trim_filter_$i.sam >$outDir/${name}_unmapped_trim_filter_noheader_$i.sam

	# step5 - merge two step alignment
	echo "step1.5 - Merge chimeric reads with mapped full-length reads."	
	awk 'NR==FNR{a[$1]=$0;next;}a[$1]{$0=a[$1]}1' $outDir/${name}_unmapped_trim_filter_noheader_$i.sam  $outDir/${name}_${i}_raw.sam  >$outDir/${name}_$i.sam
    fi
done

if [[ "$ligateSite" =~ $number ]] && [[ "$ligateSite" -eq "0" ]]; then 
    rawAlign1=$(awk '$6 !="*" && $1!="@SQ" && $1!="@PG" {print $1}' $outDir/${name}_1.sam | wc -l)
    rawAlign2=$(awk '$6 !="*" && $1!="@SQ" && $1!="@PG" {print $1}' $outDir/${name}_2.sam | wc -l)
    
    echo -e "Aligned reads total 1: "$rawAlign1 >>$summaryFile
    echo -e "Aligned reads total 2: "$rawAlign2 >>$summaryFile
else
    
    rawAlign1=$(awk '$6 !="*" && $1!="@SQ" && $1!="@PG" {print $1}' $outDir/${name}_1_raw.sam | wc -l)
    rawAlign2=$(awk '$6 !="*" && $1!="@SQ" && $1!="@PG" {print $1}' $outDir/${name}_2_raw.sam | wc -l)
    
    echo -e "First stage alignment aligned reads total 1: "$rawAlign1 >>$summaryFile
    echo -e "First stage alignment aligned reads total 2: "$rawAlign2 >>$summaryFile
    
    chimericAlign1=$(awk '$6 != "*" && $1!="@SQ" && $1!="@PG" {print $1}' $outDir/${name}_unmapped_trim_filter_noheader_1.sam | wc -l)
    chimericAlign2=$(awk '$6 != "*" && $1!="@SQ" && $1!="@PG" {print $1}' $outDir/${name}_unmapped_trim_filter_noheader_2.sam | wc -l)
    
    echo -e "Second stage alignment aligned reads total 1: "$chimericAlign1 >>$summaryFile
    echo -e "Second stage alignment aligned reads total 2: "$chimericAlign2 >>$summaryFile

    echo -e $((rawAlign1 + rawAlign2 + chimericAlign1 + chimericAlign2)) >$summaryFile.alignedN
    echo -e $((chimericAlign1 + chimericAlign2)) >$summaryFile.chimericN
fi

                                                                                              

#Remove redundant files
if [ "$saveFiles" -eq "0" ]; then
    rm -rf $outDir/*sai
    rm -rf $outDir/*unmapped*
    rm -rf $outDir/*raw*
fi
