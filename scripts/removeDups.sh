#!/bin/bash

## mHi-C
## Author: Ye Zheng
## Contact: yezheng@stat.wisc.edu
## Update: May 2018

###################################################
# step 4 remove duplicates among Valid Interactions
###################################################
validP=$1
outDir=${validP%/*}
splitByChrom=$2
summaryFile=$3
fragFreq=$4
if [ ! -d $outDir/sorttmp ]; then
    mkdir -p $outDir/sorttmp
fi

number='^[0-9]+$'
if [ "$splitByChrom" -eq 1 ];then
    ## split validPairs
    for c in ${chrList[@]} ##$(seq 1 22) X Y
    do
	if [[ $c =~ $number ]]; then
	    chrom="chr"$c
	else
	    chrom=$c
	fi
	    
	## remove PCR duplicates based on alignment chrom + position
	awk -v chrom="$chrom" '$2 == chrom {print $0}' $validP |  sort -k2,2V -k3,3n -k4 -k7,7V -k8,8n -k9 -k1 -T $outDir/sorttmp | awk -v OFS="\t" 'BEGIN{c1=0;c2=0;p1=0;p2=0;s1=0;s2=0;}(c1!=$2 || c2!=$7 || p1!=$3 || p2!=$8 || s1!=$4 || s2!=$9 ){print;c1=$2;c2=$7;p1=$3;p2=$8;s1=$4;s2=$9}' >$validP.$chrom
	
	
    done
    cat $validP.chr* >$validP.nodup
    rm -rf $validP.chr*
    
else
    sort -k2,2V -k3,3n -k4 -k7,7V -k8,8n -k9 -k1 -T $outDir/sorttmp $validP | awk -v OFS="\t" 'BEGIN{c1=0;c2=0;p1=0;p2=0;s1=0;s2=0;}(c1!=$2 || c2!=$7 || p1!=$3 || p2!=$8 || s1!=$4 || s2!=$9 ){print;c1=$2;c2=$7;p1=$3;p2=$8;s1=$4;s2=$9}' >$validP.nodup
    
fi

if [ "$fragFreq" -eq "1" ]; then
    awk '{print $5, $10}' $validP.nodup | sort -T $outDir/sorttmp | uniq -c | awk '{print $1"\t"$2"\t"$3}' >$validP.fragFreq
fi
rm -rf $outDir/sorttmp

validP_num=$(cat $validP | wc -l)
validP_nodup_num=$(cat $validP.nodup | wc -l)


echo -e "  Valid interactions count:\t"$validP_num >>$summaryFile
echo -e "  Valid interactions without duplicates count:\t"$validP_nodup_num >>$summaryFile
