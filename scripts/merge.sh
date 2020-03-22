#!/bin/bash
repName=$1
outPath=$2

## ==============================================================
## 1. Merge the restriction fragment interaciton frequency files
## ==============================================================
cat $outPath/rawDataTraining/s2_validPairs/${repName}*fragFreq | sort -k2 -k3 | awk '{print $2, $3, $1}' | awk '{a[$1" "$2]+=$3}END{for (i in a) print i,a[i]}' | sort -k1 -k2 | awk '{print $3"\t"$1"\t"$2}' >$repName.validPairs.fragFreq


## =====================================
## 2. Merge the training summary files
## =====================================
## collect all the summary files and save within the same folder
mkdir -p summaryFolder
mv ${repName}*sumFiles summaryFolder

## merge summary statistics files.
for file in summaryFolder/${repName}*sumFiles
do
    cat $file/*training.summary.alignedN >>$repName.training.summary.alignedN.raw
    cat $file/*training.summary.allMatch >>$repName.training.summary.allMatch.raw
    echo $'\r' >>$repName.training.summary.allMatch.raw
    cat $file/*training.summary.chimericN >>$repName.training.summary.chimericN.raw
    cat $file/*training.summary.revStrandN >>$repName.training.summary.revStrandN.raw
    echo $'\r' >>$repName.training.summary.revStrandN.raw
    awk -F: '{print $1"\t"$2}' $file/*training.summary.deletion >>$repName.training.summary.deletion.raw
    awk -F: '{print $1"\t"$2}' $file/*training.summary.insertion >>$repName.training.summary.insertion.raw
    awk -F: '{print $1"\t"$2}' $file/*training.summary.mismatch >>$repName.training.summary.mismatch.raw

    cat $file/*training.summary.baseQS >>$repName.training.summary.baseQS.raw

done
awk '{s+=$1} END {print s}' $repName.training.summary.alignedN.raw >$repName.training.summary.alignedN
awk '{s+=$1} END {print s}' $repName.training.summary.allMatch.raw >$repName.training.summary.allMatch
awk '{s+=$1} END {print s}' $repName.training.summary.chimericN.raw >$repName.training.summary.chimericN
awk '{s+=$1} END {print s}' $repName.training.summary.revStrandN.raw >$repName.training.summary.revStrandN
sort -k1 $repName.training.summary.deletion.raw | awk '{a[$1]+=$2}END{for (i in a) print i,a[i]}' | sort -k1 | awk '{print $1":"$2}' >$repName.training.summary.deletion
sort -k1 $repName.training.summary.insertion.raw | awk '{a[$1]+=$2}END{for (i in a) print i,a[i]}' | sort -k1 | awk '{print $1":"$2}' >$repName.training.summary.insertion
sort -k1 $repName.training.summary.mismatch.raw | awk '{a[$1]+=$2}END{for (i in a) print i,a[i]}' | sort -k1 | awk '{print $1":"$2}' >$repName.training.summary.mismatch

## Merge the sequence base quality score frequency list
python processBaseQS.py $(pwd)

## the raw folder is for backup
mkdir raw
mv *raw raw/

## ${repName}_training_sumFiles contains the merged summary files and can be used for the next simulation module.
mkdir ${repName}_training_sumFiles
mv $repName.training.summary.* ${repName}_training_sumFiles
