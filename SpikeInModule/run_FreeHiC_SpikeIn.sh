#!/bin/bash
paraFile=$1
source $paraFile

## ==========================
## print the parameter values
## ==========================
echo "projDir=${projDir}"
echo "spikeInType=${spikeInType}"
echo "dataPath=${dataPath}"
echo "resolution=${resolution}"
echo "thres=${thres}"
echo "perc=${perc}"
echo "nei=${nei}"
echo "outPath=${outPath}"

gunzip $dataPath/GM12878_rep4_chr1.fragFreq.gz
gunzip $dataPath/GM12878_rep6_chr1.fragFreq.gz

python3 $projDir/SpikeInModule/FreeHiC_SpikeIn.py -s $spikeInType -c1 $dataPath/GM12878_rep2_chr1.binPairs $dataPath/GM12878_rep4_chr1.binPairs $dataPath/GM12878_rep6_chr1.binPairs -c2 $dataPath/A549_rep1_chr1.binPairs $dataPath/A549_rep3_chr1.binPairs $dataPath/A549_rep4_chr1.binPairs -i $dataPath/GM12878_rep4_chr1.fragFreq $dataPath/GM12878_rep6_chr1.fragFreq -d 2 3 -f $dataPath/MboI_resfrag_hg19_chr1.bed -r $resolution -t $thres -sp $perc -n $nei -o $outPath
