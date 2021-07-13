#!/bin/bash

# echo "## ********************************************************"
# echo "## I. Parameters trainning from original Hi-C interactions"
# echo "## ********************************************************"

## parameters
projDir=$1 
fastqFile=$2
outDir=${3:-"${projDir}/results"}
simuName=${4:-"freeHiCsimu"}
ref=$5 

bwa=$6
samtools=$7
bedtools=$8

train=${9:-1}
simulate=${10:-1}
postProcess=${11:-1}
coreN=${12:-1}
mismatchN=${13:-3}
gapN=${14:-1}
mismatchP=${15:-404}
gapP=${16:-404}
chimericP=${17:-404}
simuN=${18}
readLen=${19}
resolution=${20}
lowerBound=${21}
refragU=${22:-800}
refragL=50 #$((seqLength * 2))
ligateSite=${23}
refrag=${24}
summaryFile=${25:-"${outDir}/summary/FreeHiC.summary"}


name=${fastqFile##*/}
mkdir -p $outDir
mkdir -p ${summaryFile%/*}
touch ${summaryFile}

## check basic input parameters
if [[ -z "$projDir" ]]; then
    echo "Please input the path to freeHiC package."
    exit 1
fi
if [[ -z "$fastqFile" ]]; then
    echo "Please input the path to the fastq file."
    exit 1
elif [[ $fastqFile == ftp* ]]
then ## NEED DOWNLOADING! Just for demo data

    mkdir -p ${projDir}/demoData

    wget ${fastqFile}_1.fastq.gz -O ${projDir}/demoData/${name}_1.fastq.gz
    wget ${fastqFile}_2.fastq.gz -O ${projDir}/demoData/${name}_2.fastq.gz

    gunzip ${projDir}/demoData/${name}_1.fastq.gz
    gunzip ${projDir}/demoData/${name}_2.fastq.gz

    fastqFile=${projDir}/demoData/${name}

elif [[ $fastqFile == http* ]]
then ## NEED DOWNLOADING! Just for demo data
    mkdir -p ${projDir}/demoData

    wget -r ${fastqFile}_1.fastq.gz -O  ${projDir}/demoData/${name}_1.fastq.gz
    wget -r ${fastqFile}_2.fastq.gz -O  ${projDir}/demoData/${name}_2.fastq.gz

    gunzip ${projDir}/demoData/${name}_1.fastq.gz
    gunzip ${projDir}/demoData/${name}_2.fastq.gz

    fastqFile=${projDir}/demoData/${name}

fi

if [[ -z "$ligateSite" ]]; then
   echo "Please input the sequence for the restriction enzyme cutting site(s)."
   exit 1
fi

## Train the input raw sequencing data
if [[ "$train" -eq "1" ]]; then
    
    ## 1. Alignment
    echo "I. Training: 1 - Alignment!"
    bash $projDir/scripts/alignment.sh "$projDir" "$fastqFile" "$outDir/rawDataTraining/s1_alignment" "$bwa" "$samtools" "$ref" "$coreN" "$mismatchN" "$gapN" "$summaryFile" "${ligateSite[@]}"

    ## 2. Read-ends Pairing
    echo "I. Training: 2 - Joining read-ends!"
    ## if -t 1, training of mismatch, gap and base quality are also implemented
    python3 $projDir/scripts/joinEnds.py -r1 ${outDir}/rawDataTraining/s1_alignment/${name}_1.sam -r2 ${outDir}/rawDataTraining/s1_alignment/${name}_2.sam -o ${outDir}/rawDataTraining/s1_alignment/${name}.sam -t 1 -sf $summaryFile


    ## 3. Validation
    echo "I. Training: 3 - Categorize read-pairs!"
    python3 $projDir/scripts/categorizePairs.py -f ${refrag} -r ${outDir}/rawDataTraining/s1_alignment/${name}.sam -o ${outDir}/rawDataTraining/s2_validPairs -l $refragL -u $refragU -d $lowerBound -m "window" -b $resolution -sf $summaryFile

    ## 4. Remove duplicates
    echo "I. Training: 4 - Remove duplicates!"
    bash $projDir/scripts/removeDups.sh "${outDir}/rawDataTraining/s2_validPairs/${name}.validPairs" "0" "$summaryFile" "1"

    ## 5. Bining
    echo "I. Training: 5 - Binning!"
    mkdir -p "${outDir}/rawDataTraining/s3_binPairs"
    mkdir -p ${outDir}/sorttmp
    awk '{print $2, $6, $7, $11}' ${outDir}/rawDataTraining/s2_validPairs/${name}.validPairs.nodup | sort -T ${outDir}/sorttmp | uniq -c | awk -v OFS="\t" '{print $2, $3, $4, $5, $1}' | sort -k1,1V -k2,2n -k3,3V -k4,4n -T ${outDir}/sorttmp >${outDir}/rawDataTraining/s3_binPairs/${name}.binPairs
    rm -rf ${outDir}/sorttmp

fi

echo "## *****************************************"
echo "## II. Simulating interaction read sequences"
echo "## *****************************************"

if [[ "$simulate" = "1" ]]; then
    echo "II. Simulating!"
    mkdir -p ${outDir}/${simuName}
    python3 $projDir/scripts/freeHiC_main.py -f ${refrag} -i "${outDir}/rawDataTraining/s2_validPairs/${name}.validPairs.fragFreq" -o "${outDir}/${simuName}/simuSequence" -fn "${name}" -n "${simuN}" -m "${mismatchP}" -idl "${gapP}" -c "${chimericP}" -r "$readLen" -d "$refragU" -s "${bedtools}" -g "${ref%%.gz}" -sf "$summaryFile"
    
    rm -rf ${outDir}/${simuName}/simuSequence/${name}.readPos*
    rm -rf ${outDir}/${simuName}/simuSequence/${name}.readqsOri*
    rm -rf ${outDir}/${simuName}/simuSequence/${name}.readSeq*
    
    if [[ "$postProcess" -eq "1" ]]; then

	## 1. Alignment
	echo "II. Simulation Post-processing: 1 - Alignment!"
	echo "\nSimulation Post-processing!\n" >>$summaryFile
	
	fastqFile="${outDir}/${simuName}/simuSequence/${name}"

	bash $projDir/scripts/alignment.sh "$projDir" "$fastqFile" "${outDir}/${simuName}/simuProcess/s1_alignment" "$bwa" "$samtools"  "$ref" "$coreN" "$mismatchN" "$gapN" "$summaryFile" "0"


	## match sam queryname from both ends %%%% NEW SECTION %%%%%
	$samtools view ${outDir}/${simuName}/simuProcess/s1_alignment/${name}_1.sam | awk '{print $1}' >${outDir}/${simuName}/simuProcess/s1_alignment/queryName_1
	$samtools view ${outDir}/${simuName}/simuProcess/s1_alignment/${name}_2.sam | awk '{print $1}' >${outDir}/${simuName}/simuProcess/s1_alignment/queryName_2
	cat ${outDir}/${simuName}/simuProcess/s1_alignment/queryName_1 ${outDir}/${simuName}/simuProcess/s1_alignment/queryName_2 | sort | uniq -c | awk '$1 == 2 {print $2}' >${outDir}/${simuName}/simuProcess/s1_alignment/queryName.common
	queryLen1=$(cat ${outDir}/${simuName}/simuProcess/s1_alignment/queryName_1 | wc -l)
	queryLen2=$(cat ${outDir}/${simuName}/simuProcess/s1_alignment/queryName_2 | wc -l)
	queryLenC=$(cat ${outDir}/${simuName}/simuProcess/s1_alignment/queryName.common | wc -l)

	if [[ "$queryLen1" != "$queryLenC" ]] || [[ "$queryLen2" != "$queryLenC" ]] || [[ "$queryLen1" != "$queryLen2" ]]; then
	    echo "Matching simulation query names between end 1 and end2......"
	    ${samtools} view -H ${outDir}/${simuName}/simuProcess/s1_alignment/${name}_1.sam >${outDir}/${simuName}/simuProcess/s1_alignment/${name}_1.sam.header
	    ${samtools} view -H ${outDir}/${simuName}/simuProcess/s1_alignment/${name}_2.sam >${outDir}/${simuName}/simuProcess/s1_alignment/${name}_2.sam.header

	    
	    awk 'NR==FNR{a[$0]=1;next}a[$1]' ${outDir}/${simuName}/simuProcess/s1_alignment/queryName.common ${outDir}/${simuName}/simuProcess/s1_alignment/${name}_1.sam  >>${outDir}/${simuName}/simuProcess/s1_alignment/${name}_1.sam.header
	    awk 'NR==FNR{a[$0]=1;next}a[$1]' ${outDir}/${simuName}/simuProcess/s1_alignment/queryName.common ${outDir}/${simuName}/simuProcess/s1_alignment/${name}_2.sam  >>${outDir}/${simuName}/simuProcess/s1_alignment/${name}_2.sam.header
	    mv  ${outDir}/${simuName}/simuProcess/s1_alignment/${name}_1.sam.header ${outDir}/${simuName}/simuProcess/s1_alignment/${name}_1.sam
	    mv  ${outDir}/${simuName}/simuProcess/s1_alignment/${name}_2.sam.header ${outDir}/${simuName}/simuProcess/s1_alignment/${name}_2.sam

	fi

	rm -rf ${outDir}/${simuName}/simuProcess/s1_alignment/queryName* 
	
	## sort SAM file by read id
	${samtools} sort -n -@ $coreN -o ${outDir}/${simuName}/simuProcess/s1_alignment/${name}_1.sorted.sam ${outDir}/${simuName}/simuProcess/s1_alignment/${name}_1.sam
	${samtools} sort -n -@ $coreN -o ${outDir}/${simuName}/simuProcess/s1_alignment/${name}_2.sorted.sam ${outDir}/${simuName}/simuProcess/s1_alignment/${name}_2.sam
	rm -rf ${outDir}/${simuName}/simuProcess/s1_alignment/${name}_1.sam
	rm -rf ${outDir}/${simuName}/simuProcess/s1_alignment/${name}_2.sam
	mv ${outDir}/${simuName}/simuProcess/s1_alignment/${name}_1.sorted.sam ${outDir}/${simuName}/simuProcess/s1_alignment/${name}_1.sam
	mv ${outDir}/${simuName}/simuProcess/s1_alignment/${name}_2.sorted.sam ${outDir}/${simuName}/simuProcess/s1_alignment/${name}_2.sam
	
	## 2. Read-ends Pairing
	echo "II. Simulation Post-processing: 2 - Joining read-ends!"
	## if -t 1, training of mismatch, gap and base quality are also implemented
	python3 $projDir/scripts/joinEnds.py -r1 ${outDir}/${simuName}/simuProcess/s1_alignment/${name}_1.sam -r2 ${outDir}/${simuName}/simuProcess/s1_alignment/${name}_2.sam -o ${outDir}/${simuName}/simuProcess/s1_alignment/${name}.sam -t 0 -sf $summaryFile
	
	
	## 3. Validation
	echo "II. Simulation Post-processing: 3 - Categorize read-pairs!"
	python3 $projDir/scripts/categorizePairs.py -f ${refrag} -r ${outDir}/${simuName}/simuProcess/s1_alignment/${name}.sam -o ${outDir}/${simuName}/simuProcess/s2_validPairs -l $refragL -u $refragU -d $lowerBound -m "window" -b $resolution -sf $summaryFile
	
	## 4. Binning
	echo "II. Simulation Post-processing: 4 - Binning!"
	
	mkdir -p ${outDir}/${simuName}/simuProcess/s3_binPairs
	mkdir -p ${outDir}/${simuName}/sorttmp
	
	awk '{print $2, $6, $7, $11}' ${outDir}/${simuName}/simuProcess/s2_validPairs/${name}.validPairs | sort -T ${outDir}/${simuName}/sorttmp | uniq -c | awk -v OFS="\t" '{print $2, $3, $4, $5, $1}' | sort -k1,1V -k2,2n -k3,3V -k4,4n -T ${outDir}/${simuName}/sorttmp >${outDir}/${simuName}/simuProcess/s3_binPairs/${name}.binPairs
	
	rm -rf ${outDir}/${simuName}/sorttmp
    fi
fi



