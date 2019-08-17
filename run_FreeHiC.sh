#!/bin/bash
paraFile=$1
source $paraFile

## ==========================
## print the parameter values
## ==========================
echo "projDir=${projDir}"
echo "fastqFile=${fastqFile}"
echo "ref=${ref}"
echo "refrag=${refrag}"
echo "outDir=${outDir}"
echo "simuName=${simuName}"
echo "summaryFile=${summaryFile}"

echo "bwa=${bwa}"
echo "samtools=${samtools}"
echo "bedtools=${bedtools}"

echo "train=${train}"
echo "simulate=${simulate}"
echo "postProcess=${postProcess}"
echo "coreN=$coreN"
echo "mismatchN=${mismatchN}"
echo "gapN=${gapN}"
echo "mismatchP=${mismatchP}"
echo "gapP=${gapP}"
echo "chimericP=${chimericP}"
echo "simuN=${simuN}"
echo "readLen=${readLen}"
echo "resolution=${resolution}"
echo "lowerBound=${lowerBound}"
echo "refragU=${refragU}"
echo "ligateSite=${ligateSite}"

## =========================
## Check the input arguments
## =========================
echo "Pre-checking the input parameters......"
if [[ -z "$projDir" ]]; then
    echo "Please provide the path to FreeHi-C package, namely the path to the FreeHi-C repository after downloanding or cloning."
    exit 1
elif [[ ! -d "$projDir" ]]; then
    echo "The FreeHi-C repository does not exist. Please check the input path!"
    exit 1
else
    echo "  1. FreeHi-C package directory exist."
fi

if [ -z "$fastqFile" ] && [ "$train" = "0" ]; then
    echo "  2. No path to the fastq files is provided and no training is specified. Continue......"
elif [ -z "$fastqFile" ] && [ "$train" = "1" ]; then
    echo "Please provide the path to the fastq files."
    exit 1
elif [ -z "$fastqFile" ] && [ "$train" = "" ]; then
    echo "Please provide the path to the fastq files."
    exit 1
elif [ $fastqFile != ftp* ] && [ $fastqFile == http* ]
then
    if [[ ! -s "${fastqFile}_1.fastq" ]]; then
	echo "${fastqFile}_1.fastq does not exist or is empty!"
	exit 1
    elif [[ ! -s "${fastqFile}_2.fastq" ]]; then
	echo "${fastqFile}_2.fastq does not exist or is empty!"
	exit 1
    else
	echo "  2. ${fastqFile}_1.fastq and ${fastqFile}_2.fastq exist and are not empty."
    fi

else
    echo "  2. ${fastqFile} will be downloaded."
fi

if [[ -z "$ref" ]]; then
    echo "Please provide the path to the reference genome file!"
    exit 1
elif [[ ! -s "$ref" ]]; then
    echo "Reference genome file, $ref, does not exist or is empty!"
    exit 1
else
    echo "  3. Reference genome file exists and is not empty."
fi


if [[ -z "$refrag" ]]; then
    echo "Please provide the path to the restriction enzyme fragment file!"
    exit 1
elif [[ ! -s "$ref" ]]; then
    echo "Restriction enzyme fragment file, $refrag, does not exist or is empty!"
    exit 1
else
    echo "  4. Restriction enzyme fragment file exists and is not empty."
fi

if [[ -z "$simuName" ]]; then
    echo "  5. No name for this simulation run is provided. FreeHi-C will use the default name 'freeHiCsimu' which will be utilized as subfolder name to store all the simulation results. Please check if there are duplicated runs under the same name!"
    
else
    echo "  5. Simulation run name is provided."
fi

if [[ -z "$outDir" ]]; then
    echo "  6. No path to the output directory is provided. FreeHi-C will use the default path, ${projDir}/results!"
else
    echo "  6. Output directory is provided."
fi

if [[ -z "$summaryFile" ]]; then
    echo "  7. No path to the simulation summary file is provided. FreeHi-C will use the default path, ${outDir}/summary/FreeHiC.summary!"
else
    echo "  7. Summary file path is provided."
fi

if [[ -z "$bwa" ]]; then
    echo "Please provide the path to the BWA executable file! If such path has been added to the shell environment variable, PATH variable, please set bwa='bwa'."
    exit 1
else
    echo "  8. Path to BWA executable file is provided. Checking bwa help manual......"
    $bwa
    # if [[ $? = "1" ]]; then
    # 	echo "BWA test run failed! Please check the path to BWA executable file."
    # 	exit 1
    # fi
fi


if [[ -z "$samtools" ]]; then
    echo "Please provide the path to the samtools executable file! If such path has been added to the shell environment variable, PATH variable, please set samtools='samtools'."
    exit 1
else
    echo "  9. Path to samtools executable file is provided. Checking samtools help manual......"
    $samtools --help
    if [[ $? -ne 0 ]]; then
	echo "Samtools test run failed! Please check the path to Samtools executable file."
	exit 1
    fi
fi

if [[ -z "$bedtools" ]]; then
    echo "Please provide the path to the bedtools executable file! If such path has been added to the shell environment variable, PATH variable, please set bedtools='bedtools'."
    exit 1
else
    echo "  10. Path to bedtools executable file is provided. Checking bedtools manual......"
    $bedtools --help
    if [[ $? -ne 0 ]]; then
	echo "Bedtools test run failed! Please check the path to Bedtools executable file."
	exit 1
    fi

fi

if [[ -z "$train" ]]; then
    echo "  11. No action for deciding whether to train the raw input Hi-C sequencing data. FreeHi-C will use the default action, 1, meaning training the raw input data."
elif [ "$train" = "1" ] || [ "$train" = "0" ]; then
    echo "  11. Action is defined for training raw input Hi-C sequencing data."
else
    echo "Please provide 1 for training raw input Hi-C sequencing data or 0 for not training."
    exit 1
fi

if [[ -z "$simulate" ]]; then
    echo "  12. No action for deciding whether to simulate. FreeHi-C will use the default action, 1, meaning simulate."
    simulate=1
elif [ "$simulate" = "1" ] || [ "$simulate" = "0" ]; then
    echo "  12. Action is defined for simulate Hi-C sequencing data."
else
    echo "Please provide 1 for simulate Hi-C sequencing data or 0 for not simulate."
    exit 1
fi

if [[ "$simulate" = "1" ]]; then
    if [[ -z "$postProcess" ]]; then
	echo "  13. No action for deciding whether to process the simulated Hi-C sequencing data. FreeHi-C will use the default action, 1, meaning process the simulated Hi-C sequencing data."
    elif [ "$postProcess" = "1" ] || [ "$postProcess" = "0" ]; then
	echo "  13. Action is defined for whether to process simulated Hi-C sequencing data."
    else
	echo "Please provide 1 for processing simulated Hi-C sequencing data or 0 for not processing."
	exit 1
    fi
    
    if [[ -z "$coreN" ]]; then
	echo "  14. No number of processing core is specified. FreeHi-C will use the default number 1."
    elif [[ $coreN =~ ^[\-0-9]+$ ]] && [[ $coreN -gt 0 ]]; then
	echo "  14. Number of processing core is specified."
    else
	echo "Please set the coreN to be a positive integer!"
	exit 1
    fi
    
    if [[ -z "$mismatchN" ]]; then
	echo "  15. No number of maximum alignment mismatches allowed is specified. FreeHi-C will use the default number 3."
    elif [[ $mismatchN =~ ^[\-0-9]+$ ]] && [[ $mismatchN -gt 0 ]]; then
	echo "  15. Number of maximum alignment mismatches allowed is specified."
    else
	echo "Please set the mismatchN to be a positive integer!"
	exit 1
    fi

    if [[ -z "$gapN" ]]; then
	echo "  16. No number of maximum alignment gaps allowed is specified. FreeHi-C will use the default number 1."
    elif [[ $gapN =~ ^[\-0-9]+$ ]] && [[ $gapN -gt 0 ]]; then
	echo "  16. Number of maximum alignment gaps allowed is specified."
    else
	echo "Please set the gapN to be a positive integer!"
	exit 1
    fi
    
    if [[ -z "$mismatchP" ]]; then
	echo "  17. No number for the percentage of simulated reads that have mismatches is specified. FreeHi-C will use the default percentage that is equal to that of the input raw Hi-C sequencing data."
    elif [[ $mismatchP =~ ^[\-0-9]+$ ]] && [[ $mismatchP -ge 0 ]] && [[ $mismatchP -le 100 ]]; then
	echo "  17. Number for the percentage of simulated reads that have mismatches is specified."
    else
	echo "Please set mismatchP to be a non-negative integer between 0~100!"
	exit 1
    fi
    
    if [[ -z "$gapP" ]]; then
	echo "  18. No number for the percentage of simulated reads that have gap is specified. FreeHi-C will use the default percentage that is equal to that of the input raw Hi-C sequencing data."
    elif [[ $gapP =~ ^[\-0-9]+$ ]] && [[ $gapP -ge 0 ]] && [[ $gapP -le 100 ]]; then
	echo "  18. Number for the percentage of simulated reads that have gaps is specified."
    else
	echo "Please set gapP to be a non-negative integer between 0~100!"
	exit 1
    fi
    
    if [[ -z "$chimericP" ]]; then
	echo "  19. No number for the percentage of simulated reads that are chimeric reads is specified. FreeHi-C will use the default percentage that is equal to that of the input raw Hi-C sequencing data."
    elif [[ $chimericP =~ ^[\-0-9]+$ ]] && [[ $chimericP -ge 0 ]] && [[ $chimericP -le 100 ]]; then
	echo "  19. Number for the percentage of simulated reads that are chimeric reads is specified."
    else
	echo "Please set chimericP to be a non-negative integer between 0~100!"
	exit 1
    fi
    
    if [[ -z "$simuN" ]]; then
	echo "Please provide the total number of reads to simulate!"
	exit 1
    elif [[ $simuN =~ ^[\-0-9]+$ ]] && [[ $simuN -gt 0 ]]; then
	echo "  20. Sequencing depth of simulation is specified."
    else
	echo "Please set simuN to be a positive integer!"
	exit 1
    fi
    
    if [[ -z "$readLen" ]]; then
	echo "Please provide the read length of the simulated read end! This number cannot exceed the input raw Hi-C reads."
	exit 1
    elif [[ $readLen =~ ^[\-0-9]+$ ]] && [[ $readLen -gt 0 ]]; then
	echo "  21. Read length of simulation is specified. Please make sure it does not exceed the read length of input raw Hi-C reads."
    else
	echo "Please set readLen to be a positive integer and less or equal to the read length of input raw Hi-C reads!"
	exit 1
    fi
fi

if [[ -z "$resolution" ]]; then
    echo "Please provide the resolution to process the simulated reads! "
    exit 1
elif [[ $resolution =~ ^[\-0-9]+$ ]] && [[ $resolution -gt 0 ]]; then
    echo "  22. Resolution for processing Hi-C reads is specified."
else
    echo "Please set resolution to be a positive integer!"
    exit 1
fi

if [[ -z "$lowerBound" ]]; then
    echo "Please provide the lower bound for the valid long-range interactions. We recommend using 2 times the resolution to remove these short-range interactions "
    exit 1
elif [[ $lowerBound =~ ^[\-0-9]+$ ]] && [[ $lowerBound -gt 0 ]]; then
    echo "  23. Lower bound for valid long-range interactions is specified."
else
    echo "Please set lowerBound to be a positive integer!"
    exit 1
fi

if [[ -z "$refragU" ]]; then
    echo "  24. No upper bound for the read pair distances summation to its assigned restriction enzyme cutting site is set. FreeHi-C will use the default upper bound 800 meaning 800bp! "
elif [[ $refragU =~ ^[\-0-9]+$ ]] && [[ $refragU -gt 0 ]]; then
    echo "  24. Upper bound for the distance to the assigned restriction enyzme cutting site is specified."
else
    echo "Please set the refragU to be a positive integer!"
    exit 1
fi

if [[ -z "$ligateSite" ]]; then
    echo "Please provide the sequences at the ligation sites! Please note, it is not the restriction enzyme cuting site but the sequences at the ligation sites. For example, if the restriction enzyme is HindIII which recognize 'AAGCTT', the ligation site sequences should be 'AAGCTAGCTT' for HindIII. Simlarly, ligateSite='GATCGATC' for MboI."
    exit 1
else
    echo "  25. Sequences at the ligation site is provided."
fi

python3 -c "import numpy"
if [[ $? = "1" ]]; then
    echo "numpy module is not installed for python3!"
    exit 1
fi
python3 -c "import scipy"
if [[ $? = "1"  ]]; then
    echo "scipy module is not installed for python3!"
    exit 1
fi
python3 -c "import pysam"
if [[ $? = "1"  ]]; then
    echo "pysam module is not installed for python3!"
    exit 1
fi
python3 -c "from bx.intervals.intersection import Intersecter, Interval"
if [[ $? = "1"  ]]; then
    echo "bx-python module is not installed for python3!"
    exit 1
fi
python3 -c "import pyximport;pyximport.install()"
if [[ $? = "1"  ]]; then
    echo "cython module is not installed for python3!"
    exit 1
fi

echo "Input parameter checking is successfull!!"
echo "Start FreeHi-C......"
## =================
## Start of FreeHi-C
## =================
bash freeHiC.sh "$projDir" "$fastqFile" "$outDir" "$simuName" "$ref" "$bwa" "$samtools" "$bedtools" "$train" "$simulate" "$postProcess" "$coreN" "$mismatchN" "$gapN" "$mismatchP" "$gapP" "$chimericP" "$simuN" "$readLen" "$resolution" "$lowerBound" "$refragU" "$ligateSite" "$refrag" "$summaryFile"
