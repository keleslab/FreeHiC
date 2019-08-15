# FreeHi-C: high fidelity Hi-C data simulation for benchmarking and data augmentation

Ye Zheng and Sündüz Keleş. FreeHi-C: high fidelity Hi-C data simulation for benchmarking and data augmentation. BioRxiv 2019.

The pipeline is developed in Keles Research Group in University of Wisconsin - Madison and please contact Ye Zheng (yezheng@stat.wisc.edu) or open an issue in the github repository for any question and suggestion. 

## What is FreeHi-C?

FreeHi-C is short for **Fr**agment  interactions **e**mpirical **e**stimation for fast simulation of **Hi-C** data. It is a data-driven Hi-C data simulator for simulating and augmenting Hi-C datasets.  FreeHi-C employs a non-parametric strategy for estimating an interaction distribution of genome fragments and simulates Hi-C reads from interacting fragments. Data from FreeHi-C exhibit higher fidelity to the biological Hi-C data. FreeHi-C not only can be used to study and benchmark a wide range of Hi-C analysis methods but also boosts power and enables false discovery rate control for differential interaction detection algorithms through data augmentation.

FreeHi-C is designed for studies that are prone to simulate Hi-C interactions from the real data and add deviations from the true ones. Therefore, FreeHi-C requires real Hi-C sequencing data (FASTQ format) as input along with user-defined simulation parameters. FreeHi-C will eventually provide the simulated genomics contact counts in a sparse matrix format (BED format) which is compatible with the standard input of downstream Hi-C analysis.

![FreeHi-C pipeline diagram](/figures/FreeHiC_pipeline.png)

## FreeHi-C Usage

### 1. Preparation

    git clone https://github.com/yezhengSTAT/FreeHiC

FreeHi-C installation is finished once you successsfully git clone the repository. The default setting of this FreeHi-C pipeline is a demo run utilizing a real but small Hi-C data: [Plasmodium falciparum genome Trophozoite stage](https://noble.gs.washington.edu/proj/plasmo3d). The raw input data will be downloaded automatically. In preparation for such run, you will need to install
-   BWA: [BWA installation](http://bio-bwa.sourceforge.net/)  (>=0.5.9)
-   samtools: [samtools installation](http://samtools.sourceforge.net/) (>=1.3)
-   bedtools: [bedtools installation](https://bedtools.readthedocs.io/en/stable/content/installation.html) (>=2.27.0)
-   GNU C++ compiler (>= 4.8.1)
-   python3 with corresponding modules required: numpy (>= 1.13.1), scipy (>= 0.19.1), pysam (>= 0.12.0), bx-python (>= 0.5.0), Cython (>= 0.27.3). For quick python module installation, python-requirements.txt is provided in this repository. Run

		pip install -r python-requirements.txt

Subsequently, set the paths to the software executable file in the parameter file (FreeHiC_parameters) accordingly. Other parameters have been set for the demo data run but they can always be customized for you own usage.

#### 1.1 [Alternative] Creating environment using conda.

1. Install [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

2. Build conda environment:

```
conda env create -f FreeHiC_conda_environment.yml
```

3. Active conda environment for FreeHi-C:

```
conda activate FreeHiC
```

4. Test the FreeHi-C enviroment:

```
bwa
bedtools
```

5. For linux, if you get the following error messages:

```
>bedtools                                                                                                                                                          (FreeHiC) 
bedtools: */envs/FreeHiC/bin/../lib/libstdc++.so.6: version `GLIBCXX_3.4.20' not found (required by bedtools)
bedtools: */envs/FreeHiC/bin/../lib/libstdc++.so.6: version `CXXABI_1.3.8' not found (required by bedtools)
bedtools: */envs/FreeHiC/bin/../lib/libstdc++.so.6: version `GLIBCXX_3.4.21' not found (required by bedtools)
```

Try the following to update libgcc:
```
conda install -c anaconda libgcc=5.2.0
```

6. Run FreeHi-C following the next steps.

### 2. Setting the parameters in the "FreeHiC_parameters'' file

    Path to general folders and necessary genomic files.
    1. projDir        : Path to the FreeHi-C repository. 
    2. fastqFile      : Path to the raw squencing fastq file. The two ends fastq file should have been split with matching read ID name and use the name format, *_1.fastq and *_2.fastq, to indicate the two ends, where * indicate the same fastq file prefix name for two ends. Fastq files should not be zipped. For example, in the demo run, the input fastq files are GSM1215593_trimmedAndFiltered-TROPHOZOITES-XL-AGGG-L2_1.fastq and GSM1215593_trimmedAndFiltered-TROPHOZOITES-XL-AGGG-L2_2.fastq. Therefore, the fastqFile = "GSM1215593_trimmedAndFiltered-TROPHOZOITES-XL-AGGG-L2", namely the * part.
    3. ref            : Path to the reference genome file.
    4. refrag         : Restriction enzyme fragment file. The defailed format and generation is provided in section 5 "Creating simulation using your own data".
    5. simuName       : Name of simulation run. Default is "freeHiCsimu".
    6. outDir         : Path to the output files. Default is ${projDir}/results.
    7. summaryFile    : Path to the summary file and we suggest that a summary folder is specified first, for example, the default is "${outDir}/summary/${simuName}_FreeHiC.summary".
    
    Path to software.
    8. bwa            : Path to the executable bwa aligner file, e.g./path/to/folder/where/bwa/installed/bwa. If the /path/to/folder/where/bwa/installed has been added to environmental variable, $PATH, set bwa=bwa.
    9. samtools       : Path to the executable samtools file, e.g./path/to/folder/where/samtools/installed/samtools. If the /path/to/folder/where/samtools/installed has been added to the environmental variable, $PATH, set samtools=samtools.
    10. bedtools      : Path to the executable bedtools file, e.g./path/to/folder/where/bedtools/installed/bin/bedtools. If the /path/to/folder/where/bedtools/installed/bin has been added to the environmental variable, $PATH, set bedtools=bedtools.
    
    Parameters set for simulation.
    11. train         : Indicator (1 for training or 0 for not training) to specify whether the raw sequencing data need processing and training. For simulation of the same targeted biological replicate but different simulation parameters, if the data has already be trained and the training summary files are ready, there is no need to train again. Default value is 1.
    12. simulate      : Indicator (1 for simulating or 0 for not simulating) to specify whether to simulate the Hi-C sequencing reads. Default is 1.
	13. postProcess   : Indicator (1 for processing or 0 for not processing) to specify whether the simulated sequencing reads need processing, namely alignment, joining the read ends, validation and binning. Default value is 1.
    14. coreN         : Number of cores for parallel alignment. Default value is 1.
    15. mismatchN     : Number of mismatches allowed in alignment and simulations. Default value is 3.
    16. gapN          : Number of gaps allowed in alignment and simulations. Default is 1.
    17. mismatchP     : Percentage of simulated reads that have mismatches (0~100). Default is the same percentage in the input raw data and input "" to specify the default setting.
    18. gapP          : Percentage of simuated reads that have gaps, including insertions and deletions (0~100). Default is the same percentage in the input raw data and input "" to specify the default setting.
    19. chimericP     : Percentage of simulatted reads that are chimeric reads (0~100). Default is the same percentage in the input raw data and input "" to specify the default setting.
    20. simuN         : Total number of reads to be simulated. Due to the validation filtering, we recommend input 140-150% of the targeted sequencing depth.
    21. readLen       : Length of the simulated reads.
    22. resolution    : The window size of binning.
    23. lowerBound    : The lower bound for valid long-range interactions. We recommend using 2 times the resolution to remove short-range interactions.
    24. refragU       : Upper bound of the read pair distances summation from read end alignment position to its assigned RE fragment cutting site. Default value is 800.
    25. ligateSite    : Sequences at ligation sites. This is for rescuing chimeric reads, namely reads that span the ligation sites. Please note, it is not the restriction enzyme cuting site but the sequences at the ligation sites. For example, if the restriction enzyme is HindIII which recognize "AAGCTT", the ligation site sequences should be "AAGCTAGCTT" for HindIII. Simlarly, ligateSite = "GATCGATC" for MboI.  If multiple cutters are utilized, they can be listed in an array: ligateSite=("AAGCTAGCTT" "GATCGATC" "..."). If multiple cutting sites occur within one read, the shortest trimmed reads will be kept. If you do not want to rescue the chimeric reads, simply set it to be 0: ligateSite=0
	

### 3. Running FreeHi-C<a id="sec-1-2-3" name="sec-1-2-3"></a>

    bash run_FreeHiC.sh FreeHiC_parameters

### 4. Output from FreeHi-C<a id="sec-1-2-4" name="sec-1-2-4"></a>

Under output directory:
-   results
    -   rawDataTraining (raw input sequencing data processing and training results)
        -   s1_alignment
            -   repName_1.sam
            -   repName_2.sam
            -   repName.sam
        -   s2_validPairs
            -   repName.validPairs
            -   repName.validPairs.nodup
            -   repName.validPairs.fragFreq (fragment interaction frequency)
        -   s3_binPairs
            -   repName.binPairs
    -   nameOfSimulation (simulation results)
        -   simuSequence (simulated sequencing data)
            -   repName_1.fastq
            -   repName_2.fastq
        -   simuProcess (simulated sequences processing)
            -   s1_alignment
                -   repName_1.sam
                -   repName_2.sam
                -   repName.sam
            -   s2_validPairs
                -   repName.validPairs
            -   s3_binPairs
                -   repName.binPairs

Under summary directory:
-   summary
    -   nameOfSimu_FreeHiC.summary (raw data processing and simulation related number summary)
    -   nameOfSimu_FreeHiC.summary.[alignedN/allMatch/baseQS/chimericN/deletion/insertion/mismatch/revStrandN] (raw input sequencing data parameter training results)

Under the directory of input parameter file:
- FreeHiC_parameters (This is your input parameter file)

### 5. Creating simulation with your own data

1. Follow section 1 "Preparation" to download FreeHi-C repository and install essential software.

2. Download or prepare your own Hi-C sequencing data (FASTQ files). Hi-C reads are paired-end reads, therefore there will be two fastq files. The two ends fastq file should have been split with matching read ID name and use the name format, *_1.fastq and *_2.fastq, to indicate the two ends. * indicate the same fastq file prefix name for two ends. Fastq files should already be unzipped. For example, in the demo run, the input fastq files are GSM1215593_trimmedAndFiltered-TROPHOZOITES-XL-AGGG-L2_1.fastq and GSM1215593_trimmedAndFiltered-TROPHOZOITES-XL-AGGG-L2_2.fastq. Therefore, the fastqFile = "GSM1215593_trimmedAndFiltered-TROPHOZOITES-XL-AGGG-L2", namely the * part.

3. Download or prepare the corresponding reference genome (FASTA file) for alignment.

4. Download or prepare the restriction enzyme fragment (BED file). HiCPro provides scripts and detailed instruction to generate such files ([http://nservant.github.io/HiC-Pro/UTILS.html#digest-genome-py](http://nservant.github.io/HiC-Pro/UTILS.html#digest-genome-py)). The restriction enzyme fragment file has six columns and provides genomic coordinates for the starting and ending positions of each fragment.

```
chr1    0       16007   HIC_chr1_1      0       +
chr1    16007   24571   HIC_chr1_2      0       +
chr1    24571   27981   HIC_chr1_3      0       +
chr1    27981   30429   HIC_chr1_4      0       +
chr1    30429   32153   HIC_chr1_5      0       +
chr1    32153   32774   HIC_chr1_6      0       +
chr1    32774   37752   HIC_chr1_7      0       +
chr1    37752   38369   HIC_chr1_8      0       +
chr1    38369   38791   HIC_chr1_9      0       +
chr1    38791   39255   HIC_chr1_10     0       +
chr1    39255   43602   HIC_chr1_11     0       +
```
5. Set the name for this simulation run. Such name will be used to create a subfolder under outDir to save all the simulation results. Default name is "freeHiCsimu" but please make sure every run has its own name. Otherwise, the results will be overlapped.

6. Set the path to the software: bwa, samtools and bedtools. Please note, these path should be set to the software executable file. For example, ```bedtools=/project/tools/bedtools/bin/bedtools```. If the software installation folder has been added to the environmental variable, $PATH, please set ```bedtools=bedtools```.

7. Set action indicators (1 for implement and 0 for not implement) for ```train```, ```simulate```, and ```postProcess``` . For simulation of the same targeted biological replicate but using different simulation parameters (such as different sequencing depths), if the data has already be trained and the training summary files are ready, there is no need to train again. If users want to process the simulated sequencing data in their own way, please set ```postProcess=0``` so that FreeHi-C will stop after generating the simulated FASTQ files.

8. If users have no special requirement for the mismatches, gap and chimeric reads, please denote "" to the parameters, for example, mismatchN="", so that FreeHi-C will use the default setting. The default setting is mismatchN=3, gapN=1, mismatchP, gapP, chimericP will be equal to the same percentage of reads that have mismatches, gaps, are chimeric reads after aligning the input raw Hi-C reads. Please note, mismatchP, gapP, chimericP value range is 0~100.

9. ```simuN``` records the number of reads that users want to simulate, namely the sequencing depth of simulation data. Please note, the final number of simulated valid read pairs will be smaller than this number because after adding noise and deviations to the simulation, not all the simulated reads are alignable and valid read pairs. Please set this number to be ~120% more than your expectation sequencing depth. This number can go beyond the original sequencing depth of input Hi-C data.

10. ```readLen``` records the read length of the simulated read ends. Please set a number that is small or equal to the read length of the input Hi-C read ends.

11. ```resolution``` is the length of window to bin the genome. For high sequencing depth Hi-C data, it can be 5000 or 10000. For low sequencing depth Hi-C data, it can be 40000 or even larger such as 1000000.

12. ```lowerBound``` is the lower bound for defining a valid long-range interaction read pair. Hi-C is designed to study long-range genomic contact however due to the structure of DNA fiber, many random interactions occur between genome regions that are linearly close to each other. These random short-range interactions are not the target of Hi-C studies hence should be excluded from the downstread analysis. We recommend setting it to be 2 times the resolution.

13. ```refragU``` is the upper bound to filter and save valid read pairs. According to the Hi-C experimental protocol, the read ends should originate from regions near the restriction enzyme cutting site. If the distance from the read ends alignment position to the nearest restriction enyzme cutting site is too far to be true. Such read pairs should be discarded. This upper bound is the maximum distance summation from two ends to their own assigned restriction enyzme cutting site. Default is set to be 800 meaning 800bp.

14. ```ligateSite``` records the sequences at ligation sites. This is for rescuing chimeric reads, namely reads that span the ligation sites hence contains genome sequences from two non-adjacent regions. Please note, it is not the restriction enzyme cuting site but the sequences at the ligation sites. For example, if the restriction enzyme is HindIII which recognize "AAGCTT" and have cut at each recognized site, after ligation, the sequence at each ligation site should be "AAGCTAGCTT" for HindIII. Simlarly, ligateSite = "GATCGATC" for MboI.

### 6. Docker image for FreeHi-C

If you have problem with dependencies and software versions, container technology [Docker](https://www.docker.com/) can help address these issue. 

1. Create an account at [Docker Hub](https://hub.docker.com).

2. Install docker: 
   - macOS: Try [Docker Mac Desktop](https://hub.docker.com/editions/community/docker-ce-desktop-mac) first. If failed, try [Docker Toolbox (Mac)](https://docs.docker.com/toolbox/toolbox_install_mac/).
   - Windows: Try [Docker Windows Desktop](https://hub.docker.com/editions/community/docker-ce-desktop-windows) first. If failed, tryp [Docker Toolbox (Windows)](https://docs.docker.com/toolbox/toolbox_install_windows/).
   - Linux: Look for the suitable [Docker Engine](https://hub.docker.com/search/?type=edition&offering=community).
   
3. After installing Docker, test the Docker installation by running 

```
docker run hello-world
```
Expect something like the following:

```
Hello from Docker!
This message shows that your installation appears to be working correctly.
```
4. Pull FreeHi-C Docker image and re-name such image.

```
docker pull yezheng/freehic_docker
docker tag yezheng/freehic_docker freehic_docker
```
5. Modify the parameter values in the parameter file for docker run ```FreeHiC_parameters_docker```.

6. Run FreeHi-C Docker image. 
   -  If all the input data, including raw sequencing data (fastq files), reference genome (fasta file) and restriction fragment file (bed file), are saved within the same data folder.
   
	   ```
	   docker run -v "/path/to/parameter/file/FreeHiC_parameters_docker:/FreeHiC/FreeHiC_parameters" -v "/path/to/input/data/folder:/FreeHiC/data" -v "/path/to/results/folder:/FreeHiC/results" freehic_docker bash run_FreeHiC.sh FreeHiC_parameters
	   ```

		where ```-v "/path/to/parameter/file/FreeHiC_parameters_docker:/FreeHiC/FreeHiC_parameters"``` is to pass user defined parameter file, ```FreeHiC_parameters_docker```, to the Docker container.

		```-v "/path/to/input/data/folder:/FreeHiC/data"``` is to tell the FreeHi-C Docker container that data folder is "/path/to/input/data/folder". 
		```-v "/path/to/results/folder:/FreeHiC/results" ``` is to tell the FreeHi-c Docker container that results should be saved at "/path/to/results/folder".

   - If input data are saved under different path, the absolute path to the input file should be given. The file name should be consistent with that in the parameter file.
   
	   ```
	   docker run -v "/path/to/parameter/file/FreeHiC_parameters_docker:/FreeHiC/FreeHiC_parameters" -v "/path/to/demoData/demoRep_1.fastq:/FreeHiC/data/demoRep_1.fastq" -v "/path/to/demoData/demoRep_2.fastq:/FreeHiC/data/demoRep_2.fastq" -v "/path/to/ref.fasta:/FreeHiC/data/ref.fasta" -v "/path/to/restrictFrag.bed:/FreeHiC/data/restrictFrag.bed" -v "/path/to/results:/FreeHiC/results" freehic_docker bash run_FreeHiC.sh FreeHiC_parameters
	   ```
All the raw data processing and simulation results are saved at ```/path/to/results```.

### 7. Parallel running for large Hi-C data

1. Parallel processing and training of input Hi-C sequencing data:

FreeHi-C utilizes BWA to align the input Hi-C data. Users can make use of the multi-threading mode in BWA by setting ```coreN``` to accelerate the alignment process. Besides, if high-throughput computing resources are available, for example [Center for High Throughput Computing (CHTC)](chtc.cs.wisc.edu), users can first split the input Hi-C sequencing fastq file into smaller files and do it to both ends respectively and correspondingly. Then align and train the smaller sequencing file independently. Please note, the ```lineN``` parameter should take a quadruple value.

```
lineN=4*(num of reads in each small file)

split -l $lineN -d -a 3 /path/to/HiC/sequence/data/${name}_1.fastq $outPath/splitHiC_end1/$prefixname\_  ## The smaller sequencing file after split is save as $outPath/splitHiC_end1/${prefixname}_000 $outPath/splitHiC_end1/${prefixname}_001 $outPath/splitHiC_end1/${prefixname}_002 ...


for f in $outPath/splitHiC_end1/*;
do
    mv $f $f\_1.fastq
done
## rename the smaller sequencing file: $outPath/splitHiC_end1/${prefixname}_000_1.fastq $outPath/splitHiC_end1/${prefixname}_001_1.fastq
## do the same processing on read end 2.
```

2. Run FreeHi-C to align and train each split read pairs:

For each pair of smaller Hi-C sequencing data, ```$outPath/splitHiC_end1/${prefixname}_001_1.fastq``` and ```$outPath/splitHiC_end1/${prefixname}_001_2.fastq```, run FreeHi-C training module only by setting ```train=1```, ```simulate=0```, ```postProcess=0``` along with ```fastqFile=$outPath/splitHiC_end1/${prefixname}_001``` etc.

```
bash run_freehic.sh FreeHiC_parameters_onlyTraining
```

3. Merge the input Hi-C data processing and training results:

3.1. Merge the valid read pairs from each small run under ```$outPath/rawDataTraining/*.validPairs.fragFreq```. First collect all the ${prefixname}_00*.validPairs.fragFreq files and save within the same folder. Then merge by grouping by the first and second column and sum the third column. 

```
cat *validParis.fragFreq | sort -k2 -k3 | awk '{print $2, $3, $1}' | awk '{a[$1" "$2]+=$3}END{for (i in a) print i,a[i]}' | sort -k1 -k2 | awk '{print $3, $1, $2}' >merged.validPairs.fragFreq
```

3.2. Merge the training summary output. According to the ```$summaryFile``` value, all the training summary output is save as ```$summaryFile.*```. For example, ```$summaryFile.chimercN``` contains one number which is the total number of chimeric reads found in this split Hi-C input fastq file. To merge, first collect all the summary files and save within the same folder. The detailed merging command can be found under ```scripts/merge.sh```. Please adjust the paths accordingly.

```
bash merge.sh $prefixName $outPath
```

4. Parallel simulation:

Simulation can be paralleled by 1) simulate by chromosome 2) simulate a small number of reads but implement it several times simualtenously.

4.1. Simulate by chromosome. Separate the merged.validPairs.fragFreq by chromosome. FreeHi-C will empirically simulate by this file therefore replace the complete fragment frequency list by such smaller fragment frequency list only involving the target chromosome.

```
awk -v OFS="\t" '$2 == "chr1" && $3 == "chr1" {print $0}' merged.validPairs.fragFreq >merged_chr1.validPairs.fragFreq
```


4.2. Simulate many runs. For example, set ```simuN``` to be (total number of reads to be simulated)/100 but generate 100 FreeHiC_parameters_* files with different simuName. Run simulation and postProcess module for each run. Then merge the simulated bin-pairs.

```
cat *.binPairs | awk '{a[$1" "$2" "$3" "$4]+=$5}END{for (i in a) print i,a[i]}' | sort >merge.binPairs
```

### 8. Run time
- The whole procedure, including raw data processing, simulation and post-processing, should be expected to finish within 4 hours running the demo data using a single-core on a normal computer or server. The following are reference summary of runtime and memory using GM12878, A549 and P.falciparum with respect to the sequencing depth and number of running cores. For more details, please refer to the manuscript.

![Plasmodium_runtime1](/figures/Plasmodium_runtime1.png)
![Plasmodium_runtime4](/figures/Plasmodium_runtime4.png)
![GM12878_runtime](/figures/GM12878_runtime.png)
![A549_runtime](/figures/A549_runtime.png)
![A549_memory](/figures/memory.png)
