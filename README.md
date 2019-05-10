# FreeHi-C: high fidelity Hi-C data simulation for benchmarking and data augmentation

Ye Zheng and Sündüz Keleş. FreeHi-C: high fidelity Hi-C data simulation for benchmarking and data augmentation. BioRxiv 2019.

The pipeline is developed in Keles Research Group in University of Wisconsin - Madison and please contact Ye Zheng (yezheng@stat.wisc.edu) or open an issue in the github repository for any question and suggestion. 

## What is FreeHi-C?

FreeHi-C is short for **Fr**agment  interactions **e**mpirical **e**stimation for fast simulation of **Hi-C** data. It is a data-driven Hi-C data simulator for simulating and augmenting Hi-C datasets.  FreeHi-C employs a non-parametric strategy for estimating an interaction distribution of genome fragments and simulates Hi-C reads from interacting fragments. Data from FreeHi-C exhibit higher fidelity to the biological Hi-C data. FreeHi-C not only can be used to study and benchmark a wide range of Hi-C analysis methods but also boosts power and enables false discovery rate control for differential interaction detection algorithms through data augmentation.

![FreeHi-C pipeline diagram](/figures/FreeHiC_pipeline.png)

## FreeHi-C Usage

### 1. Preparation

    git clone https://github.com/yezhengSTAT/FreeHiC

FreeHi-C installation is finished once you successsfully git clone the repository. The default setting of this FreeHi-C pipeline is a demo run utilizing a real but small Hi-C data: [Plasmodium falciparum genome Trophozoites stage](https://noble.gs.washington.edu/proj/plasmo3d). The raw input data will be downloaded automatically. In preparation for such run, you will need to install
-   BWA: [BWA installation](http://bio-bwa.sourceforge.net/)  (>=0.5.9)
-   samtools: [samtools installation](http://samtools.sourceforge.net/) (>=1.3)
-   bedtools: [bedtools installation](https://bedtools.readthedocs.io/en/stable/content/installation.html) (>=2.27.0)
-   python3 with corresponding modules required: numpy (>= 1.13.1), scipy (>= 0.19.1), pysam (>= 0.12.0), bx-python (>= 0.5.0)

Subsequently, set the software paths in the parameter file accordingly. Other parameters have been set for the demo data run but they can always be customized for you own usage.

### 2. Setting the parameters in the "FreeHiC_parameters'' file

    Path to general folders and necessary genomic files.
    1. projDir        : Path to the FreeHi-C repository.
    2. fastqFile      : Path to the raw squencing fastq file. The two ends fastq file should have been split and use the name format, *_1.fastq and *_2.fastq, to indicate the two ends.
    3. ref            : Path to the reference genome file.
    4. refrag         : Restriction enzyme fragment file.
    5. simuName       : Name of simulation run.
    6. outDir         : Path to the output files.
    7. summaryFile    : Path to the summary file and we suggest that a summary folder is specified first, for example, "${projDir}/summary/${simuName}_FreeHiC.summary".
    
    Path to software.
    8. bwaDir         : Path to the BWA aligner.
    9. samtoolsDir    : Path to the samtools.
    10. bedtoolsDir   : Path to bedtools.
    
    Parameters set for simulation.
    11. train         : Indicator (1 for training or 0 for not training) to specify whether the raw sequencing data need processing and training. For simulation of the same targeted biological replicate but different simulation parameters, if the data has already be trained and the training summary files are ready, there is no need to train again. Default value is 1.
    12. postProcess   : Indicator (1 for processing or 0 for not processing) to specify whether the simulated sequencing reads need processing, namely alignment, joining the read ends, validation and binning. Default value is 1.
    13. coreN         : Number of cores for parallel alignment. Default value is 1.
    14. mismatchN     : Number of mismatches allowed in alignment and simulations. Default value is 3.
    15. gapN          : Number of gaps allowed in alignment and simulations. Default is 1.
    16. mismatchP     : Proportion of simulated reads that have mismatches. Default is the same proportion in the input raw data and input "" to specify the default setting.
    17. gapP          : Proportion of simuated reads that have gaps, including insertions and deletions. Default is the same proportion in the input raw data and input "" to specify the default setting.
    18. chimericP     : Proportion of simulatted reads that are chimeric reads. Default is the same proportion in the input raw data and input "" to specify the default setting.
    19. simuN         : Total number of reads to be simulated. Due to the validation filtering, we recommend input 140-150% of the targeted sequencing depth.
    20. readLen       : Length of the simulated reads.
    21. resolution    : The window size of binning.
    22. lowerBound    : The lower bound for valid long-range interactions. We recommend using 2 times the resolution to remove short-range interactions.
    23. refragU       : Upper bound of the read pair distances summation from read end alignment position to its assigned RE fragment cutting site. Default value is 800.
    24. cutsite       : Restriction enzyme cutting site. This is for rescuing chimeric reads, namely reads that span the restriction enzyme cutting site. For example, "AAGCTAGCTT" for HindIII and "GATCGATC" for MboI. If multiple cutters are utilized, they can be listed in an array: seqLength=("AAGCTAGCTT" "GATCGATC" "..."). If multiple cutting sites occur within one read, the shortest trimmed reads will be kept. If you do not want to rescue the chimeric reads, simply set it to be 0: cutsite=0

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
    
Run time:
- The whole procedure, including raw data processing, simulation and post-processing, should be expected to finish within 4 hours running the demo data using a single-core on a normal computer or server.
