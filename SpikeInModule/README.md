# FreeHi-C Spike-in Simulations for Benchmarking Differential Chromatin Interaction Detection

Ye Zheng, Peigen Zhou, Sündüz Keleş. [FreeHi-C Spike-in Simulations for Benchmarking Differential Chromatin Interaction Detection]() bioRxiv (2020).

The pipeline is developed in Keles Research Group at University of Wisconsin - Madison. Please contact Ye Zheng (yezheng@stat.wisc.edu) or open an issue in the github repository for any question and suggestion. 

## What is FreeHi-C Spike-in Module?

FreeHi-C (**Fr**agment  interactions **e**mpirical **e**stimation for fast simulation of **Hi-C** data) version 2.0 enables spike-ins of differential chromatin interactions that can be data-driven or user-specified. The ability to simulate true differential chromatin interactions allows us to study operating characteristics of state-of-the-art differential interaction detection methods with simulations where both the background and the DCI signals mimic actual data.


FreeHi-C Spike-in module replies on the FreeHi-C training module to process the raw data and obtain the fragment-pair interaction frequency file and bin-pair interaction frequency file. The spike-in module adds data-driven or user-defined spike-ins on the fragment-pair interaction frequency file. Based on this updated fragment-pair interaction frequency file, FreeHi-C simulation module will simulate and get the contact count matrix.

![FreeHi-C Spike-in Module diagram](/figures/FreeHiC_SpikeIn.png)

## FreeHi-C Spike-in Module Usage

### 1. Preparation

In preparation for running spike-in module, apart from the preparation for the FreeHi-C simulation pipeline, you will need to install python3 with corresponding modules required: 

- pandas (0.25.1)
- sklearn (0.19.1)
- statsmodels (0.10.1). 

The versions being tested are specified in the parenthesis.

Subsequently, set the paths to the software executable file in the parameter file (FreeHiC_SpikeIn_parameters) accordingly. Other parameters have been set for the demo data run but they can always be customized for you own usage.


### 2. Parameters explanation

```
spikein          (-s/--spikein)        : Spike-in type. Two options: "data" for data-driven, "user" for user-define. The default is data-driven.
spikeinfile      (-sf/--spikeinfile)   : User defined spike-ins file (tsv format). The format is: chrA, binA, chrB, binB, fold-change. 
                                         Fold change is the interaction frequency of condition2/condition1. No header is needed. If "data-driven" option is selected, this parameter can be left empty.
condition1       (-c1/--condition1)    : A series of bin-pairs interaction files (tsv format) from condition that the spike-ins will be added on. 
                                         For example, spike-ins, determined by comparing A549 and GM12878, are planned to be added on replicates of GM12878. Full path to the bin-pairs from GM12878 should be listed here as condition 1. Multiple files can be separated by space.
condition2       (-c2/--condition2)    : A series of bin-pairs interaction files (tsv format) from condition that the spike-ins signals are calculated. 
                                         For example, spike-ins, determined by comparing A549 and GM12878, are planned to be added on replicates of GM12878. Full path to the bin-pairs from A549 should be listed here as condition 2. Multiple files can be separated by space.
interaction      (-i/--interaction)    : Valid fragment-pairs interaction frequency file(s) (tsv format). A single or a series of files can be input. 
                                         Full path should be given.
index            (-d/--index)          : The index of the valid fragment-pairs interaction frequency file(s) that comes from the same raw sample of the bin-pairs provided in "--condition1". 
                                         A single or a series of index can be input. The length of this parameter should be the same as that of "--interaction".
resfrag          (-f/--resfrag)        : Restriction enzyme fragment file (tsv format). Full path should be given.
outPath          (-o/--outPath)        : Directory to save the spike-in fragment interaction frequency output.
resolution       (-r/--resolution)     : Resolution, i.e. the size of bin.
threshold        (-t/--threshold)      : Threshold of absoluate log2 fold-changes to determine the true differential signals in comparing bin-pairs from condition 1 and condition 2. Default is 2.
percentage       (-sp/--percentage)    : The percentage of spike-ins (0-100) randomly sampled from the true differential set. 100 means the entire bin-pairs in the true set will be added as spike-ins. Default is 100.
number           (-sn/--number)        : [Optional] The number of spike-ins randomly sampled from the true differential set. 
                                         If this number is larger than the total number of bin-pairs in true set, the whole true set will be used.
smooth           (-sm/--smooth)        : [Optional] Indicator for smoothing the neighbors of spike-in bin-pairs. 1 for smooth and 0 for not smooth. 
                                         We recommend smooth and 1 is default value.
neighbor         (-n/--neighbor)       : [Optional] The number of neighbors to be smoothed. Default is 2 meaning 2 bins distance from the target spike-in bin-pairs.
verbose          (-v/--verbose)        : [Optional] Verbose. Default is true.
```

### 3. Running FreeHi-C<a id="sec-1-2-3" name="sec-1-2-3"></a>

- 3.1 Run FreeHi-C pipeline to train the raw data by setting ```train=1```, ```simulate=0```, ```postProcess=0```.

- 3.2 Generate the spike-in fragment-pair interaction frequency file.

```
    bash /path/to/run_FreeHiC_SpikeIn.sh FreeHiC_SpikeIn_parameters
```

- 3.3 Replace the ```outDir/s2_validPairs/repName.validParis.fragFreq``` generated by FreeHi-C training module by this ```*fragFreq.spikeIn```.

- 3.4 Continue to run FreeHi-C pipeline to simulate by setting ```train=0```, ```simulate=1```, ```postProcess=1```.

### 4. Output from FreeHi-C Spike-in Module <a id="sec-1-2-4" name="sec-1-2-4"></a>

Under output directory:
-   results
	- *.fragFreq.spikeIn
	
```*fragFreq.spikeIn``` shares the same prefix file name as the input fragment-pair interaction frequency file(s) through ```--interaction```. Replace the FreeHiC/outDir/s2_validPairs/repName.validParis.fragFreq by this ```*fragFreq.spikeIn``` and continue the simulation module of FreeHi-C, you will get the simulated contact matrices of the spike-in samples.

