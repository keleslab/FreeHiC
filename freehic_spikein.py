import numpy as np
import pandas as pd
import argparse
import os, sys
from sklearn.neighbors import KDTree
from statsmodels.nonparametric.kernel_regression import KernelReg

def get_args():
    '''get arguments'''
    parser = argparse.ArgumentParser(description = '------------Usage Start------------',
                                     epilog = '------------Usage End------------')
    parser.add_argument('-s', '--spikein', help = 'Type of determining the spike-in. Two options: "data" for data-driven, "user" for user-define. The default is data-driven.', default = "data")
    parser.add_argument('-sf', '--spikeinfile', help = 'User defined spike-ins file (tsv format). The format is: chrA, binA, chrB, binB, fold-change. Fold change is the interaction frequency of condition2/condition1. No header is needed. If "data-driven" option is selected, this parameter can be left empty.', default = None)
    parser.add_argument('-c1', '--condition1', help = 'A series of bin-pairs interaction files (tsv format) from condition that the spike-ins will be added on. For example, spike-ins, determined by comparing A549 and GM12878, are planned to be added on replicates of GM12878. Full path to the bin-pairs from GM12878 should be listed here as condition 1. Multiple files can be separated by space.', nargs='+', required = True, default = None)
    parser.add_argument('-c2', '--condition2', help = 'A series of bin-pairs interaction files (tsv format) from condition that the spike-ins signals are calculated. For example, spike-ins, determined by comparing A549 and GM12878, are planned to be added on replicates of GM12878. Full path to the bin-pairs from A549 should be listed here as condition 2. Multiple files can be separated by space.', nargs='+', default = None)
    parser.add_argument('-i', '--interaction', help = 'Valid fragment-pairs interaction frequency file(s) (tsv format). A single or a series of files can be input. Full path should be given.', default = None, nargs='+', required = True)
    parser.add_argument('-d', '--index', help = 'The index of the valid fragment-pairs interaction frequency file(s) that comes from the same raw sample of the bin-pairs provided in "--condition1". A single or a series of index can be input. The length of this parameter should be the same as that of "--interaction".', default = None, nargs='+', required = True)
    parser.add_argument('-f', '--resfrag', help = 'Restriction enzyme fragment file (tsv format). Full path should be given.', default = None, required = True)
    parser.add_argument('-r', '--resolution', help = 'Resolution, i.e. the size of bin.', default = None, required = True)
    parser.add_argument('-t', '--threshold', help = 'Threshold of absoluate log2 fold-changes to determine the true differential signals in comparing bin-pairs from condition 1 and condition 2. Default is 2.', default = 2)
    parser.add_argument('-sp', '--percentage', help = 'The percentage of spike-ins (0-100) randomly sampled from the true differential set. 100 means the entire bin-pairs in the true set will be added as spike-ins.', default = 100)
    parser.add_argument('-sn', '--number', help = 'The number of spike-ins randomly sampled from the true differential set. If this number is larger than the total number of bin-pairs in true set, the whole true set will be used.', default = None)
    parser.add_argument('-sm', '--smooth', help = 'Indicator for smoothing the neighbors of spike-in bin-pairs. 1 for smooth and 0 for not smooth. We recommend smooth and 1 is default value.', default = 1)
    parser.add_argument('-n', '--neighbor', help = 'The number of neighbors to be smoothed. Default is 2 meaning 2 bins distance from the target spike-in bin-pairs.', default = 2)
    parser.add_argument('-o', '--outPath', help = 'Directory to save the spike-in fragment interaction frequency output.', default = None, required = True)
    parser.add_argument('-v', '--verbose', help = '(Optional) Verbose. Default is true.', default = True)

    args = parser.parse_args()

    # print(args)
    
    if args.spikein != "data" and args.spikein != "user":
        print("There are two options: 'data' for data-driven, 'user' for user-defined. Please choose from these two!")
        sys.exit()
        
    if args.spikeinfile is not None:
        if not os.path.exists(args.spikeinfile):
            print(args.spikeinfile + "does not exist!")
            sys.exit()    
                
    for condFile in args.condition1:
        if not os.path.exists(condFile):
            print(condFile + "from condition 1 does not exist!")
            sys.exit()
    
    if args.condition2 is None:
        if args.spikein == "data":
            print("Please provide bin-pairs for condition2!")
            sys.exit()
    else:
        for condFile in args.condition2:
            if not os.path.exists(condFile):
                print(condFile + "from condition 2 does not exist!")
                sys.exit()
            
    for interFile in args.interaction:
        if not os.path.exists(interFile):
            print(interFile + "from condition 1 does not exist!")
            sys.exit()
    
    if len(args.interaction) != len(args.index):
        print("The length of --interaction and --index does not match!")
        sys.exit()
    
    if not os.path.exists(args.resfrag):
        print("Restriction enzyme fragment file does not exist!")
        sys.exit()
        
    if float(args.threshold) < 2:
        print("Warning: We recommend setting the threshold to determine the true differential bin-pairs using at least 2!!!å")
        
        
    if float(args.percentage) > 100 or float(args.percentage) < 0:
        print("Please provide valid spike-in percentage which ranges from 0 to 100.")
        sys.exit()
    elif float(args.percentage) < 1:
        print("Please note that the percentage data range is from 0 to 100. Current percentage is " + str(args.percentage))
        
    if args.number is not None:
        print(" '--number' is activiated. '--percentage' parameter will be ignored. If this number is larger than the size of true set, the entire true set will be used for spike-in.")
        
    if not os.path.exists(args.outPath):
        print("Directory to save the output does not exist. Creating one......")
        try:
            os.mkdir(args.outPath)
        except OSError:
            print ("Creation of the directory %s failed" % path)
            sys.exit()
        else:
            print ("Successfully created the directory %s " % path)
    
    return args



if __name__ == "__main__":

    args = get_args()
    if args.verbose:
        print("Type of determining spike-ins = ", args.spikein)
        print("List of bin-pairs for condition1 = ", args.condition1)
        print("List of bin-pairs for condition2 = ", args.condition2)
        print("Fragment interaction frequency file = ", args.interaction)
        print("Index of fragment file corresponds to the bin-pairs list of condition 1 = ", args.index)
        print("Restriction enzymen fragment file = ", args.resfrag)
        print("Resolution = ", args.resolution)
        print("Filtering threshold = ", args.threshold)
        if args.number is None:
            print("Percentage of spike-ins = ", args.percentage)
        else:
            print("Number of spike-ins = ", args.number)
        print("Number of neighbors to be smoothed = ", args.neighbor)
        print("Directory to save the spike-in output = ", args.outPath)
    
    spikeType = args.spikein
    cond1 = args.condition1
    cond2 = args.condition2
    interFilePaths = args.interaction
    condInds = args.index
    fragPath = args.resfrag
    resolution = int(args.resolution)
    thres = float(args.threshold)
    perc = float(args.percentage)
    if args.number is not None:
        number = int(args.number)
    else:
        number = None
    smooth = int(args.smooth)
    nei = int(args.neighbor)
    outPath = args.outPath
    
    
    radius = nei + 2
    radValid = nei

    keyCoord = ["chrA", "binA", "chrB", "binB"]
    keyC1 = ["cond1_" + str(i) for i in range(1, 5)]
    keyC2 = ["cond2_" + str(i) for i in range(1, 5)]

    rawFile = None
    ind = 1
    for condFile in cond1:
        if rawFile is None:
            rawFile = pd.read_csv(condFile, sep = '\t', header = None, 
                                names = ["chrA", "binA", "chrB", "binB", "cond1_" + str(ind)])
        else:
            rawFileTmp = pd.read_csv(condFile, sep = '\t', header = None, 
                                   names = ["chrA", "binA", "chrB", "binB", "cond1_" + str(ind)])
            rawFile = pd.merge(rawFile, rawFileTmp, how = "outer", on = ["chrA", "binA", "chrB", "binB"]).fillna(0)
        ind = ind + 1

    if spikeType == "data":
        ind = 1
        for condFile in cond2:
            rawFile = pd.merge(rawFile, pd.read_csv(condFile, sep = '\t', header = None, 
                                                names = ["chrA", "binA", "chrB", "binB", "cond2_" + str(ind)]), how = "outer", on = ["chrA", "binA", "chrB", "binB"]).fillna(0)
            ind = ind + 1



        ## normalize by the sequencing depth
        seqDep = rawFile[["cond1_" + str(i) for i in range(1, 5)] + ["cond2_" + str(i) for i in range(1, 5)]].sum(axis = 0)
        normFile = rawFile[["cond1_" + str(i) for i in range(1, 5)] + ["cond2_" + str(i) for i in range(1, 5)]]/seqDep * seqDep.max()

        ## Filter and get the full true set list
        medFile = pd.concat([rawFile[keyCoord], normFile[keyC1], normFile[keyC1].median(axis=1).rename("condition1"), normFile[keyC2].median(axis=1).rename("condition2")], axis=1).copy()
        medFile[["condition1", "condition2"]] = medFile[["condition1", "condition2"]].replace(0, 0.45)
        medFile["log2FC"] = np.log2(medFile.condition2/medFile.condition1)

        chromeSome = medFile.chrA.unique()
        chromeDist = {chrome: (i + 1) * resolution * radius * 10 for i, chrome in enumerate(chromeSome)}
        medFile.loc[:, 'chromDist'] = medFile.chrA.replace(chromeDist)

        trueSet = medFile[medFile.log2FC.abs() > 2]

        ## Get the spike-in bin-pairs list
        if number is not None:
            if number > trueSet.shape[0]:
                spSelect = trueSet
            else:
                spSelect = trueSet.sample(number)
        elif perc < 100:
            spSelect = trueSet.sample(frac = perc/100)
        else:
            spSelect = trueSet

    elif spikeType == "user":

        ## normalize by the sequencing depth
        seqDep = rawFile[["cond1_" + str(i) for i in range(1, 5)]].sum(axis = 0)
        normFile = rawFile[["cond1_" + str(i) for i in range(1, 5)]]/seqDep * seqDep.max()

        ## Filter and get the full true set list
        medFile = pd.concat([rawFile[keyCoord], normFile[keyC1], normFile[keyC1].median(axis=1).rename("condition1")], axis=1).copy()

        chromeSome = medFile.chrA.unique()
        chromeDist = {chrome: (i + 1) * resolution * radius * 10 for i, chrome in enumerate(chromeSome)}
        medFile.loc[:, 'chromDist'] = medFile.chrA.replace(chromeDist)

        spikeinFile = pd.read_csv(args.spikeinfile, sep = '\t', header = None, names = ["chrA", "binA", "chrB", "binB", "fc"])
        spikeinFile["condition1"] = 1
        spikeinFile["condition2"] = spikeinFile["fc"]

        spSelect = pd.merge(spikeinFile, medFile.reset_index(), how = "left", on = keyCoord).set_index('index').fillna(0)


    ## Build KDTree for fast neighbor bin-pairs extraction
    ## TODO: Add chrA checking, Only same chromesome can applied this
    ## Solve: Group by chrA and chrB, use a dictionary to save KDTree
    if smooth == 1:
        distVars = ['chromDist', 'binA', 'binB']
        X = medFile[distVars].values
        tree = KDTree(X, metric="manhattan")
        inds, dists = tree.query_radius(spSelect[distVars], r=resolution * radius, return_distance=True, sort_results=True)

        ## Spike-in and smoothing ---!!!! super slow
        for i in range(0, spSelect.shape[0]): #
            ## Spike-ins
            fc = spSelect.condition2.iloc[i]/spSelect.condition1.iloc[i]
            target = medFile.loc[spSelect.index[i], keyC1]
            if fc >1:
                target = target.replace(0, 0.45) * fc
            else:
                target *= fc
            medFile.loc[spSelect.index[i], keyC1] = target

            ## Smooth neighbors of target bin-pairs for spike-ins
            x = dists[i].reshape(-1, 1)
            updateInd = inds[i][x.reshape(inds[i].shape) <= resolution *radValid]
            for keyTmp in ["cond1_" + str(i) for i in range(1, 5)]:
                y = medFile.loc[inds[i], keyTmp].values.reshape(-1, 1)
                smOBJ = KernelReg(y, x, var_type = ["c"], bw = [resolution/2])
                yhat, xhat = smOBJ.fit(np.linspace(x[0], x[-1], num = len(x)))
                medFile.loc[updateInd, keyTmp] = yhat[x.reshape(inds[i].shape) <= resolution *2]

    ## Reinforce the spike-in signals 
    rate = spSelect.condition2/spSelect.condition1
    medFile.loc[spSelect.index, keyC1] = normFile.loc[spSelect.index, keyC1].mul(rate, axis = 0)

    ## Convert bin-pairs changes into fragment pairs changes
    fragPath="/p/keles/yezheng/volumeA/HiC_essentialData/MboI_resfrag_hg19_short.bed"
    fragFile = pd.read_csv(fragPath, sep = '\t', header = None, 
                        names = ["chrom", "start", "end", "name", "label", "strand"])

    fragFile["binS"] = np.floor(fragFile.start/resolution) * resolution + resolution/2
    fragFile["binE"] = np.floor(fragFile.end/resolution) * resolution + resolution/2

    
    for i in range(0, len(interFilePaths)):
        interFilePath = interFilePaths[i]
        interFileName = os.path.basename(interFilePath)
        print("Processing " + interFileName)
        # interFilePath = "/p/keles/scrna-seq/volumeC/freeHiC/GM12878/rep2/chr1/s1_training/validPairs/rep2_chr1.validPairs.fragFreq.tmp"

        interFreq = pd.read_csv(interFilePath, sep = "\t", header = None, names = ["freq", "fragA", "fragB"])

        ## Adjust fragment pairs which span multiple bins - split the fragment interaction frequency
        fragF = pd.merge(pd.merge(interFreq, fragFile[["name", "chrom", "binS", "binE"]], how = "left", left_on = "fragA", right_on = "name"), fragFile[["name", "chrom", "binS", "binE"]], how = "left", left_on = "fragB", right_on = "name")
        idx1 = (fragF.binS_x == fragF.binE_x) & (fragF.binS_y == fragF.binE_y)
        idx2 = (fragF.binS_x != fragF.binE_x) & (fragF.binS_y == fragF.binE_y)
        idx3 = (fragF.binS_x == fragF.binE_x) & (fragF.binS_y != fragF.binE_y)
        idx4 = (fragF.binS_x != fragF.binE_x) & (fragF.binS_y != fragF.binE_y)
        fragF1Clean = fragF.loc[idx1, ["freq", "fragA", "fragB", "chrom_x", "binS_x", "chrom_y", "binS_y"]].copy().rename(columns = {"chrom_x":"chrA", "binS_x":"binA", "chrom_y":"chrB", "binS_y":"binB"})

        fragF2 = fragF[idx2].copy()
        fragF2tmp1 = pd.concat([pd.DataFrame(fragF2["freq"]/2), fragF2[["fragA", "fragB", "chrom_x", "binS_x", "chrom_y", "binS_y"]].rename(columns = {"chrom_x":"chrA", "binS_x":"binA", "chrom_y":"chrB", "binS_y":"binB"})], axis = 1)
        fragF2tmp2 = pd.concat([pd.DataFrame(fragF2["freq"]/2), fragF2[["fragA", "fragB", "chrom_x", "binE_x", "chrom_y", "binS_y"]].rename(columns = {"chrom_x":"chrA", "binE_x":"binA", "chrom_y":"chrB", "binS_y":"binB"})], axis = 1)
        fragF2Clean = fragF2tmp1.append(fragF2tmp2)


        fragF3 = fragF[idx3].copy()
        fragF3tmp1 = pd.concat([pd.DataFrame(fragF3["freq"]/2), fragF3[["fragA", "fragB", "chrom_x", "binS_x", "chrom_y", "binS_y"]].rename(columns = {"chrom_x":"chrA", "binS_x":"binA", "chrom_y":"chrB", "binS_y":"binB"})], axis = 1)
        fragF3tmp2 = pd.concat([pd.DataFrame(fragF3["freq"]/2), fragF3[["fragA", "fragB", "chrom_x", "binS_x", "chrom_y", "binE_y"]].rename(columns = {"chrom_x":"chrA", "binS_x":"binA", "chrom_y":"chrB", "binE_y":"binB"})], axis = 1)
        fragF3Clean = fragF3tmp1.append(fragF3tmp2)


        fragF4 = fragF[idx4].copy()
        fragF4tmp1 = pd.concat([pd.DataFrame(fragF4["freq"]/4), fragF4[["fragA", "fragB", "chrom_x", "binS_x", "chrom_y", "binS_y"]].rename(columns = {"chrom_x":"chrA", "binS_x":"binA", "chrom_y":"chrB", "binS_y":"binB"})], axis = 1)
        fragF4tmp2 = pd.concat([pd.DataFrame(fragF4["freq"]/4), fragF4[["fragA", "fragB", "chrom_x", "binS_x", "chrom_y", "binE_y"]].rename(columns = {"chrom_x":"chrA", "binS_x":"binA", "chrom_y":"chrB", "binE_y":"binB"})], axis = 1)
        fragF4tmp3 = pd.concat([pd.DataFrame(fragF4["freq"]/4), fragF4[["fragA", "fragB", "chrom_x", "binE_x", "chrom_y", "binS_y"]].rename(columns = {"chrom_x":"chrA", "binE_x":"binA", "chrom_y":"chrB", "binS_y":"binB"})], axis = 1)
        fragF4tmp4 = pd.concat([pd.DataFrame(fragF4["freq"]/4), fragF4[["fragA", "fragB", "chrom_x", "binE_x", "chrom_y", "binE_y"]].rename(columns = {"chrom_x":"chrA", "binE_x":"binA", "chrom_y":"chrB", "binE_y":"binB"})], axis = 1)
        fragF4Clean = fragF4tmp1.append([fragF4tmp2, fragF4tmp3, fragF4tmp4])

        fragFClean = fragF1Clean.append([fragF2Clean, fragF3Clean, fragF4Clean]).groupby(["fragA", "fragB", "chrA", "binA", "chrB", "binB"]).sum().reset_index()


        ## Compare the old and new bin-pairs interaction frequency
        condInd = condInds[i]
        diffInd = normFile["cond1_" + str(condInd)] != medFile["cond1_" + str(condInd)]
        fcFile = pd.concat([medFile.loc[diffInd, keyCoord + ["cond1_" + str(condInd)]].rename(columns = {"cond1_" + str(condInd): "new"}), 
                            normFile.loc[diffInd, ["cond1_" + str(condInd)]].rename(columns = {"cond1_" + str(condInd): "old"})], axis = 1)

        ## originally there is non-zero fragment interactions
        fragFCleanPos = pd.merge(fragFClean, fcFile[fcFile.old > 0], how = "inner", on = keyCoord)
        fragFCleanPos["freq"] = fragFCleanPos["freq"] * fragFCleanPos["new"] / fragFCleanPos["old"]
        fragFCleanPosFrag = fragFCleanPos[["fragA", "fragB", "freq"]].groupby(["fragA", "fragB"]).sum().reset_index()

        ## originally there is zero fragment interaction - random select one frag-pair and assign the new frequency to it
        fragFCleanZero = pd.merge(fragFClean, fcFile[fcFile.old == 0], how = "inner", on = keyCoord).groupby(keyCoord).agg(np.random.choice).reset_index()
        fragFCleanZero["freq"] = fragFCleanZero["new"]
        fragFCleanZeroFrag = fragFCleanZero[["fragA", "fragB", "freq"]].groupby(["fragA", "fragB"]).sum().reset_index()
        fragFCleanZeroPosFrag = fragFCleanPosFrag.append(fragFCleanZeroFrag).groupby(["fragA", "fragB"]).sum().reset_index()

        ## merge with other fragment pairs that have no change
        inInd = pd.merge(interFreq.reset_index(), fragFCleanZeroPosFrag, how = "inner", on = ["fragA", "fragB"]).set_index('index')
        fragFinal = interFreq[~interFreq.index.isin(inInd.index)].append(fragFCleanZeroPosFrag)
        fragFinal

        fragFinal.to_csv(outPath + "/" + interFileName + ".spikeIn", sep = "\t", header = False, index = False)
        print("Spike-in for " + interFileName + " has been generated!")
