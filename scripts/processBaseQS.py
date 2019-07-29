#!/usr/bin/env python
import sys

if __name__ == "__main__":
    # file = "/home/yzheng74/freeHiC/jobs_A549/rep2/rep2.training.summary.baseQS.raw"
    file = sys.argv[1] 
    
    baseDic={}
    with open(file) as baseQSFile:
        for line in baseQSFile:
            if ':' in line:
                
                key, val = line.split(':', 1)
                valL = len(val)
                valList = val[1:(valL-2)].split(', ')
                valListN = [int(x) for x in valList]
                if key not in baseDic.keys():
                    baseDic[key] = valListN
                else:
                    baseDic[key] = [sum(x) for x in zip(baseDic[key], valListN)]
    
        
            else:
                pass
    outFile = open(file[:-4], "w+")
    for k in range(1, len(baseDic.keys())+1):
        outFile.write(str(k) + ":" + str(baseDic[str(k)]) + "\n")
        
    outFile.close()
        
