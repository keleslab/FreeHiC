cimport cython
import numpy as np
from subprocess import call
from libc.stdlib cimport qsort
from cpython cimport array
from libc.stdlib cimport malloc
from libc.stdlib cimport rand, RAND_MAX


cdef char *reverse_seq(char *seq_src, int seq_len):

    cdef char *seq_dest = <char *>malloc(seq_len)
    #cdef bytes py_bytes = seq.encode('ascii')
    #cdef char *seq_src = py_bytes
    cdef int i = 0
    cdef int d = 0
    
    # seq_dest[seq_len] = '\0'
    for i in xrange(seq_len):
        d = seq_len - i - 1
        if seq_src[i] == 65:
            seq_dest[d] = 84
        elif seq_src[i] == 71:
            seq_dest[d] = 67
        elif seq_src[i] == 67:
            seq_dest[d] = 71
        elif seq_src[i] == 84:
            seq_dest[d] = 65
    return seq_dest#[:seq_len].decode('ascii')

cdef void sort_c(int *a, int SIZE):
    cdef int i, passed, hold

       # ------do the sorting...ascending-------------
       # for every array elements do this...
    for passed in range(SIZE): #(passed = 1; passed <= (SIZE-1); passed++):
        # for every 2 array elements comparison do
        # the comparison and swap...
        for i in xrange(SIZE -1): #(i = 0; i <= (SIZE-2); i++):# set the condition...
            if a[i] > a[i + 1]:# put the a[i] in temporary variable hold...
                hold = a[i]
                     # put the a[i + 1] in a[i]
                a[i] = a[i + 1]
                     # put the hold in a[i + 1], one swapping is
                     #completed...and repeat for other elements...
                a[i + 1] = hold
            elif a[i] == a[i+1]:
                a[i+1] += 1   



def add_mismatch(str seq, long revStrand, str mutation):

    cdef bytes seq_bytes = seq.encode('ascii')
    cdef char *seq_src = seq_bytes
    cdef bytes mut_bytes = mutation.encode('ascii')
    cdef char *mut_src = mut_bytes
    cdef int seqN = len(seq), i, j, flag = 0
    cdef long mutationN = len(mutation), ind
    cdef int *mutIndex= <int *>malloc(mutationN)
    
    if not mutIndex:
        raise MemoryError()

    if revStrand:
        mut_src_tmp = reverse_seq(mut_src, mutationN)
        mut_src = mut_src_tmp
    
    for j in xrange(100):
        mutI = np.random.choice(seqN, size = 3, replace = False)
        for i in xrange(mutationN):
            mutIndex[i] = mutI[i] #rand() % seqN #int((rand()/RAND_MAX)*seqN) # random integer #mutIndex = np.random.choice(seqN, size = mutationN)
        sort_c(mutIndex, mutationN)
        flag = 0
        for i in xrange(mutationN):
            if seq_src[mutIndex[i]] == mut_src[i]:
                flag += 1
        if flag == 0:
            break
            
    for i in xrange(mutationN):
        seq_src[mutIndex[i]] = mut_src[i]


    return seq_src[:seqN].decode('ascii')


@cython.profile(True)
def add_insert(str seq, long insertionN):
    
    cdef bytes seq_bytes = seq.encode('ascii')
    cdef char *seq_src = seq_bytes
    cdef int seqN = len(seq),  i, c, base, tmpBase, baseInd
    cdef char *seq_new = <char *>malloc(seqN + insertionN)
    cdef int *coord= <int *>malloc(insertionN)

    if not coord:
        raise MemoryError()
    
    for i in xrange(seqN):
        tmpBase = seq_src[i]
        seq_new[i] = tmpBase
    
    for i in xrange(insertionN):
        coord[i] = rand() % seqN #int((rand()/RAND_MAX)*seqN)
    sort_c(coord, insertionN)

    for i in xrange(insertionN):
        baseInd = rand() % 4
        if baseInd == 0:
            base = 65 #A
        elif baseInd == 1:
            base = 84 #T
        elif baseInd == 2:
            base = 67 #C
        else:
            base = 71 #G
        for c in xrange(seqN + i - coord[i]):
            tmpBase = seq_new[seqN + i - c -1]
            seq_new[seqN + i -c] = tmpBase
            
        seq_new[coord[i]] = base
    
        
    return seq_new[:(seqN + insertionN)].decode('ascii')



@cython.profile(True)
def add_delete(str seq, long deletionN):
    cdef bytes seq_bytes = seq.encode('ascii')
    cdef char *seq_src = seq_bytes
    cdef int seqN = len(seq), i, c, base, tmpBase, baseInd
    cdef int *coord= <int *>malloc(deletionN)
    
    if not coord:
        raise MemoryError()

    for i in xrange(deletionN):
        coord[i] = rand() % seqN
    sort_c(coord, deletionN)

    for i in xrange(deletionN):
        for c in range(coord[i]+1, seqN):
            tmpBase = seq_src[c]
            seq_src[c-1] = tmpBase
        
    return seq_src[:(seqN - deletionN)].decode('ascii')

@cython.profile(True)
def read_totalN(str summaryFile, long allMatch, long revStrandN, long chimericN):

    cdef str line
    
    ## read in training parameters files
    with open(summaryFile + ".allMatch") as allMatchFile:
        for line in allMatchFile:
            allMatch = int(line)
    with open(summaryFile + ".revStrandN") as revStrandFile:
        for line in revStrandFile:
            revStrandN = int(line)

    with open(summaryFile + ".chimericN") as chimericFile:
        for line in chimericFile:
            chimericN = int(line)

    return (allMatch, revStrandN, chimericN)

cdef struct fragStruct:
        int chrom
        int start
        int end
        int fragLen



        
@cython.cdivision(True)
cdef void simulation(str fragment, str interaction, int number, str summaryFile, long allMatch, long revStrandN, long chimericN, int mutationP, int indelP, int distance, int readLen, str outdir, str bedtools, str fileName, str genome):
    
    ## read in training data
    cdef list freq = [], frag1 = [], frag2 = [], mismatchValues = [], mismatchKeys = [], insertValues = [], insertKeys = [], deleteValues = [], deleteKeys = []
    cdef list interInfo, mismatchProb, insertProb, deleteProb, valList, fragProp, qprob

    cdef char ** baseQualSeq1 = <char **>malloc(number * sizeof(char*))
    cdef char ** baseQualSeq2 = <char **>malloc(number * sizeof(char*))
    cdef int *readLenList = <int *>malloc(number * sizeof(int))
    if not baseQualSeq1 or not baseQualSeq2 or not readLenList:
        raise MemoryError()

    cdef dict resFragIndex = {}
    cdef list resFragChrom = [], resFragStart = [], resFragEnd = []
    cdef long lineN = 0, contactN = 0, interactionN = 0, mismatchIndex, insertIndex, deleteIndex, gapLen
    cdef long mismatchN = 0, insertN = 0, insertTmpN = 0, deleteTmpN = 0, deleteN = 0, mismatchLen = 0, insertLen = 0, deleteLen = 0, tmpSum = 0
    cdef long tmpN, p, pos, qLen, qSum, i, end, keyInt
    cdef int chrom, fCount = 0, qsLen = 42, baseQSsum = 0, readLenCount = 0, rLen
    cdef long[:] qscore, index
    cdef str line, seq, key, value, baseQualSeqStr, readPos

    print("Read in restriction fragment dictionary!")
    ## construct restriction fragment dictionary
    with open(fragment) as fragFile:
        for line in fragFile:            
            fragInfo = line.split()
            resFragIndex[fragInfo[3]] = fCount
            resFragChrom.append(fragInfo[0][3:])
            resFragStart.append(fragInfo[1])
            resFragEnd.append(fragInfo[2])
            fCount += 1

    print("Read in interaction frequency!")
    ## read in fragment interactions
    with open(interaction, "r") as interFile:
        for line in interFile:
            interInfo = line.split()
            tmpN = int(interInfo[0])
            freq.append(tmpN)
            contactN = contactN + tmpN
            fragInterIndex = resFragIndex[interInfo[1]]
            frag1.append(fragInterIndex)
            fragInterIndex = resFragIndex[interInfo[2]]
            frag2.append(fragInterIndex)
            interactionN += 1
    ## refresh number
    if number == 0:
        number = contactN

    print("Convert dictionary and list into C structure!")
    ## convert restrction fragment dictionary into C structure
    cdef fragStruct *resFragC = <fragStruct *> malloc(fCount*sizeof(fragStruct))
    cdef fragStruct resFragTmp
    
    if not resFragC:
        raise MemoryError()

    for i in xrange(fCount):        
        resFragC[i].chrom = int(resFragChrom[i])
        resFragC[i].start = int(resFragStart[i])
        resFragC[i].end = int(resFragEnd[i])
        resFragC[i].fragLen = resFragC[i].end - resFragC[i].start
    

    ## convert fragment interaction and sample index into C structure
    cdef int *frag1C = <int *>malloc(interactionN * cython.sizeof(int))
    cdef int *frag2C = <int *>malloc(interactionN * cython.sizeof(int))
    cdef long *indexC = <long *>malloc(number * cython.sizeof(long))
    if not frag1C or not frag2C or not indexC:
        raise MemoryError()


    for i in xrange(interactionN):
        frag1C[i] = frag1[i]
        frag2C[i] = frag2[i]
        

    print("Sample fragment interaction!")
    ## sample fragment interaction index
    fragProp = list(np.array(freq)/contactN)
    index = np.random.choice(interactionN, size = number, replace = True, p = list(fragProp))
    for i in xrange(number):
        indexC[i] = index[i]

    
    ## sample restriction fragment start/end part, sequence strand, distrance
    cdef int *seg1 = <int *>malloc(number*cython.sizeof(int))
    cdef int *seg2 = <int *>malloc(number*cython.sizeof(int))
    cdef int *dist1 = <int *>malloc(number*cython.sizeof(int))
    cdef int *dist2 = <int *>malloc(number*cython.sizeof(int))
    cdef int *revStrand1 = <int *>malloc(number*cython.sizeof(int))
    cdef int *revStrand2 = <int *>malloc(number*cython.sizeof(int))
    cdef int *baseQS = <int *>malloc(qsLen*cython.sizeof(int))
    if not seg1 or not seg2 or not dist1 or not dist2 or not revStrand1 or not revStrand2 or not baseQS:
        raise MemoryError()

    

    for i in xrange(number):
        seg1[i] = rand() % 2
        seg2[i] = rand() % 2
        dist1[i] = rand() % distance
        dist2[i] = rand() % distance
        if np.random.choice(allMatch, size = 1)[0] < revStrandN: #int(rand()/RAND_MAX*allMatch) < revStrandN: #rand() % allMatch < revStrandN:
            revStrand1[i] = 1
        else:
            revStrand1[i] = 0

        if np.random.choice(allMatch, size = 1)[0] < revStrandN: #int(rand()/RAND_MAX*allMatch) < revStrandN: #rand() % allMatch < revStrandN:
            revStrand2[i] = 1
        else:
            revStrand2[i] = 0
    
    print("Read in trained parameters!")
    with open(summaryFile + ".mismatch") as misFile:
        for line in misFile:
            if ':' in line:
                key, value = line.split(':', 1)
                mismatchValues.append(int(value))
                mismatchKeys.append(key)
                mismatchN = mismatchN + int(value)
                mismatchLen = mismatchLen + 1
            else:
                pass
        


    with open(summaryFile + ".insertion") as misFile:
        for line in misFile:
            if ':' in line:
                key, value = line.split(':', 1)
                insertValues.append(int(value))
                insertKeys.append(int(key))
                insertN = insertN + int(value)
                insertLen = insertLen + 1
            else:
                pass

    with open(summaryFile + ".deletion") as misFile:
        for line in misFile:
            if ':' in line:
                key, value = line.split(':', 1)
                deleteValues.append(int(value))
                deleteKeys.append(int(key))
                deleteN = deleteN + int(value)
                deleteLen = deleteLen + 1
            else:
                pass

    for i in range(number):
        baseQualSeq1[i] = <char *>malloc(readLen * sizeof(char))
        baseQualSeq2[i] = <char *>malloc(readLen * sizeof(char))

    with open(summaryFile + ".baseQS") as baseQSFile:
        for line in baseQSFile:
            if ':' in line:
                if readLenCount >= readLen:
                    break
                else:

                    key, val = line.split(':', 1)
                    valL = len(val)
                    valList = val[1:(valL-2)].split(', ')
                    baseQSsum = 0
                    for v in xrange(qsLen):
                        baseQS[v] = int(valList[v])
                        baseQSsum += baseQS[v]
        
                    qprob = [float(baseQS[v])/baseQSsum for v in range(qsLen)]
                    qscore = np.random.choice(qsLen, size = number, p = qprob)
                    for i in range(number):
                        baseQualSeq1[i][readLenCount] = qscore[i] + 33
                    qscore = np.random.choice(qsLen, size = number, p = qprob)
                    for i in range(number):
                        baseQualSeq2[i][readLenCount] = qscore[i] + 33
                    readLenCount += 1
            else:
                pass
    

    ## sample interaction coordinate and strand
    mismatchProb = [mis/mismatchN for mis in mismatchValues]
    insertProb = [ins/insertN for ins in insertValues]
    deleteProb = [dele/deleteN for dele in deleteValues]

    cdef long randNum
    cdef str strand="+-"
    cdef int rid

    if indelP != 404:
        insertTmpN = int(insertN * indelP * allMatch / (100.0 * (insertN + deleteN)))
        deleteTmpN = int(deleteN * indelP * allMatch / (100.0 * (insertN + deleteN)))
        insertN = insertTmpN
        deleteN = deleteTmpN

    if mutationP != 404:
        mismatchN = int(mutationP * allMatch / 100.0)

    print(allMatch)
    print(mismatchN)
    print(insertN)
    print(deleteN)
# 71891702
# 5016295
# 301627
# 196755

    print("Simulate end 1!")
    ## for end 1
    readPosOri1 = open(outdir + "/" + fileName + ".readPosOri1", "w+")
    readPosMut1 = open(outdir + "/" + fileName + ".readPosMut1", "w+")
    readPosIns1 = open(outdir + "/" + fileName + ".readPosIns1", "w+")
    readPosDel1 = open(outdir + "/" + fileName + ".readPosDel1", "w+")

    readqsOri1 = open(outdir + "/" + fileName + ".readqsOri1", "w+")
    # readqsMut1 = open(outdir + "/" + fileName + ".readqsMut1", "w+")

    for i in xrange(number):
        if i % 500000 == 0:
            print(i)
        fi = frag1C[indexC[i]]
        resFragTmp = resFragC[fi]
        chrom = resFragTmp.chrom
        if np.random.choice(allMatch, size = 1)[0] < chimericN:
            rLen = np.random.choice(range(25, readLen), size = 1)[0]
            pos = resFragTmp.end - revStrand1[i] * resFragTmp.fragLen + (revStrand1[i] - 1)*rLen
            end = resFragTmp.end - revStrand1[i] * resFragTmp.fragLen + revStrand1[i] * rLen
            
        else:
            rLen = readLen
            pos = resFragC[fi].start + dist1[i] + seg1[i] * (resFragC[fi].fragLen - distance - readLen) ## simulate reads at the beginning or end of RE + reads direction
            end = pos + readLen

        readLenList[i] = rLen
        randNum = np.random.choice(allMatch, size = 1)[0] #int(rand()/RAND_MAX*allMatch)  #rand() % allMatch
        if randNum < mismatchN: ## mutation: including mutation only, mutation + insertion, mutation + deletion
            readPosMut1.write('chr' + str(chrom) + '\t' + str(pos) + '\t' + str(end) + "\t" + str(i) + "\t100\t" + strand[revStrand1[i]] + "\n")
        elif randNum < (insertN + mismatchN): ## insertion
            readPosIns1.write('chr' + str(chrom) + '\t' + str(pos) + '\t' + str(end) + "\t" + str(i) + "\t100\t" + strand[revStrand1[i]] + "\n")
        elif randNum < (insertN + deleteN + mismatchN): ## deletion
            readPosDel1.write('chr' + str(chrom) + '\t' + str(pos) + '\t' + str(end) + "\t" + str(i) + "\t100\t" + strand[revStrand1[i]] + "\n")
        else: ## original
            readPosOri1.write('chr' + str(chrom) + '\t' + str(pos) + '\t' + str(end) + "\t" + str(i) + "\t100\t" + strand[revStrand1[i]] + "\n")
            readqsOri1.write('+\n' + baseQualSeq1[i][:rLen].decode('ascii') + '\n')


    readPosOri1.close()
    readPosMut1.close()
    readPosIns1.close()
    readPosDel1.close()
    readqsOri1.close()
    # readqsMut1.close()

    print("Get fasta from the reference genome!")
    call(bedtools + ' getfasta -fi ' + genome + ' -bed ' + outdir + '/' + fileName + '.readPosOri1 -name -s -fo ' + outdir + "/" + fileName + '.readSeqOri1', shell = True)
    call(bedtools + ' getfasta -fi ' + genome + ' -bed ' + outdir + '/' + fileName + '.readPosMut1 -name -s -tab -fo ' + outdir + "/" + fileName + '.readSeqMut1', shell = True)
    call(bedtools + ' getfasta -fi ' + genome + ' -bed ' + outdir + '/' + fileName + '.readPosIns1 -s -name -tab -fo ' + outdir + "/" + fileName + '.readSeqIns1', shell = True)
    call(bedtools + ' getfasta -fi ' + genome + ' -bed ' + outdir + '/' + fileName + '.readPosDel1 -s -name -tab -fo ' + outdir + "/" + fileName + '.readSeqDel1', shell = True)
    
    call(r"""awk '{{ split($0, id, "("); name = id[1]; getline; getline x<"{0}/{1}.readqsOri1"; getline y<"{0}/{1}.readqsOri1"; print name "\n" $0 "\n" x "\n" y;}}' {0}/{1}.readSeqOri1 > {0}/{1}_1.fastq""".format(outdir, fileName), shell = True)
    endFile1 = open(outdir + "/" + fileName + "_1.fastq", "a+")

    print("Add mutation!")
    ## add mutation
    cdef long mutCount = 1
    with open(outdir + "/" + fileName + '.readSeqMut1') as mutFile:
        for line in mutFile:
            if mutCount % 500000 == 0:
                print(mutCount)
            mutCount += 1
            ridTmp, seq = line.strip("\n").split("\t")
            try:
                rid = int(ridTmp)
            except:
                rid = int(ridTmp[0:-3])
                          
            mismatchIndex = np.random.choice(mismatchLen, size = 1, replace = True, p = mismatchProb)[0]
            seq = add_mismatch(seq, revStrand1[rid],  mismatchKeys[mismatchIndex])

            ##add gap
            randNum = np.random.choice(allMatch, size = 1)[0] #int(rand()/RAND_MAX*allMatch)  #rand() % allMatch
            if randNum < insertN: #insert

                insertIndex = np.random.choice(insertLen, size = 1, replace = True, p = insertProb)[0]
                gapLen = insertKeys[insertIndex]
                seq = add_insert(seq, gapLen)
                baseQualSeqStr = baseQualSeq1[rid][:readLenList[rid]].decode('ascii') + baseQualSeq1[rid][(readLenList[rid]-2):(readLenList[rid] - 1)].decode('ascii') * gapLen

            elif randNum < (insertN + deleteN): # delete
                deleteIndex = np.random.choice(deleteLen, size = 1, replace = True, p = deleteProb)[0]
                gapLen = deleteKeys[deleteIndex]
                seq = add_delete(seq, gapLen)
                baseQualSeqStr = baseQualSeq1[rid][:(readLenList[rid] - gapLen)].decode('ascii')
            else:
                baseQualSeqStr = baseQualSeq1[rid][:readLenList[rid]].decode('ascii')

            endFile1.write("@" + str(rid) + "\n" + seq + "\n+\n" + baseQualSeqStr + "\n")

    print("Add insertion!")
    ## add insertion
    with open(outdir + "/" + fileName + '.readSeqIns1') as insFile:
        for line in insFile:
            ridTmp, seq = line.strip("\n").split("\t")
            try:
                rid = int(ridTmp)
            except:
                rid = int(ridTmp[0:-3])
            insertIndex = np.random.choice(insertLen, size = 1, replace = True, p = insertProb)[0]
            gapLen = insertKeys[insertIndex]
            seq = add_insert(seq, gapLen)
            baseQualSeqStr = baseQualSeq1[rid][:readLenList[rid]].decode('ascii') + baseQualSeq1[rid][(readLenList[rid]-2):(readLenList[rid] - 1)].decode('ascii') * gapLen

            endFile1.write("@" + str(rid) + "\n" + seq + "\n+\n" + baseQualSeqStr + "\n")

    print("Add deletion")
    ## add deletion
    with open(outdir + "/" + fileName + '.readSeqDel1') as insFile:
        for line in insFile:
            ridTmp, seq = line.strip("\n").split("\t")
            try:
                rid = int(ridTmp)
            except:
                rid = int(ridTmp[0:-3])
            deleteIndex = np.random.choice(deleteLen, size = 1, replace = True, p = deleteProb)[0]
            gapLen = deleteKeys[deleteIndex]
            seq = add_delete(seq, gapLen)
            baseQualSeqStr = baseQualSeq1[rid][:(readLenList[rid] - gapLen)].decode('ascii')

            endFile1.write("@" + str(rid) + "\n" + seq + "\n+\n" + baseQualSeqStr + "\n")
    endFile1.close()

    ## for end 2
    print("Simulate end 2!")
    readPosOri2 = open(outdir + "/" + fileName + ".readPosOri2", "w+")
    readPosMut2 = open(outdir + "/" + fileName + ".readPosMut2", "w+")
    readPosIns2 = open(outdir + "/" + fileName + ".readPosIns2", "w+")
    readPosDel2 = open(outdir + "/" + fileName + ".readPosDel2", "w+")

    readqsOri2 = open(outdir + "/" + fileName + ".readqsOri2", "w+")
    # readqsMut2 = open(outdir + "/" + fileName + ".readqsMut2", "w+")
    for i in xrange(number):
        if i % 500000 == 0:
            print(i)
        fi = frag2C[indexC[i]]
        resFragTmp = resFragC[fi]
        chrom = resFragTmp.chrom
        if np.random.choice(allMatch, size = 1)[0] < chimericN:
            rLen = np.random.choice(range(25, readLen), size = 1)[0]
            pos = resFragTmp.end - revStrand2[i] * resFragTmp.fragLen + (revStrand2[i] - 1)*rLen
            end = resFragTmp.end - revStrand2[i] * resFragTmp.fragLen + revStrand2[i] * rLen
            
        else:
            rLen = readLen
            pos = resFragC[fi].start + dist2[i] + seg2[i] * (resFragC[fi].fragLen - distance - readLen) ## simulate reads at the beginning or end of RE + reads direction
            end = pos + readLen

        readLenList[i] = rLen
        randNum = np.random.choice(allMatch, size = 1)[0] #int(rand()/RAND_MAX*allMatch)  #rand() % allMatch
        if randNum < mismatchN: ## mutation: including mutation only, mutation + insertion, mutation + deletion
            readPosMut2.write('chr' + str(chrom) + '\t' + str(pos) + '\t' + str(end) + "\t" + str(i) + "\t100\t" + strand[revStrand2[i]] + "\n")
        elif randNum < (insertN + mismatchN): ## insertion
            readPosIns2.write('chr' + str(chrom) + '\t' + str(pos) + '\t' + str(end) + "\t" + str(i) + "\t100\t" + strand[revStrand2[i]] + "\n")
        elif randNum < (insertN + deleteN + mismatchN): ## deletion
            readPosDel2.write('chr' + str(chrom) + '\t' + str(pos) + '\t' + str(end) + "\t" + str(i) + "\t100\t" + strand[revStrand2[i]] + "\n")
        else: ## original
            readPosOri2.write('chr' + str(chrom) + '\t' + str(pos) + '\t' + str(end) + "\t" + str(i) + "\t100\t" + strand[revStrand2[i]] + "\n")
            readqsOri2.write('+\n' + baseQualSeq2[i][:rLen].decode('ascii') + '\n')

    readPosOri2.close()
    readPosMut2.close()
    readPosIns2.close()
    readPosDel2.close()
    readqsOri2.close()
    # readqsMut2.close()

    print("Get fasta from reference genome!")
    call(bedtools + ' getfasta -fi ' + genome + ' -bed ' + outdir + '/' + fileName + '.readPosOri2 -name -s -fo ' + outdir + "/" + fileName + '.readSeqOri2', shell = True)
    call(bedtools + ' getfasta -fi ' + genome + ' -bed ' + outdir + '/' + fileName + '.readPosMut2 -name -s -tab -fo ' + outdir + "/" + fileName + '.readSeqMut2', shell = True)
    call(bedtools + ' getfasta -fi ' + genome + ' -bed ' + outdir + '/' + fileName + '.readPosIns2 -s -name -tab -fo ' + outdir + "/" + fileName + '.readSeqIns2', shell = True)
    call(bedtools + ' getfasta -fi ' + genome + ' -bed ' + outdir + '/' + fileName + '.readPosDel2 -s -name -tab -fo ' + outdir + "/" + fileName + '.readSeqDel2', shell = True)
    #call(r"""awk '{{ if(substr($0, length($0)-2, length($0)) == ")") {{ name = substr($0, 2, length($0)-3)}} else {{name = $0}}; getline; getline x<"{0}/{2}.readqsOri2"; getline y<"{0}/{2}.readqsOri2"; print name "\n" $0 "\n" x "\n" y;}}' {0}/{2}.readSeqOri2 > {0}/{2}_2.fastq""".format(outdir, fileName), shell = True)
    call(r"""awk '{{ split($0, id, "("); name = id[1]; getline; getline x<"{0}/{1}.readqsOri2"; getline y<"{0}/{1}.readqsOri2"; print name "\n" $0 "\n" x "\n" y;}}' {0}/{1}.readSeqOri2 > {0}/{1}_2.fastq""".format(outdir, fileName), shell = True)
    endFile2 = open(outdir + "/" + fileName + "_2.fastq", "a+")

    print("Add mutation!")
    ## add mutation
    mutCount = 1
    with open(outdir + "/" + fileName + '.readSeqMut2') as mutFile:
        for line in mutFile:
            if mutCount % 500000 == 0:
                print(mutCount)
            mutCount += 1
            ridTmp, seq = line.strip("\n").split("\t")
            try:
                rid = int(ridTmp)
            except:
                rid = int(ridTmp[0:-3])
                          
            mismatchIndex = np.random.choice(mismatchLen, size = 1, replace = True, p = mismatchProb)[0]
            seq = add_mismatch(seq, revStrand2[rid],  mismatchKeys[mismatchIndex])
            ## add gap
            randNum = np.random.choice(allMatch, size = 1)[0] #int(rand()/RAND_MAX*allMatch)  #rand() % allMatch
            if randNum < insertN: #insert

                insertIndex = np.random.choice(insertLen, size = 1, replace = True, p = insertProb)[0]
                gapLen = insertKeys[insertIndex]
                seq = add_insert(seq, gapLen)
                baseQualSeqStr = baseQualSeq2[rid][:readLenList[rid]].decode('ascii') + baseQualSeq2[rid][(readLenList[rid]-2):(readLenList[rid] - 1)].decode('ascii') * gapLen

            elif randNum < (insertN + deleteN): # delete
                deleteIndex = np.random.choice(deleteLen, size = 1, replace = True, p = deleteProb)[0]
                gapLen = deleteKeys[deleteIndex]
                seq = add_delete(seq, gapLen)
                baseQualSeqStr = baseQualSeq2[rid][:(readLenList[rid] - gapLen)].decode('ascii')
            else:
                baseQualSeqStr = baseQualSeq2[rid][:readLenList[rid]].decode('ascii')
            ## baseQualSeqStr = baseQualSeq2[rid][:readLenList[rid]].decode('ascii')
            endFile2.write("@" + str(rid) + "\n" + seq + "\n+\n" + baseQualSeqStr + "\n")

    print("Add insertion!")
    ## add insertion
    with open(outdir + "/" + fileName + '.readSeqIns2') as insFile:
        for line in insFile:
            ridTmp, seq = line.strip("\n").split("\t")
            try:
                rid = int(ridTmp)
            except:
                rid = int(ridTmp[0:-3])
            insertIndex = np.random.choice(insertLen, size = 1, replace = True, p = insertProb)[0]
            gapLen = insertKeys[insertIndex]
            seq = add_insert(seq, gapLen)
            baseQualSeqStr = baseQualSeq2[rid][:readLenList[rid]].decode('ascii') + baseQualSeq2[rid][(readLenList[rid]-2):(readLenList[rid] - 1)].decode('ascii') * gapLen

            endFile2.write("@" + str(rid) + "\n" + seq + "\n+\n" + baseQualSeqStr + "\n")

    print("Add deletion!")
    ## add deletion
    with open(outdir + "/" + fileName + '.readSeqDel2') as insFile:
        for line in insFile:
            ridTmp, seq = line.strip("\n").split("\t")
            try:
                rid = int(ridTmp)
            except:
                rid = int(ridTmp[0:-3])
            deleteIndex = np.random.choice(deleteLen, size = 1, replace = True, p = deleteProb)[0]
            gapLen = deleteKeys[deleteIndex]
            seq = add_delete(seq, gapLen)
            baseQualSeqStr = baseQualSeq2[rid][:(readLenList[rid] - gapLen)].decode('ascii')

            endFile2.write("@" + str(rid) + "\n" + seq + "\n+\n" + baseQualSeqStr + "\n")
    endFile2.close()

    


def main(str fragment, interaction, outdir, fileName, number, mutationP, indelP, chimericP, readLen, distance, bedtools, genome, summaryFile, verbose):

    cdef long allMatch = 0, revStrandN = 0, chimericN = 0
    (allMatch, revStrandN, chimericN) = read_totalN(summaryFile, allMatch, revStrandN, chimericN)
    if chimericP != 404:
        chimericN = int(chimericP * allMatch /100.0)
    simulation(fragment, interaction, number, summaryFile, allMatch,  revStrandN, chimericN, mutationP, indelP, distance, readLen, outdir, bedtools, fileName, genome)

