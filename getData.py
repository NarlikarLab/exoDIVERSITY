import numpy as np
import re,itertools,os
from config import *
#import pprint

def printHelp():
    '''
    print the whole process to provide input
    '''
    print 'exoDiversity [options]'
    print '\t -f: Input fasta file'
    print '\t -r: Reads file either in BAM (sorted) or bedGraph format'
    print '\t -format: Format of the reads file BAM or BED'
    print '\t -g: genome file containing sizes of call chromosomes'
    print '\t -ctrl: Control reads file in the same format as the reads file'
    print '\t -o: Output directory'
    print '\t -rev: is 1 if reverse complement is to be considered; otherwise 0. Default 1'
    print '\t -mask: is 1 if repeats are to be masked; otherwise 0. Default 0'
    print '\t -initialWidth: The width of the motifs at starting point'
    print '\t -minMode: Minimum number of modes in which data should be divided. Default 1'
    print '\t -maxMode: Maximum number of modes in which the data should be divided. Default 10'
    print '\t -rWidth: The width of the read windows for both positive and negative strand. Default 5'
    print '\t -nproc: The number of processors to be used for computation. Default is the number of cores the system has'
    print '\t -bin: Binarize read counts based on median, first quartile or third quartile or keep when file is already in binary form {median,Q1,Q3,keep}. Default median'
    print '\t -ntrials: Number of trials for each model. Default 5'
    print '\t -pcZeros: Pseudo count for 0s in reads data. Default 0.9'
    print '\t -pcOnes: Pseudo count for 1s in reads data. Default 0.1'
    print '\t -twobit: 2bit file (from UCSC browser) for sequence alignment plots'
    print 'In case sequence wise + and - ve read counts are present'
    print '\t -p: The positive strand reads file'
    print '\t -n: The negative strand reads file'

def checkValidFile(s):
    if s == '':
        print 'One or more input files are missing'
        printHelp()
        exit()
    if not os.path.isfile(s):
        print 'Could not open file' + s
        exit()
    else:
        return s

def checkFormat(s):
    if s=='BED' or s=='BAM':
        return s
    else:
        print 'Reads file format',s,'is not supported. \n Please enter BAM (sorted) or BED format reads.'
        printHelp()

def getOutdir(s):
    if not(os.path.isdir(s)):
        os.system('mkdir -p '+s)
    else:
        print 'Directory already exists, overwrite contents? (y/n) '
        inp = raw_input()
        if inp == 'y':
            if len(os.listdir(s))==0:
                pass
            else:
                os.system(' rm -r '+s+'/*')
        else:
            print 'Please enter a different directory name'
            exit()
    return s

def getFlag(s,opt):
    s1 = re.search('(0|1)\Z',s)
    if s1 is None:
        print 'Invalid value for ' + opt + ' :' + s +'.Please enter 0 or 1 as flag value.'
        printHelp()
        exit()
    s = s1
    s = int(s.group(0))
    return s

def getValidPosNum(s,opt):
    s1 = re.search('^[1-9][0-9]*\Z',s)
    if s1 is None:
        print 'Invalid value for' +opt +': ' + s + '.Please enter a valid positive number.'
        printHelp()
        exit()
    s = s1
    s = int(s.group(0))
    return s

def getValidRealNum(s,opt):
    s1 = re.search('^[0-9]+?(\.[0-9]+)?\Z', s)
    if s1 is None:
        print 'Invalid value for '+opt+': '+s+'.Please enter a valid real number between 0 and 1'
        printHelp()
        exit()
    s = s1
    s = float(s.group(0))
    if s > 1 and opt not in ['-pcZeros','-pcOnes']: 
        printHelp()
        exit()
    return s

def getValidBinarizingParam(s,opt):
    if s in ["median","Q1","Q3","keep"]:
        return s
    else:
        print "Invalid input for binarizing parameter "+ opt + " : "+ s
        printHelp()
        exit()

def fetchParameters(options):
    # options are sys.argv
    options = options[1:]
    if len(options) < 1:
        printHelp()
        exit()

    d = {'-f':'','-r':'','-format':'','-g':'','-ctrl':'','-p':'','-n':'','-twobit':'','-rev':defaultRevFlag,'-o':'','-mask':defaultMaskFlag,'-initialWidth':defaultInitialWidth,'-minMode':defaultminMode,'-maxMode':defaultmaxMode,'-rWidth':defaultReadWidth,'-nproc':defaultNProc,'-bin':defaultBinarizingParam, '-ntrials':defaultnoOfModels, '-pcZeros':defaultPseudocountZeros, '-pcOnes':defaultPseudocountOnes,'-seed':defaultSeed}
    validOptions = d.keys()
    dF = {'-f':checkValidFile,'-twobit':checkValidFile,'-r':checkValidFile,'-g':checkValidFile,'-ctrl':checkValidFile,'-p':checkValidFile,'-n':checkValidFile,'-o':getOutdir,'-format':checkFormat,'-rev':getFlag,'-mask':getFlag,'-initialWidth':getValidPosNum,'-rWidth':getValidPosNum,'-minMode':getValidPosNum,'-maxMode':getValidPosNum,'-nproc':getValidPosNum,'-bin':getValidBinarizingParam,'-ntrials':getValidPosNum, '-pcZeros':getValidRealNum, '-pcOnes':getValidRealNum,'-seed':getValidPosNum}
    for i in range(len(options)):
        if (i%2 == 0):
            if (options[i] not in validOptions):
                print options[i] + ' is not a valid option'
                printHelp()
                exit()
        else:
            if options[i-1] in ['-f','-r','-g','-ctrl','-twobit','-p','-n','-o']:
                op = dF[options[i-1]](options[i])
            elif options[i-1] == '-format':
                op = dF[options[i-1]](options[i])
            else:
                op = dF[options[i-1]](options[i],options[i-1])
            d[options[i-1]] = op

    if d['-minMode'] > d['-maxMode']:
        print '-maxMode cannot be less than -minMode'
        exit()
    if d['-r']!='':
        if d['-format']=='':
            print "Format is empty!\n"
            exit()
        if d['-format']=='BED' and d['-g']=='':
            print "Please enter a valid genome file containign sizes of all chromosomes"
            exit()

    if d['-f']=='':
        print 'Input sequence fasta file is compulsory'
        exit()
    if d['-r']=='' and (d['-p']=='' or d['-n']==''):
        print 'Input reads file is compulsory or please provide sequence wise reads for + and - strands separately'
        exit()
    return d

def readSeqs(sequencefile,outdir,maskReps):
    outfile = outdir + '/headers.txt'
    sf = open(sequencefile,'r')
    of = open(outfile,'w')
    countseqs = 0
    for l in sf:
        if l[0] == '>':
            of.write(l[1:])
            continue
        elif (l == ''):
            continue
        else:
            countseqs += 1
    of.close()
    sf.seek(0)
    seqs = np.zeros(countseqs,dtype=object)
    ctr = 0
    seqflag = 0
    #maskflag = 0
    for l in sf:
        if l[0] == '>':
            seqflag = 1
            continue
        elif l == '':
            seqflag = 0
            continue
        elif seqflag == 1:
            s = l.strip('\n')
            if not(maskReps): # maskReps != 1
                sl_replaceA = re.sub(r'[aA]',r'0',s)
                sl_replaceC = re.sub(r'[cC]',r'1',sl_replaceA)
                sl_replaceG = re.sub(r'[gG]',r'2',sl_replaceC)
                sl_replaceT = re.sub(r'[tT]',r'3',sl_replaceG)
                if 'N' in s:
                    sl_replaceN = re.sub(r'[N]',r'4',sl_replaceT)
                    #maskflag = 1
                    sl_final = sl_replaceN
                else:
                    sl_final = sl_replaceT
            elif maskReps:
                sl_replaceA = re.sub(r'A',r'0',s)
                sl_replaceC = re.sub(r'C',r'1',sl_replaceA)
                sl_replaceG = re.sub(r'G',r'2',sl_replaceC)
                sl_replaceT = re.sub(r'T',r'3',sl_replaceG)
                sl_replaceRepeats = re.sub(r'[acgtN]',r'4',sl_replaceT)
                sl_final = sl_replaceRepeats
		#maskflag = 1
            sarr = np.array([int(ele) for ele in sl_final],dtype=int)
            seqs[ctr] = sarr
            ctr += 1
            seqflag = 0
    sf.close()
    return seqs


def getSeqWithRevcomp(sequences):
    complement = {0:3,1:2,2:1,3:0,4:4}
    size = np.shape(sequences)
    new_seqarr = np.zeros(size[0],dtype=object)
    for i in range(size[0]):
        new_seq = np.zeros(2*len(sequences[i]),dtype=int)
        j = 0
        for k in range(len(sequences[i])):
            new_seq[j] = sequences[i][k]
            j+=1
        for ele in reversed(sequences[i]):
            new_seq[j] = complement[ele]
            j+=1
        new_seqarr[i] = new_seq
    return new_seqarr

def appendRevReads(posreads,negreads):
    noOfseqs = len(posreads)
    posreadswithrc = np.zeros(noOfseqs,dtype=object)
    negreadswithrc = np.zeros(len(negreads),dtype=object)
    for i in range(noOfseqs):
        #for positive reads the rev com will be flipped neg reads
        posreadswithrc[i] = np.zeros((len(posreads[i]) + len(negreads[i])),dtype=int)
        posreadswithrc[i][:len(posreads[i])] = posreads[i]
        posreadswithrc[i][len(posreads[i]):] = np.flipud(negreads[i])
        
        negreadswithrc[i] = np.zeros((len(negreads[i]) + len(posreads[i])),dtype=int)
        negreadswithrc[i][:len(negreads[i])] = negreads[i]
        negreadswithrc[i][len(negreads[i]):]= np.flipud(posreads[i])

    return (posreadswithrc,negreadswithrc)

def getMarkovProb(seq,order,L):
    freq = dict()
    for o in range(order+1):
        tempdict=dict()
        for i in xrange(o,L):
            k=tuple(seq[(i-o):(i+1)])
            if 4 in k:
                continue
            if k not in tempdict:
                tempdict[k]=1
            else:
                tempdict[k]+=1
        freq.update(tempdict)
        del(tempdict)

    markovparams=dict()
    for o in range(order+1):
        allnthorderkeys = filter((lambda x: len(x)==(o+1)),freq.keys())
        if o == 0:
            denom=sum([freq[k2] for k2 in allnthorderkeys])
            for k1 in allnthorderkeys:
                markovparams[k1]=float(freq[k1])/denom
        else:
            lookback = o
            #print allnthorderkeys
            for k1 in allnthorderkeys:
                denom = freq[k1[:lookback]]
                markovparams[k1] = float(freq[k1])/denom
    del(freq)
    return markovparams


def getMarkovVal(s,ind,markovdict,order):
    if ind < order:
        val = tuple(s[:ind+1])
        i=ind
    else:
        val = tuple(s[ind-order:ind+1])
        i=order
    finalval = []
    while i >= 0:
        if val[i]==4:
            break
        else:
            finalval[:0]= [val[i]]
        i=i-1
    if finalval == []:
        #print ind,finalval,-1
        return -1
    else:
        #print ind,finalval,markovdict[tuple(finalval)]
        return markovdict[tuple(finalval)]
    

def getBackground(seqarr,order):
    ### sequence wise background calculation
    background = np.zeros(len(seqarr),dtype=object)
    for i in range(len(seqarr)):
        s=seqarr[i]
        L=len(s)
        background[i]=np.zeros(L,dtype='Float64')
        markovdict = getMarkovProb(s,order,L)
        
        #pp = pprint.PrettyPrinter(indent=4)
        #pp.pprint(markovdict)
        for j in xrange(L):
            val = getMarkovVal(s,j,markovdict,order)
            background[i][j]=val
        #print background[i]
        #exit()
        del(markovdict)

    return background

def storeSettings(paramVals):
    outfile = paramVals['-o']+'/settings.txt'
    of = open(outfile,'w')
    of.write('Input fasta file: '+paramVals['-f']+'\n')
    if paramVals['-r'] !='':
        of.write('Reads file: '+paramVals['-r']+'\n')
    if paramVals['-ctrl']!='':
        of.write('Control reads file: '+paramVals['-ctrl']+'\n')
    of.write('Format of reads file: '+paramVals['-format'])
    of.write('Positive reads file: '+paramVals['-p']+'\n')
    of.write('Negative reads file:'+paramVals['-n']+'\n')
    of.write('Output directory: '+paramVals['-o']+'\n')
    of.write('Pseudo count for 0s: '+str(paramVals['-pcZeros'])+' and that for 1s: '+str(paramVals['-pcOnes'])+'\n')
    of.write('Minimum number of modes: '+str(paramVals['-minMode'])+'\n')
    of.write('Maximum number of modes: '+str(paramVals['-maxMode'])+'\n')
    of.write('Reverse strand: '+str(paramVals['-rev'])+'\n')
    of.write('Mask repeats: '+str(paramVals['-mask'])+'\n')
    of.write('Initial width: '+str(paramVals['-initialWidth'])+'\n')
    of.write('Number of trials per model: '+str(paramVals['-ntrials'])+'\n')
    of.write('Binarize read counts based on: '+paramVals['-bin']+'\n')
    of.close()

def getReadFrequencies(readfreqfile,noOfSeqs):
    rf = open(readfreqfile,'r')
    rline = rf.readline().strip()
    readfreqs = np.zeros(noOfSeqs,dtype=object)
    i=0
    while rline:
        rlist = map(float,rline.split('\t')[3:])
        #rlist = map(float,rline.split('\t')[2:])
        for j in range(len(rlist)):
            if(not(int(rlist[j])==0 or int(rlist[j])==1)):
                print j,rlist[j]
                print "One of the reads file is not in binary format"
                exit()
        readfreqs[i] = np.array(rlist)
        rline = rf.readline().strip()
        i+=1
    rf.close()
    return readfreqs

def makeLookaheadTable(sequences):
    noOfSeqs = len(sequences)
    lookahead = np.zeros(noOfSeqs,dtype=object)
    for i in range(noOfSeqs):
        lookahead[i] = np.zeros(len(sequences[i]),dtype=int)
        j=len(sequences[i])-1
        count=1
        while j>=0:
            if sequences[i][j]==4:
                lookahead[i][j]=0
                count=1
                j-=1
            else:
                lookahead[i][j] = count
                count+=1
                j-=1
    return lookahead
  
def getReadsBackground(reads,pcZeros,pcOnes):
    noOfSeqs = len(reads)
    rb = np.zeros(noOfSeqs,dtype=object)
    for i in range(noOfSeqs):
        L = len(reads[i])
        rb[i] = np.zeros(L,dtype=float)
        totalreads = L
        totalOnes = sum(reads[i])
        totalZeros = totalreads - totalOnes
        totalOnes += pcOnes
        totalZeros += pcZeros
        totalreads += (pcZeros + pcOnes)
        prob = [float(totalZeros)/totalreads,float(totalOnes)/totalreads]
        for j in range(L):
            rb[i][j] = prob[int(reads[i][j])]
    del prob
    return rb
