import sys,ctypes,os
from cstructures import *
import numpy as np
import getData as gd
from config import *
import saveFiles as sv
import readsProcessing as rp
import postExecWrapper as pw
import multiprocessing as mp
import time

def getValues(paramVals):
    allseqs = gd.readSeqs(paramVals['-f'],paramVals['-o'],paramVals['-mask'])
        
    noOfSeqs = len(allseqs)
    readfreq_pos = gd.getReadFrequencies(paramVals['-p'],noOfSeqs)
    readfreq_neg = gd.getReadFrequencies(paramVals['-n'],noOfSeqs)
    if not(all(map((lambda x:x==noOfSeqs),[len(readfreq_pos),len(readfreq_neg)]))):
        print "The number of sequences do not match the number of read sequences"
        exit()
    if paramVals['-rev'] == 1:
        allseqs = gd.getSeqWithRevcomp(allseqs)
        (readfreq_pos,readfreq_neg) = gd.appendRevReads(readfreq_pos,readfreq_neg)


    background = gd.getBackground(allseqs,2)
    ##  When we were taking sequence wise background freqs
    #posreadsback = gd.getBackground(readfreq_pos,0)
    #negreadsback = gd.getBackground(readfreq_neg,0)
    ## When we are taking sequence wise background freqs with pseudocounts
    
    posreadsback = gd.getReadsBackground(readfreq_pos,paramVals['-pcZeros'],paramVals['-pcOnes'])
    negreadsback = gd.getReadsBackground(readfreq_neg,paramVals['-pcZeros'],paramVals['-pcOnes'])

    return [allseqs,readfreq_pos,readfreq_neg,background,posreadsback,negreadsback]

def startTraining(args):
    try:
        sequences = args[0]
        posreads = args[1]
        negreads = args[2]
        outdir = args[3]
        twobit = args[4]
        mode = args[6].value
        rev = args[19].value
        seed = args[9].value
        args = args[5:]
        to = datalib.trainData(*args)
        noOfSeqs = len(sequences)
        pyTo = getTrainOut(to,mode,noOfSeqs)
        sv.saveDetails(sequences,posreads,negreads,outdir,twobit,rev,pyTo)
        print "Mode "+str(mode)+" with seed "+str(seed)+" completed"
        x=datalib.freeTo(to)
    except (KeyboardInterrupt, SystemExit):
        exit(1)
        
def callTraindata_parallel(pyvars):    
    params = pyvars[6]
    data = makeC2darray(pyvars[0],ctypes.c_int)
    lookahead = gd.makeLookaheadTable(pyvars[0])
    features = getCvar(map(len,pyvars[0]),ctypes.c_int)
    posreads = makeC2darray(pyvars[1],ctypes.c_double)
    negreads = makeC2darray(pyvars[2],ctypes.c_double)
    
    lookahead = makeC2darray(lookahead,ctypes.c_int)
    featureValues = 4
    noOfSeqs = len(pyvars[0])
    ds = DataSet(data,features,posreads,negreads,lookahead,featureValues,noOfSeqs)
    alpha = ctypes.c_float(defaultAlpha)
    pcReads = getCvar([params['-pcZeros'],params['-pcOnes']],c_float)
    background = makeC2darray(pyvars[3],ctypes.c_double)
    posreadsback = makeC2darray(pyvars[4],ctypes.c_double)
    negreadsback = makeC2darray(pyvars[5],ctypes.c_double)
    minWidth = defaultminWidth
    maxWidth = defaultmaxWidth
    
    noOfProc = params['-nproc']
    if noOfProc < 0:
        noOfProc = mp.cpu_count()

    joined = 0
    started = 0

    rWidth = params['-rWidth']
    motWidths = [getCvar([params['-initialWidth']]*mode,c_int) for mode in range(params['-minMode'],params['-maxMode']+1)]
    ########## Start from ntrials different starting positions
    #seeds = [params['-seed']+x for x in range(params['-ntrials'])]
    allprocs = sum([[(x,y) for x in range(1,params['-ntrials']+1)] for y in range(params['-minMode'],params['-maxMode']+1)],[])
    noOfModels = len(allprocs)
    pc = min(noOfModels,noOfProc)
    procArr = [0]*pc
    p = [0]*pc
    try:
        while joined < noOfModels:
            time.sleep(0.5)
            for j in range(pc):
                if procArr[j]!=0:
                    if not p[j].is_alive():
                        p[j].join()
                        procArr[j]=0
                        joined = joined+1
                if procArr[j]==0 and started < noOfModels:
                    trial = '/run'+str(allprocs[started][0])
                    modedir = '/'+str(allprocs[started][1])+'modes'
                    outdir = params['-o']+modedir+trial
                    print "Starting mode ",allprocs[started][1]," with seed: ",allprocs[started][0]+params['-seed']-1
                    if params['-v'] == 1:
                        likefile = outdir + '/likelihood.txt'
                    else:
                        likefile= ''
                    bestmodelfile = outdir + '/bestModelParams.txt'
                    os.makedirs(outdir)
                    # allprocs[started][0] ===> seed
                    # allprocs[started][1] ===> mode
                    mWidth = motWidths[allprocs[started][1] - params['-minMode']]
                    seed = ctypes.c_int(allprocs[started][0]+params['-seed']-1)
                    pr = mp.Process(target=startTraining, args=([pyvars[0],pyvars[1],pyvars[2],outdir,params['-twobit'],ctypes.byref(ds),ctypes.c_int(allprocs[started][1]),alpha,pcReads,seed,background,posreadsback,negreadsback,mWidth,ctypes.c_int(rWidth),ctypes.c_int(minWidth),ctypes.c_int(maxWidth),ctypes.c_int(defaultPositiveOffset),ctypes.c_int(defaultNegativeOffset),ctypes.c_int(params['-rev']),ctypes.c_int(params['-gobeyond']),likefile,bestmodelfile],))
                    p[j]=pr
                    p[j].start()
                    procArr[j]=allprocs[joined][0]
                    started = started +1
    except(KeyboardInterrupt, SystemExit):
        for i in p:
            if type(i) is int: continue
            i.terminate()
            i.join()
        exit(1)

    bestmodel = pw.fetchBestModel(params['-o'],params['-minMode'],params['-maxMode'],params['-ntrials'])
    if(params['-twobit']!=''):
        sv.createRPlots(params['-o'],params['-minMode'],params['-maxMode'],params['-o']+"/reads/posreads", params['-o']+"/reads/negreads",params['-twobit'])

    ##### plot seq scores in prob space
    alldata = [pyvars[0],pyvars[1],pyvars[2],pyvars[3],pyvars[4],pyvars[5],lookahead]
    bmfile = params['-o']+'/'+bestmodel+'/bestModelParams.txt'
    headerfile = params['-o']+'/headers.txt'
    scoreOutfile = params['-o']+'/'+bestmodel+'/seqScoresPerMode_normalized.txt'
    sv.seqsInProbSpace(alldata,bmfile,scoreOutfile,params['-rev'],params['-mask'],headerfile,params['-o']+'/'+bestmodel)
    #####

    ##### save HTML output
    sv.createHTML(params['-o'],params['-minMode'],params['-maxMode'],params['-twobit'],params['-f'],bestmodel)
    #####
    
    ##### delete files
    if(os.path.isdir(params['-o']+"/reads/posreads") and os.path.isdir(params['-o']+"/reads/negreads")):
        os.system("rm -r "+params['-o']+"/reads/posreads "+params['-o']+"/reads/negreads ")
    if params['-ctrl']!='':
        os.system("rm -r "+params['-o']+"/reads/control_posreads "+params['-o']+"/reads/control_negreads")
    #####
    print "All models completed"

def savetxt(outfile,arr):
    of = open(outfile,'w')
    for i in range(len(arr)):
        of.write('\t'.join(map(str,arr[i]))+'\n')
    of.close()

def callTraindata(pyvars):    
    params = pyvars[6]
    data = makeC2darray(pyvars[0],ctypes.c_int)
    lookahead = gd.makeLookaheadTable(pyvars[0])
    features = getCvar(map(len,pyvars[0]),ctypes.c_int)
    posreads = makeC2darray(pyvars[1],ctypes.c_double)
    negreads = makeC2darray(pyvars[2],ctypes.c_double)
    
    lookahead = makeC2darray(lookahead,ctypes.c_int)
    featureValues = 4
    noOfSeqs = len(pyvars[0])
    ds = DataSet(data,features,posreads,negreads,lookahead,featureValues,noOfSeqs)
    alpha = ctypes.c_float(defaultAlpha)
    pcReads = getCvar([params['-pcZeros'],params['-pcOnes']],c_float)
    background = makeC2darray(pyvars[3],ctypes.c_double)
    posreadsback = makeC2darray(pyvars[4],ctypes.c_double)
    negreadsback = makeC2darray(pyvars[5],ctypes.c_double)
    mWidth = getCvar([params['-initialWidth']]*params['-minMode'],c_int)
    minWidth = defaultminWidth
    maxWidth = defaultmaxWidth

    for s in range(9,10):
        seed=s
        print 'Seed ',seed
        od = params['-o'] +'/'+str(params['-minMode'])+'modes/run'+str(seed)
        os.makedirs(od)
        likefilename = od +'/likelihood.txt'
        bestmodelparams = od + '/bestModelParams.txt'
        savetxt(od+'/background.txt',pyvars[3])
        savetxt(od+'/numSeqs.txt',pyvars[0])
        to = datalib.trainData(ctypes.byref(ds),ctypes.c_int(params['-minMode']),alpha,pcReads,ctypes.c_int(seed),background,posreadsback,negreadsback,mWidth,ctypes.c_int(minWidth),ctypes.c_int(maxWidth),ctypes.c_int(defaultPositiveOffset),ctypes.c_int(defaultNegativeOffset),ctypes.c_int(params['-rev']),likefilename,bestmodelparams)
        pyTo = getTrainOut(to,params['-minMode'],noOfSeqs)
        print 'Now converting trainout in to python variables'
        sv.saveDetails(pyvars[0],pyvars[1],pyvars[2],od,params['-rev'],pyTo)
        x=datalib.freeTo(to)

if __name__=='__main__':
    try:
        paramVals = gd.fetchParameters(sys.argv)
        print "Loading the data ..."
        rp.getReadFiles(paramVals)
        if not(os.path.isfile(paramVals['-p'])) or not(os.path.isfile(paramVals['-n'])):
            print "One or both of the read files didn't get created"
            exit()
        print "Reads data loaded"
        gd.storeSettings(paramVals)
        pyvars = getValues(paramVals)
        pyvars.append(paramVals)
        callTraindata_parallel(pyvars)
        #callTraindata(pyvars)
    except (KeyboardInterrupt, SystemExit):
        os.system('setterm -cursor on')
        exit(2)
    exit(0)
