### To copy the files from the run with highest likelihood in the parent directory

### To find the best model and create a link

import sys
import os
import re 
import numpy as np
import collections
import math
from loadModel import *

def chooseBestTrial(outdir,minMode,maxMode,noOfmodels):
    maxpostDetails = dict()
    for i in range(minMode,maxMode+1):
        if not(os.path.isdir(outdir+'/'+str(i)+'modes')):
            continue
        maxpost = -9999999999
        bestrun = 0
        lenOflabels = -1
        for j in range(1,noOfmodels+1):
            bestmodelfile = outdir+'/'+str(i)+'modes/run'+str(j)+'/bestModel.out' 
            if (not(os.path.isfile(bestmodelfile))):
                continue
            best = open(bestmodelfile)
            l = best.readline()
            post = float (l.strip().split(' ')[1])
            labels = best.readline()
            labels = set(best.readline().strip().split('\t'))
            if '-1' in labels:
                lenOflabels = len(labels)-1
            else:
                lenOflabels = len(labels)
            if lenOflabels != i:
                continue
            if post > maxpost:
                maxpost = post
                bestrun = j
            best.close()
            likefile = outdir + '/' + str(i)+'modes/run'+str(j)+'/likelihood.txt'
            imagefile = outdir + '/' + str(i)+'modes/run'+str(j)+'/likelihood.png'
            if not(os.path.isfile(imagefile)):
               if (not(os.path.isfile(likefile))):
                   pass
               else:
                   os.system('Rscript plotLikelihood.r '+likefile+' '+imagefile)
                   os.system('rm '+ likefile)
        if bestrun == 0:
            print "Check the bestModel files for mode: "+str(i)
            os.system('rm '+outdir+'/'+str(i)+'modes/*.png '+outdir+'/'+str(i)+'modes/*.txt')
            continue
        else:
            os.system("cp "+outdir+'/'+str(i)+'modes/run'+str(bestrun)+'/* ' + outdir+'/'+str(i)+'modes')
            rundirs = [outdir+"/"+str(i)+"modes/run"+str(x) for x in range(1,noOfmodels+1)]
            os.system("rm -r "+" ".join(rundirs))
        maxpostDetails[i] = (bestrun,maxpost)
    return maxpostDetails

def getBestModelBIC(outdir,minMode,maxMode,noOfmodels,maxpostdetails):
    BICscores = dict()
    minscore = 999999999
    bestmodel = 0
    for m in sorted(maxpostdetails.keys()):
        bt = outdir + '/' +str(m) +'modes/bestModel.out'
        if not(os.path.isfile(bt)):
            continue
        f = open(bt,'r')
        l = f.readline()
        l1 = l.strip().split(': ')
        loglike = float(l1[1])
        labels = f.readline()
        labels = map(int,f.readline().strip().split('\t'))
        noOfSeqs = len(labels)
        modes = collections.Counter(labels)
        noOfmodes = len(modes.keys())
        modesPrior = 1
        modespostterm = sum([modes[i]*(math.log(modes[i]+modesPrior)-math.log(noOfSeqs + modesPrior*noOfmodes)) for i in modes.keys()])
        bestmodelparams = outdir+"/"+str(m)+"modes/bestModelParams.txt"
        model = loadmodel(bestmodelparams)
        
        for _ in range(3):
            l = f.readline()
        motiflengths = []
        l = f.readline()
        while l:
            l1 = l.strip().split(': ')
            if len(l1)<2:
                break
            motiflengths.append(int(l1[1]))
            l = f.readline()
        readlengths = []
        l = f.readline()
        while l:
            if l[0] == 'N':
                break
            l1 = l.strip().split(': ')
            readlengths.append(int(l1[1]))
            l = f.readline()
        f.close()
        noOfMotifs = len(motiflengths)
        ## BIC: k*ln(n) - 2*ln(L)
        # k: the number of parameters estimated by the model
        # n: the sample size
        noOfMotifParams = 3
        noOfReadParams = 1
        k = 0
        for i in range(noOfMotifs):
            if i in modes:
                k += (motiflengths[i]*noOfMotifParams + 2*readlengths[i] + 2)
        
        score = k * math.log(noOfSeqs) - 2*loglike
        print 'Mode: ',m,'\tlikelihood: ',loglike, '\tscore: ',score
        
        BICscores[m] = score
        if score < minscore:
            minscore = score
            bestmodel = str(m)+'modes'

    print 'Min BIC score: ',minscore,'Best Model: ',bestmodel
    if bestmodel == 0:
        return 0
    else:
        linkfile = 'bestModel_'+str(bestmodel)
        cwd = os.getcwd()
        os.chdir(outdir)
        os.system('ln -s '+bestmodel+' '+linkfile)
        os.chdir(cwd)
        return linkfile

def fetchBestModel(outdir,minMode,maxMode,noOfmodels):
    maxpostdetails = chooseBestTrial(outdir,minMode,maxMode,noOfmodels)
    bestmodel = getBestModelBIC(outdir,minMode,maxMode,noOfmodels,maxpostdetails)
    return bestmodel
