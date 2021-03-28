import numpy as np
import math

modesPrior = 1.0

def scoreReads(model,data,mode,start,posreldist,negreldist,index):
    tot = model['modeSeqCounts'][mode]+sum(model['pcReads'])
    s1 = 1.0
    j = start-posreldist
    i = 0
    r1 = model['allpospwms'][mode]
    while i < model['posreadsWidths'][mode]:
        rc = int(data[1][index][j])
        s1 = s1*((r1[rc][i] + model['pcReads'][rc])/tot)/data[4][index][j]
        j = j+1
        i = i+1
    s2 = 1.0
    r2 = model['allnegpwms'][mode]
    j = start + negreldist - model['negreadsWidths'][mode]
    i = 0
    while i < model['negreadsWidths'][mode]:
        rc = int(data[2][index][j])
        s2 = s2*((r2[rc][i] + model['pcReads'][rc])/tot)/data[5][index][j]
        j = j+1
        i = i+1
    s = s1*s2 
    return s

def scoreMotif(model,data,mode,start,index):
    s = 1.0
    m1 = model['allpwms'][mode]
    for i in range(model['motifWidths'][mode]):
        s = s*((m1[ data[0][index][start+i] ][i]+model['alpha'])/(model['modeSeqCounts'][mode] + model['features']*model['alpha']))
        if (data[3][index][start+i] < 0.00001): continue
        s = s/data[3][index][start+i]
    return s

def scoreSequences(data,model,revFlag):
    ####### data = [seqs,readfreq_pos,readfreq_neg,background,posreadsback,negreadsback,lookahead]
    ####### model = [noOfmodes,features,totseqs,allpwms,allpospwms,allnegpwms,modeSeqCounts,motifWidths,posreadsWidths,negreadsWidths,reldist,alpha,pcReads]
    seqscores = np.zeros(len(data[0]),dtype=object)
    for index in range(len(data[0])): # over all seqs in data
        modeWiseScores = np.zeros(model['noOfmodes'],dtype=float)
        for i in range(model['noOfmodes']): # over all modes
            if model['modeSeqCounts'][i] == 0:
                continue
            sumscore = 0.0
            L = len(data[0][index])

            preldist = model['reldist'][i][0]
            nreldist = model['reldist'][i][1]

            if (nreldist < model['motifWidths'][i]):
                rhs = model['motifWidths'][i]
            else:
                rhs = nreldist
            if (preldist<=0):
                j = 0
            else:
                j = preldist
            while j < (L-rhs+1):
                if (revFlag):
                    if ((preldist > 0) and (j > L/2 -rhs) and (j < L/2 + preldist)):
                        j=j+1
                        continue
                    if ((preldist<=0) and (j > L/2-rhs) and (j < L/2)):
                        j = j+1
                        continue
                if data[6][index][j] < model['motifWidths'][i]:
                    v = data[6][index][j]
                    j = j+v+1
                    continue
                rval = scoreReads(model,data,i,j,preldist,nreldist,index)
                motval = scoreMotif(model,data,i,j,index)
                currscore = motval*rval
                sumscore += currscore
                # if currscore > maxscore:
                #     maxmotval = motval
                #     maxrval = rval
                #     maxscore = currscore
                #     maxindex = j
                j = j+1
            modeWiseScores[i] = ((model['modeSeqCounts'][i]+modesPrior)/(model['totseqs']+model['noOfmodes']*modesPrior))*(1.0/L)*sumscore
        seqscores[index] = modeWiseScores
    return seqscores

def scoreSeqBack(dnaseq,posreadsseq,negreadsseq):
    logsum = 0.0
    for i in range(len(dnaseq)):
        if dnaseq[i]<=0:
            continue
        dna = math.log(dnaseq[i])
        pr = math.log(posreadsseq[i])
        nr = math.log(negreadsseq[i])
        logsum += (dna+pr+nr)
    return logsum

def scoreWithBackground(data,revFlag):
    ####### data = [seqs,readfreq_pos,readfreq_neg,background,posreadsback,negreadsback,lookahead]
    seqscores = np.zeros(len(data[0]),dtype=float)
    for index in range(len(data[0])): # over all seqs in data
        totback = scoreSeqBack(data[3][index],data[4][index],data[5][index])
        seqscores[index] = totback
    return seqscores

