import numpy as np
import re
import commands
import os

def revcomp(seq):
    comp = {'A':'T','a':'t','C':'G','c':'g','G':'C','g':'c','T':'A','t':'a','N':'N'}
    rcseq = [comp[x] for x in reversed(seq)]
    rcseq =''.join(rcseq)
    return rcseq

def getSeq(chrom,regionStart,regionEnd,twobitfile,stpref,reportedSt,seqlen):
    (status,op) = commands.getstatusoutput('twoBitToFa -seq='+chrom+' -start='+str(regionStart)+' -end='+str(regionEnd)+' '+twobitfile+' stdout')
    seq = op.split('\n')
    header = seq[0][1:]
    seq = ''.join(seq[1:]) # made single line
    if len(seq)!=seqlen:
        print 'Sequence length is not equal to ',seqlen
        print header,seq
        exit()
    if stpref == '+':
        if reportedSt == '+':
            pass
        else:
            seq = ''.join(revcomp(seq))
    elif stpref =='-':
        if reportedSt == '+':
            seq = ''.join(revcomp(seq))
        else:
            pass
    return seq

def getInfo(infofile):
    info = open(infofile,'r')
    motifmidCoords = {}
    motifWidths = {}
    l = info.readline()
    counter = 0
    for l in info:
        l = l.strip().split('\t')
        header = re.split(r'[:|-]',l[0])
        header[1] = int(header[1])
        header[2] = int(header[2])
        mode = int(l[1])
        motstart = header[1] + int(l[2])
        if len(l[-1]) < 6:
            w =  int(l[-1])
        else:
            w = len(l[-1])
        motifmid = motstart + w/2
        strand = l[5]
        if mode not in motifWidths:
            motifWidths[mode] = w
            motifmidCoords[mode] = [[header[0],motifmid,strand,counter]]
        else:
            motifmidCoords[mode].append([header[0],motifmid,strand,counter])
        counter += 1
    info.close()
    return (motifmidCoords,motifWidths,counter)

def createPlots(outdir,modeOrder):
    if modeOrder == []:
        os.system('Rscript draw1.r '+outdir + '/fastaMatrix.txt '+outdir+'/alignedMotifs.png')
    else:
        #print 'Rscript draw1.r '+outdir+'/fastaMatrix.txt '+outdir+'/alignedMotifs.png '+' '.join(map(str,modeOrder))
        os.system('Rscript draw1.r '+outdir+'/fastaMatrix.txt '+outdir+'/alignedMotifs.png '+' '.join(map(str,modeOrder)))

def getRegions(infofile,shifts,refMode,side,twobitfile,outdir):
    modeWiseMotifmid, motifWidths,noOfSeqs = getInfo(infofile)
    ## Shift the motif mids accordingly
    refpwmWidth = motifWidths[refMode]
    for m in modeWiseMotifmid.keys():
        if m == refMode:
            continue
        if m == -1: continue
        if motifWidths[m] < refpwmWidth:
            if shifts[m][0] == 'L': # For left shift
                dist = abs(abs(shifts[m][1]) - motifWidths[m]/2)
                moveBy = (refpwmWidth/2) - dist
                if shifts[m][2] == '-':
                    moveBy = -1*moveBy
                for i in range(len(modeWiseMotifmid[m])):
                    if modeWiseMotifmid[m][i][2] == '+':
                        modeWiseMotifmid[m][i][1] += moveBy
                    elif modeWiseMotifmid[m][i][2] == '-':
                        modeWiseMotifmid[m][i][1] -= moveBy
                        
            elif shifts[m][0] == 'R': # For right shift
                dist = abs(shifts[m][1]) + motifWidths[m]/2
                moveBy = (refpwmWidth/2) - dist
                if shifts[m][2] == '-':
                    moveBy = -1*moveBy
                for i in range(len(modeWiseMotifmid[m])):
                    if modeWiseMotifmid[m][i][2] == '+':
                        modeWiseMotifmid[m][i][1] += moveBy
                    elif modeWiseMotifmid[m][i][2] == '-':
                        modeWiseMotifmid[m][i][1] -= moveBy
                        
        elif motifWidths[m] > refpwmWidth:
            if shifts[m][0] == 'L': # For left shift
                dist = abs(motifWidths[m]/2 - refpwmWidth/2)
                moveBy = abs(shifts[m][1]) - dist
                if shifts[m][2] == '-':
                    moveBy = -1*moveBy
                for i in range(len(modeWiseMotifmid[m])):
                    if modeWiseMotifmid[m][i][2] == '+':
                        modeWiseMotifmid[m][i][1] += moveBy
                    elif modeWiseMotifmid[m][i][2] == '-':
                        modeWiseMotifmid[m][i][1] -= moveBy
                        
            elif shifts[m][0] == 'R': # For right shift
                dist = abs(motifWidths[m]/2 - refpwmWidth/2)
                moveBy =  abs(shifts[m][1]) + dist
                if shifts[m][2] =='-':
                    moveBy = -1*moveBy
                for i in range(len(modeWiseMotifmid[m])):
                    if modeWiseMotifmid[m][i][2] == '+':
                        modeWiseMotifmid[m][i][1] -= moveBy
                    elif modeWiseMotifmid[m][i][2] == '-':
                        modeWiseMotifmid[m][i][1] += moveBy
        else: # motifWidths[m] == refpwmwidth
            if shifts[m][2] == '-':
                shifts[m][1] *= -1
            if shifts[m][0] == 'L':
                for i in range(len(modeWiseMotifmid[m])):
                    if modeWiseMotifmid[m][i][2] == '+':
                        modeWiseMotifmid[m][i][1] += shifts[m][1]
                    elif modeWiseMotifmid[m][i][2] == '-':
                        modeWiseMotifmid[m][i][1] -= shifts[m][1] 
            elif shifts[m][0] == 'R':
                for i in range(len(modeWiseMotifmid[m])):
                    if modeWiseMotifmid[m][i][2] == '+':
                        modeWiseMotifmid[m][i][1] -= shifts[m][1]
                    elif modeWiseMotifmid[m][i][2] == '-':
                        modeWiseMotifmid[m][i][1] += shifts[m][1]
    '''
    If strand preference is +, keep + strand regions as it is and flip the - strand regions
    If strand preference is -, flip the + strand regions and keep the - strand regions as it is
    '''
    alignedSeqs = {}
    fmVec = np.zeros(noOfSeqs,dtype=object)
    d = {'A':'0','a':'0','C':'1','c':'1','G':'2','g':'2','T':'3','t':'3','N':'4'}
    for m in modeWiseMotifmid.keys():
        if m==-1: continue
        motifwidth = motifWidths[m]
        st = shifts[m][2]
        for i in range(len(modeWiseMotifmid[m])):
            chrom = modeWiseMotifmid[m][i][0]
            motifmid = modeWiseMotifmid[m][i][1]
            reportedSt = modeWiseMotifmid[m][i][2]
            counter = modeWiseMotifmid[m][i][3]
            if motifwidth%2==0:
                regionStart = motifmid-side
                regionEnd = motifmid+side
                seqlen = 2*side
                seq = getSeq(chrom,regionStart,regionEnd,twobitfile,st,reportedSt,seqlen)
            else:
                regionStart = motifmid-side
                regionEnd = motifmid+side+1
                seqlen = side*2 + 1
                seq = getSeq(chrom,regionStart,regionEnd,twobitfile,st,reportedSt,seqlen)
                seq = seq[:-1]
            header = chrom + '\t' +str(regionStart)+'\t'+str(regionEnd)
            alignedSeqs[counter] = (header,str(m),reportedSt,str(motifwidth),seq)
            fmVec[counter] = '\t'.join([str(m)]+[d[x] for x in seq])+'\n'
 
    out = open(outdir+'/alignedMotifs.txt','w')
    fm = open(outdir+'/fastaMatrix.txt','w')
    for i in range(noOfSeqs): 
        out.write('\t'.join(alignedSeqs[i])+'\n')
        fm.write(fmVec[i])
    out.close()
    fm.close()

