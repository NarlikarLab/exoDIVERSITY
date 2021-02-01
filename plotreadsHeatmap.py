import sys
import os
import re
import numpy as np
from itertools import izip

def getReadsHeatmap(outdir,posreadsDir,negreadsDir,aligndetails,side,modeOrder=[]):
    ad = open(aligndetails,'r')
    modeWiseStrands = {}
    for l in ad:
        if l == '\n':
            continue
        l = l.strip().split('\t')
        modeWiseStrands[int(l[0])] = l[1]
    ad.close()
    ## get and strore the chrom, regstart, regend, strand, motwidth
    alignedMotifs = 'alignedMotifs.txt'
    am = open(outdir + '/' + alignedMotifs,'r')
    count = 0
    for l in am:
        count+=1
    am.seek(0)
    seqregions = np.zeros(count,dtype=object)
    motifWidths = {}
    i = 0
    for l in am:
        l = l.strip().split('\t')
        chrom = l[0]
        regstart = int(l[1])
        regend = int(l[2])
        mode = int(l[3])
        st = l[4]
        if modeWiseStrands[mode] == '-' :
            if st == '+':
                st = '-'
            else:
                st = '+'
        motwidth = int(l[5])
        motifWidths[mode] = motwidth
        seqregions[i] = (chrom,regstart,regend,mode,st,motwidth)
        i+=1
    am.close()
    
    ## Writing the seqregions to a bed file
    regionfile = outdir + '/regionsForReadsplot.txt'
    rf = open(regionfile,'w')
    for i in range(count):
        line = seqregions[i][0]+'\t'+str(seqregions[i][1])+'\t'+str(seqregions[i][2])+'\t'+str(seqregions[i][3])+'\n'
        rf.write(line)
    rf.close()

    posreadsfile = outdir + '/posreadsForPlot.txt'
    negreadsfile = outdir + '/negreadsForPlot.txt'
    os.system('python getreads.py '+regionfile+' '+posreadsDir+' '+posreadsfile)
    os.system('python getreads_neg.py '+regionfile+' '+negreadsDir+' '+negreadsfile)
    '''
    ######### Plotting binary reads
    posreadsbin = outdir + '/posreadsForPlot_bin.txt'
    negreadsbin = outdir + '/negreadsForPlot_bin.txt'
    os.system('python makebinary_wrtMean.py '+posreadsfile+' '+posreadsbin)
    os.system('python makebinary_wrtMean.py '+negreadsfile+' '+negreadsbin)
    posreadsfile = posreadsbin
    negreadsfile = negreadsbin
    ########## plotting binary reads
    '''

    pr = open(posreadsfile,'r')
    nr = open(negreadsfile,'r')
    allposreads = np.zeros(count,dtype=object)
    allnegreads = np.zeros(count,dtype=object)
    i = 0
    for (prl,nrl) in izip(pr,nr):
        prl = prl.strip().split('\t')[3:]
        nrl = nrl.strip().split('\t')[3:]
        allposreads[i] = prl
        allnegreads[i] = nrl
        i+=1
    pr.close()
    nr.close()
    
    # According to mode and strand and motwidth accumulate all posreads and negreads
    modeWisePosreads = {}
    modeWiseNegreads = {}
    for m in motifWidths.keys():
        modeWisePosreads[m] = []
        modeWiseNegreads[m] = []
    for i in range(count):
        mode = seqregions[i][3]
        st = seqregions[i][4]
        motwidth = seqregions[i][5]
        if st == '+':
            posr = allposreads[i]
            negr = allnegreads[i]
        else:
            posr = [x for x in reversed(allnegreads[i])]
            negr = [x for x in reversed(allposreads[i])]

        if motwidth%2 != 0 :
            posr = posr[:-1]
            negr = negr[:-1]
        modeWisePosreads[mode].append((i,posr))
        modeWiseNegreads[mode].append((i,negr))
    posreadsModewise = outdir + '/posreadsModeWiseForPlot.txt'
    negreadsModewise = outdir + '/negreadsModeWiseForPlot.txt'

    prmf = open(posreadsModewise,'w')
    nrmf = open(negreadsModewise,'w')
    # for m in sorted(modeWisePosreads.keys()):
    #     for (ps,ns) in zip(modeWisePosreads[m],modeWiseNegreads[m]):
    #         prmf.write(str(m)+'\t'+'\t'.join(ps)+'\n')
    #         nrmf.write(str(m)+'\t'+'\t'.join(ns)+'\n')
    # prmf.close()
    # nrmf.close()
    preadsSorted = []
    nreadsSorted = []
    for m in modeWisePosreads.keys():
        preadsSorted += [(c,m,l) for (c,l) in modeWisePosreads[m]]
        nreadsSorted += [(c,m,l) for (c,l) in modeWiseNegreads[m]]

    preadsSorted = sorted(preadsSorted,key=lambda x:x[0])
    nreadsSorted = sorted(nreadsSorted,key=lambda x:x[0])

    for i in range(count):
        mode = preadsSorted[i][1]
        ps = preadsSorted[i][2]
        ns = nreadsSorted[i][2]
        prmf.write(str(mode)+'\t'+'\t'.join(ps)+'\n')
        nrmf.write(str(mode)+'\t'+'\t'.join(ns)+'\n')
    prmf.close()
    nrmf.close()

    if modeOrder == []:
        os.system('Rscript plotreadsHeatmap.r red '+posreadsModewise+' '+outdir+'/posreadsHeatmap.png')
        os.system('Rscript plotreadsHeatmap.r blue '+negreadsModewise+' '+outdir+'/negreadsHeatmap.png')
    else:
        #print 'COMMAND:\n Rscript plotreadsHeatmap.r red '+posreadsModewise+' '+outdir+'/posreadsHeatmap.png '+' '.join(map(str,modeOrder))
        os.system('Rscript plotreadsHeatmap.r red '+posreadsModewise+' '+outdir+'/posreadsHeatmap.png '+' '.join(map(str,modeOrder)))
        os.system('Rscript plotreadsHeatmap.r blue '+negreadsModewise+' '+outdir+'/negreadsHeatmap.png '+' '.join(map(str,modeOrder)))

def main():
    outdir = sys.argv[1]
    posreadsDir = sys.argv[2]
    negreadsDir = sys.argv[3]
    side = 50
    aligndetails = 'alignmentDetails.txt'
    getReadsHeatmap(outdir,posreadsDir,negreadsDir,aligndetails,side)
#main()

            
        
