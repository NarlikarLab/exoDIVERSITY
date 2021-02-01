import sys
import numpy as np
import re
from itertools import izip
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pylab import savefig
import seaborn as sns; sns.set()

def loadInfo(infofile,dflag):
    inf = open(infofile,'r')
    startpos = []
    labels = []
    widths = dict()
    reldist = dict()
    noOfSeqs = 0
    il = inf.readline()
    for il in inf:
        noOfSeqs += 1
        il = il.strip().split('\t')
        se = map(int,re.split(r'[:|-]',il[0])[1:])
        seqlength = int(se[1]-se[0]) * 2
        label = il[1]
        labels.append(label)
        motStart = int(il[2])
        if dflag:
            strand = il[3]
            motwidth = int(il[4])
            if strand == '-':
                motStart = seqlength - motStart - motwidth
            posreadsStart = motstart - 20
            negreadsStart = motstart + motwidth
            prd = motStart - 20
            nrd = negreadsStart+20 - motStart
        else:
            strand = il[7]
            # if strand == '-':
            #     motStart = seqlength - motStart - motwidth
            #     posreadsStart = seqlength - int(il[3]) - int(il[5])
            #     negreadsStart = seqlength - int(il[4]) - int(il[6])
            # else:
            if strand == '+':
                motwidth = len(il[-1])
                posreadsStart = int(il[3])
                negreadsStart = int(il[4])
                prd = motStart - posreadsStart
                nrd = negreadsStart + int(il[6]) - motStart
        if label not in widths.keys() and strand == '+':
            if dflag: widths[label] = (motwidth,20,20)
            else: widths[label] = (motwidth,int(il[5]),int(il[6])) # motwidth, posreadswidth, negreadwidth
            reldist[label] = (prd,nrd)
        startpos.append(motStart)

    inf.close()
    return (noOfSeqs,reldist,widths,startpos,labels)

def makeplots(info,posreadsfile,negreadsfile,outdir):
    noOfSeqs = info[0]
    reldist = info[1]
    widths = info[2]
    startpos = info[3]
    labels = info[4]
    posreads = np.zeros(noOfSeqs,dtype=object)
    negreads = np.zeros(noOfSeqs,dtype=object)
    i = 0
    with open(posreadsfile) as prf, open(negreadsfile) as nrf:
        for pl,nl in izip(prf,nrf):
            pl = [int(x) for x in pl.strip().split('\t')[3:]]
            nl = [int(x) for x in nl.strip().split('\t')[3:]]
            posreads[i] = pl + [y for y in reversed(nl)]
            negreads[i] = nl + [y for y in reversed(pl)]
            i+=1
    posreadsPerMode = dict()
    negreadsPerMode = dict()
    for i in range(noOfSeqs):
        if labels[i] < 0: continue
        (prd,nrd) = reldist[labels[i]]
        if labels[i] not in posreadsPerMode.keys():
            posreadsPerMode[labels[i]]=[posreads[i][(startpos[i]-prd):(startpos[i]+nrd)]]
            negreadsPerMode[labels[i]]=[negreads[i][(startpos[i]-prd):(startpos[i]+nrd)]]
        else:
            posreadsPerMode[labels[i]].append(posreads[i][(startpos[i]-prd):(startpos[i]+nrd)])
            negreadsPerMode[labels[i]].append(negreads[i][(startpos[i]-prd):(startpos[i]+nrd)])

    for m in reldist.keys():
        readsplot = outdir + '/OrigReads_'+str(m)+'.png'
        try:
            colmeansPos = np.mean(np.array(posreadsPerMode[m]),axis=0)
            colmeansNeg = np.mean(np.array(negreadsPerMode[m]),axis=0)
            hm = outdir + '/heatmapSox2posreads_'+str(m)+'.png'
            ax1 = plt.figure()
            ax1 = sns.heatmap(np.array(posreadsPerMode[m]),cmap="YlGnBu")
            ax1.get_figure().savefig(hm,dpi=400)
            del ax1
            hm = outdir + '/heatmapSox2negreads_'+str(m)+'.png'
            ax1 = plt.figure()
            ax1 = sns.heatmap(np.array(negreadsPerMode[m]),cmap="bwr")
            ax1.get_figure().savefig(hm,dpi=400)
            del ax1
            # print 'For mode: ',m
            # print 'Posreads max: ',np.amax(posreadsPerMode[m])
            # print 'Negreads max: ',np.amax(negreadsPerMode[m])
        except(ValueError):
            print "ERROR: Reads window size problem in mode %d ****************************"%m
            exit()
        # ymax = max(max(colmeansPos),max(colmeansNeg))
        # fig,ax = plt.subplots(nrows=1,ncols=1)
        # g = ax.plot(range(sum(reldist[m])),colmeansPos,color='red',linewidth=1.5,label="Forward strand")
        # g = ax.plot(range(sum(reldist[m])),colmeansNeg,color='blue',linewidth=1.5,label="Reverse strand")
        # a = ax.set_ylabel('Mean of read counts')
        # g = ax.set_ylim(0,ymax)
        # g = ax.set_xlim(0,sum(reldist[m]))

        # g = ax.axvline(x=reldist[m][0],color='k',linestyle='--')
        # g = ax.axvline(x=reldist[m][0]+widths[m][0],color='k',linestyle='--')
        # g = ax.axvline(x=5,color='k',linestyle='--')
        # g = ax.axvline(x=sum(reldist[m])-5,color='k',linestyle='--')

        # g = ax.legend(('Forward strand','Reverse strand'),loc='best')
        # fig.savefig(readsplot)
        # plt.close(fig)
        # del colmeansPos,colmeansNeg
        
       
if __name__=='__main__':
    infofile = sys.argv[1]
    posreadsfile = sys.argv[2]
    negreadsfile = sys.argv[3]
    outdir = sys.argv[4]
    divflg = 0 # the posreads and negreads file should have a 20 window extra at each end
    info = loadInfo(infofile,divflg)
    makeplots(info,posreadsfile,negreadsfile,outdir)
