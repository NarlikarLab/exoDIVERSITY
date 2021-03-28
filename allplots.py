import sys
import getAlignmentRegions as ar
import plotreadsHeatmap as rh
import os

def motifsAscendingOrder(infofile):
    info = open(infofile,'r')
    l = info.readline()
    modeseqs = {}
    for l in info:
        l=l.strip().split('\t')
        mode = int(l[1])
        if mode not in modeseqs:
            modeseqs[mode]=1
        else:
            modeseqs[mode]+=1
    info.close()
    ascendingOrder = sorted(modeseqs.items(),key=(lambda x:x[1]))
    return ascendingOrder


if __name__=='__main__':
    outdir = sys.argv[1]
    posreadsDir = sys.argv[2]
    negreadsDir = sys.argv[3]
    twobitfile = sys.argv[4]
    
    infofile = outdir + '/info.txt'
    aligndetails = outdir + '/alignmentDetails.txt'
    side = 25
    pwmOrder = motifsAscendingOrder(infofile)
    pwmDescOrder = [x[0] for x in reversed(pwmOrder)]
    pwmDescOrder = sorted(pwmDescOrder,reverse=True)

    refpwmmode = pwmDescOrder[0]
    st = '+'
    # shift[mode] = ('L/R',shiftBy,'+/-')
    shift = {}
    for m in pwmDescOrder:
        shift[m]=['N',0,'+']

    ad = open(aligndetails,'w')
    for m in shift.keys():
        ad.write('\t'.join([str(m),shift[m][2],str(shift[m][1]), shift[m][0]])+'\n')
    ad.close()

    ar.getRegions(infofile,shift,refpwmmode,side,twobitfile,outdir)
    ar.createPlots(outdir,pwmDescOrder)
    #print "Alignment is done"

    if not(os.path.isdir(posreadsDir)) or not(os.path.isdir(negreadsDir)):
        pass
    else:
        rh.getReadsHeatmap(outdir,posreadsDir,negreadsDir,aligndetails,side,pwmDescOrder)
    #print "Reads plotted"
