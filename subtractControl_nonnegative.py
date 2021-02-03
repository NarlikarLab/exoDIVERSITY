import sys
import numpy as np
def main():
    #for the same sequence PeakReads - scalingfactor*ControlReads
    peakreadsfile = sys.argv[1] 
    controlreadsfile = sys.argv[2]
    outfile = sys.argv[3]
    if len(sys.argv)<5:
        scalingfactor = 1
    else:
        scalingfactor = float(sys.argv[4])
    
    pf = open(peakreadsfile,'r')
    noOfSeqs = 0
    for l in pf:
        noOfSeqs+=1
    pf.close()
    allreads = np.zeros(noOfSeqs,dtype=object)
    j=0
    of = open(outfile,'w')
    pf = open(peakreadsfile,'r')
    cf = open(controlreadsfile,'r')
    for l1,l2 in zip(pf,cf):
        l1 = l1.strip('\n').split('\t')
        l2 = l2.strip('\n').split('\t')
        allreads[j]=np.zeros(len(l1),dtype='|S20')
        allreads[j][:3] = l1[:3]
        for i in range(3,len(l1)):
            val = float(l1[i])-float(l2[i])*scalingfactor
            if val < 0:
                allreads[j][i] = '0'
            else:
                allreads[j][i] = str(val)
        of.write('\t'.join(allreads[j])+'\n')
        j+=1
    of.close()
    pf.close()
    cf.close()
main()
