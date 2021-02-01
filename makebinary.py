import sys
import numpy as np

def wholedataMean(thresh,readsfile,binreadsfile):
    rf = open(readsfile,'r')
    allreadcounts = []
    for l in rf:
        l = l.strip().split('\t')
        allreadcounts += map(float,l[3:])
    rf.seek(0)
    mean = np.mean(allreadcounts)
    median = np.median(allreadcounts)
    q3 = np.percentile(allreadcounts,75)
    q1 = np.percentile(allreadcounts,25)
    if thresh == 'median':
        threshVal = median
    elif thresh == 'Q1':
        threshVal = q1
    elif thresh == 'Q3':
        threshVal = q3
    
    #print 'Median is :%f'%median,'Threshold is',threshVal, 'Mean is: ',mean
    brf = open(binreadsfile,'w')
    for rs in rf:
        rs = rs.strip().split('\t')
        i=3
        while i<len(rs):
            if rs[i]=='':
                i=i+1
                continue
            if (float(rs[i])>threshVal):
                rs[i]="1"
            else:
                rs[i]="0"
            i+=1
        brf.write('\t'.join(rs)+'\n')
    rf.close()
    brf.close()

if __name__ == '__main__':
    thresh = sys.argv[1]
    readsfile = sys.argv[2]
    binreadsfile = sys.argv[3]
    wholedataMean(thresh,readsfile,binreadsfile)

