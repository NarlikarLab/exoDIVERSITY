import sys,os
import numpy as np
def getReadCounts_chrom(chromfile):
    cf = open(chromfile,'r')
    # Make two arrays one storing the start positions
    # and the other storing the readcounts
    count=0
    for l in cf:
        count+=1
    cf.seek(0)
    startPositions = np.zeros(count,dtype=int)
    readCounts = np.zeros(count,dtype=int)
    count = 0
    for l in cf:
        l1 = l.strip('\n').split('\t')
        sp = int(l1[1])
        rc = int(l1[3])
        startPositions[count] = sp
        readCounts[count] = rc
        count +=1
    cf.close()
    return startPositions, readCounts

def interval_binsearch(low,high,num,arr):
    if num >= arr[-1]:
        return len(arr)-1
    while (low<=high):
        mid = (low+high)/2
        if (num >= arr[mid] and num < arr[mid+1]):
            return mid
        elif (num<arr[mid]):
            high = mid-1
        elif (num > arr[mid]):
            low=mid+1
    print "Couldn't get an interval"
    return -1

def fillArray(posWiseReads,startPositions,readCounts,seqstart,end):
    L=len(startPositions)
    # Assuming the startPositions are sorted
    i = interval_binsearch(0,(L-1),seqstart,startPositions)
    if i==-1:
        print seqstart,end
        exit()
    count = 0
    start = seqstart
    while i+1 < L:
        if  startPositions[i+1] <= end:
            # Fill those many positions from count to count +startPositions[i+1]-start
            for c in range(count,(count+startPositions[i+1]-start)):
                posWiseReads[c] = readCounts[i]
            count = count+(startPositions[i+1]-start)
            start = startPositions[i+1]
            i=i+1
        else:
            break
    for x in range(count,len(posWiseReads)):
        posWiseReads[x] = readCounts[i]
    #return posWiseReads

def main():
    bedfile = sys.argv[1]
    chromdir = sys.argv[2]
    outfile = sys.argv[3]
    
    bf = open(bedfile,'r')
    bl = bf.readline()
    chromWiseSeqs = dict()
    counter = 0
    while bl:
        bl = bl.strip('\n').split('\t')
        chrom = bl[0]
        start = int(bl[1])
        end = int(bl[2])
        if chrom not in chromWiseSeqs.keys():
            chromWiseSeqs[chrom]=dict()
        length = abs(end-start)
        readcounts = np.zeros(length,dtype=int)
        chromWiseSeqs[chrom][(start,end,counter)]=readcounts
        counter += 1
        bl = bf.readline()
    bf.close()
    vals=[]
    for chrom in chromWiseSeqs.keys():
        seqs = chromWiseSeqs[chrom]
        startPositions,readCounts = getReadCounts_chrom(chromdir+'/'+chrom+'.bedGraph')
        for seqdetails in seqs.keys():
            start = seqdetails[0]
            end = seqdetails[1]
            fillArray(seqs[seqdetails],startPositions,readCounts,start,end)
            vals.append((chrom,start,end,seqdetails[2]))
    
    # Write to file
    of = open(outfile,'w')
    for metadata in sorted(vals,key=lambda x: x[3]):
        chrom = metadata[0]
        start = metadata[1]
        end = metadata[2]
        counter = metadata[3]
        arr = chromWiseSeqs[chrom][(start,end,counter)]
        string = chrom+'\t'+str(start)+'\t'+str(end)+'\t'+('\t'.join(map(str,arr)))
        of.write(string+'\n')
    of.close()

main()
