#!/usr/bin/python
import sys,re
def main():
    fastafile = sys.argv[1]
    outbedfile = sys.argv[2]
    ff = open(fastafile,'r')
    of = open(outbedfile,'w')
    for l in ff:
        if l=='':
            continue
        if l[0]=='>':
            l1 = l.strip('\n')
            l1 = re.split('>|:|-',l1)
            bedline = '\t'.join(l1[1:])+'\n'
            of.write(bedline)
    ff.close()
    of.close()
main()
