import sys,os
import commands
def main():
    readsfile = sys.argv[1]
    outdir = sys.argv[2]
    if os.path.isdir(outdir):
        if(len(os.listdir(outdir))==0):
            pass
        else:
            os.system('rm -r '+outdir+'/*')
    else:
        os.makedirs(outdir)
    allchrs = dict()
    
    chromcmd = "cut -f1 "+readsfile+" | sort | uniq -c"
    st,out = commands.getstatusoutput(chromcmd)
    if st:
        print out
    else:
        out = out.split('\n')
        chroms = []
        for i in range(len(out)):
            l=filter((lambda x:x!=''),out[i].split(' '))
            chroms.append(l[-1])
    for c in chroms:
        allchrs[c] = "awk 'function abs(x){ return ((x<0.0)?-x:x) } {if($1==\""+c+"\") print $1\"\\t\"$2\"\\t\"$3\"\\t\"abs($4)}' " + readsfile + " > "+outdir + "/" + c + ".bedGraph"
        st,out = commands.getstatusoutput(allchrs[c])
        if st:
            print out
main()
