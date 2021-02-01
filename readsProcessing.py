import os,commands

def createStrandWiseReadsdir(rf,outdir,ctrl,g,fmt):
    if ctrl:
        pbg = outdir+"/reads/control_posreads.bedGraph"
        nbg = outdir+"/reads/control_negreads.bedGraph"
        pd = outdir+"/reads/control_posreads"
        nd = outdir+"/reads/control_negreads"
    else:
        pbg = outdir+"/reads/posreads.bedGraph"
        nbg = outdir+"/reads/negreads.bedGraph"
        pd = outdir+"/reads/posreads"
        nd = outdir+"/reads/negreads"
    
    if fmt == 'BAM':
        st,out = commands.getstatusoutput("genomeCoverageBed -ibam "+rf+" -bga -strand + -5 > "+pbg)
        #print "genomeCoverageBed -ibam "+rf+" -bga -strand + -5 > "+pbg
    elif fmt=='BED':
        st,out = commands.getstatusoutput("genomeCoverageBed -i "+rf+" -g "+g+" -bga -strand + -5 > "+pbg)
    if st:
        print out
        exit()

    if fmt == 'BAM':
        st,out = commands.getstatusoutput("genomeCoverageBed -ibam "+rf+" -bga -strand - -5 > "+nbg)
        #print "genomeCoverageBed -ibam "+rf+" -bga -strand - -5 > "+nbg
    elif fmt == 'BED':
        st,out =commands.getstatusoutput("genomeCoverageBed -i "+rf+" -g "+g+" -bga -strand - -5 > "+nbg)
    if st:
        print out
        exit()
    os.system("python separateChromWise.py "+pbg+" "+pd)
    os.system("python separateChromWise.py "+nbg+" "+nd)
    st,prc = commands.getstatusoutput("awk '{sum+=$4} END {print sum}' "+pbg)
    if st: print out
    st,nrc = commands.getstatusoutput("awk '{sum+=$4} END {print sum}' "+nbg)
    if st: print out
    rc = float(prc) + float(nrc)
    os.system("rm "+pbg+" "+nbg)
   
    return rc 

def subtractCtrlBinarize(peakreadsfile,controlreadsfile,outfile,scalingfactor,thresh):
    status,out = commands.getstatusoutput("python subtractControl_nonnegative.py "+peakreadsfile+" "+controlreadsfile+" "+outfile+" "+str(scalingfactor))
    if status: print out
    contents = outfile.split('.')
    binoutfile = '.'.join(contents[:-1])+'_bin_wrt'+thresh+'.'+contents[-1]
    status,out = commands.getstatusoutput("python makebinary.py "+thresh+" "+outfile+" "+binoutfile)
    if status: print out
    return binoutfile

def binarize(outfile,thresh):
    contents = outfile.split('.')
    binoutfile = '.'.join(contents[:-1])+'_bin_wrt'+thresh+'.'+contents[-1]
    status,out = commands.getstatusoutput("python makebinary.py "+thresh+" "+outfile+" "+binoutfile)
    if status: print out
    return binoutfile

def getReadFiles(paramVals):
    rf = paramVals['-r']
    ctrl = paramVals['-ctrl']
    g = paramVals['-g']
    fmt = paramVals['-format']
    outdir = paramVals['-o']
    rcexpt = 0
    rcctrl = 0
    if rf != '':
        os.makedirs(outdir+'/reads')
        if fmt == 'BAM':
            rcexpt = createStrandWiseReadsdir(rf,outdir,0,g,fmt)
            if ctrl!='':
                rcctrl = createStrandWiseReadsdir(ctrl,outdir,1,g,fmt)
        elif fmt == 'BED':
            ext = rf[:-4]
            sortedbed = rf+".sorted"+ext
            os.system("sort -k 1,1 "+rf+" > "+sortedbed)
            rcexpt = createStrandWiseReadsdir(rf,outdir,0,g,fmt)
            if ctrl!='':
                rcctrl = createStrandWiseReadsdir(ctrl,outdir,1,g,fmt)
        
        inpbed = paramVals['-f'][:-6]+".bed"
        os.system("python fastaToBed.py "+paramVals['-f']+" "+inpbed)
        os.system("python getreads.py "+inpbed+" "+outdir+"/reads/posreads "+outdir+"/reads/posreads.txt")
        os.system("python getreads.py "+inpbed+" "+outdir+"/reads/negreads "+outdir+"/reads/negreads.txt")
        posbinfile = binarize(outdir+"/reads/posreads.txt",paramVals['-bin'])
        negbinfile = binarize(outdir+"/reads/negreads.txt",paramVals['-bin'])
        if ctrl!='':
            os.system("python getreads.py "+inpbed+" "+outdir+"/reads/control_posreads "+outdir+"/reads/control_posreads.txt")
            os.system("python getreads.py "+inpbed+" "+outdir+"/reads/control_negreads "+outdir+"/reads/control_negreads.txt")
            scalingfactor = float(rcexpt)/rcctrl
            posoutfile = outdir+"/reads/posreads_controlSubtracted.txt"
            negoutfile = outdir+"/reads/negreads_controlSubtracted.txt"
            posbinfile = subtractCtrlBinarize(outdir+"/reads/posreads.txt",outdir+"/reads/control_posreads.txt",posoutfile,scalingfactor,paramVals['-bin'])
            negbinfile = subtractCtrlBinarize(outdir+"/reads/negreads.txt",outdir+"/reads/control_negreads.txt",negoutfile,scalingfactor,paramVals['-bin'])
        else:
            print 'Reads file is not in correct format'
            
    elif paramVals['-p']!='' and paramVals['-n']!='':
        if paramVals['-bin'] != 'keep':
            posbinfile = binarize(paramVals['-p'],paramVals['-bin'])
            negbinfile = binarize(paramVals['-n'],paramVals['-bin'])
        else:
            posbinfile = paramVals['-p']
            negbinfile = paramVals['-n']
    # finally keep the xxxbin file in paramVals['-xxxx']
    paramVals['-p'] = posbinfile
    paramVals['-n'] = negbinfile
