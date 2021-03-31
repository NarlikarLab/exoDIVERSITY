import sys,os
import numpy as np
import collections
#sys.path.append('/home/anushua/ReadsAndMotifs/scripts')
import getHTML as gh
import weblogoMod.weblogolib as wl
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from loadModel import loadmodel
import scoring as sc

def createLogo(sequences, filename):
    seqs = wl.read_seq_data(sequences)
    data = wl.LogoData.from_seqs(seqs)
    options = wl.LogoOptions()
    options.logo_title = str(len(sequences))+" sequences"
    options.size = "large"
    options.color_scheme = wl.colorscheme.monochrome
    formt = wl.LogoFormat(data, options)
    fout = open(filename, "w")
    wl.png_formatter(data, formt, fout)
    fout.close()

def getSeq(seqInNum):
    #convert the sequence in numeric form to string form
    d = {0:'A',1:'C',2:'G',3:'T'}
    seq = ''
    for s in seqInNum:
        seq+=d[s]
    return seq

def revcomp(seqs):
    d = {'A':'T','C':'G','G':'C','T':'A'}
    seqlen = len(seqs[0])
    noofseqs = len(seqs)
    rcseqs = np.zeros(noofseqs,dtype='S'+str(seqlen))
    for i in range(noofseqs):
        rcseqs[i] = ''.join(reversed(map((lambda x:d[x]),seqs[i])))
    return rcseqs

def createFiles(seqarr,posn,seqmodes,motifWidths,preadsStart,nreadsStart,preadsWidth,nreadsWidth,reldist,outdir,revFlag):
    opdir = '/'.join(outdir.split('/')[:-2])
    headerfile = opdir+'/headers.txt'
    hf = open(headerfile,'r')
    allheaders = np.zeros(len(seqarr),dtype=object)
    c = 0
    for h in hf:
        h = h.strip('\n')
        allheaders[c] = h
        c+=1
    hf.close()
    motifsbymode = dict()
    bedcoords = collections.OrderedDict()
    outmotifwr = outdir+"/info.txt"
    of = open(outmotifwr,'w')
    of.write('#sequence\tmode\tmotifStart\tposReadsStart\tnegReadsStart\tstrand\tmotif\n')
    for i in range(len(seqarr)):
        if posn[i]==-1 or seqmodes[i]<0:
            l = '\t'.join([allheaders[i],"-1","-1","-1","-1","-1","-1"])
            of.write(l+'\n')
            continue
        motif = seqarr[i][posn[i]:posn[i]+motifWidths[seqmodes[i]]]
        motif = getSeq(motif)            
        rd = reldist[seqmodes[i]]
      
        if revFlag and posn[i] >= (len(seqarr[i])/2):
            L = len(seqarr[i])
            motpos = L - posn[i]-motifWidths[seqmodes[i]]

            nrs = L - preadsStart[i]-preadsWidth[seqmodes[i]]
            prs = L - nreadsStart[i]-nreadsWidth[seqmodes[i]] 
            strand =  '-'
        else:
            motpos = posn[i]
            prs = preadsStart[i]
            nrs = nreadsStart[i]
            strand = '+'
        
        ss = allheaders[i].split(":")
        coords = map(int,ss[1].split("-"))
        ms = coords[0] + motpos
        me = ms + motifWidths[seqmodes[i]]
        
        if seqmodes[i] not in motifsbymode:
            motifsbymode[seqmodes[i]]=[motif]
            bedcoords[seqmodes[i]]=[(ss[0],str(ms),str(me),strand)]
            
        else:
            motifsbymode[seqmodes[i]]+= [motif]
            bedcoords[seqmodes[i]] += [(ss[0],str(ms),str(me),strand)]

        l = '\t'.join([allheaders[i],str(seqmodes[i]),str(motpos),str(prs),str(nrs),strand,motif])
        of.write(l+'\n')

    bf = open(outdir+"/events.bed","w")
    bf.write("track name=exoDiversity_output description=events in model with "+str(len(motifsbymode))+" modes\n")
    for m in sorted(bedcoords.keys()):
        for i in range(len(bedcoords[m])):
            bf.write("\t".join([bedcoords[m][i][0],bedcoords[m][i][1],bedcoords[m][i][2],str(m),'.',bedcoords[m][i][3]])+"\n")
    bf.close()

    for i in motifsbymode.keys():
        if i<0: continue
        createLogo(motifsbymode[i],outdir+"/logo_"+str(i)+".png")
        of = open(outdir + '/motif_'+str(i)+'.txt','w')
        of.write('\n'.join(motifsbymode[i]))
        of.close()
        createLogo(revcomp(motifsbymode[i]),outdir+"/logo_"+str(i)+"_rc.png")

def makeplots(posreads,negreads,preadsWidth,nreadsWidth,preadsStart,nreadsStart,labels,startPos,motifWidths,relDist,outdir):
    n = len(posreads)
    headerfile = '/'.join(outdir.split('/')[:-2]) + '/headers.txt'
    hf = open(headerfile,'r')
    allheaders = np.zeros(n,dtype=object)
    c = 0
    for h in hf:
        h = h.strip('\n')
        allheaders[c] = h
        c+=1
    hf.close()
    
    '''
    outfile = outdir +'/readsInfo.txt'
    of = open(outfile,'w')
    of.write("#sequences\tMode\tPosReadsStart\tNegReadStart\tPosReadsWidth\tNegReadsWidth\n")
    for i in range(n):
        line = allheaders[i]+"\t"+str(labels[i])+"\t"+str(preadsStart[i])+"\t"+str(nreadsStart[i])+"\t"+str(preadsWidth[labels[i]])+"\t"+str(nreadsWidth[labels[i]])+"\n"
        of.write(line)
    of.close()
    '''
    posreadsPerMode = dict()
    negreadsPerMode = dict()
    totlength = dict()
    plotMotStart={}

    for i in range(n):
        if labels[i]<0 or startPos[i]<0: continue
        s = 0
        if relDist[labels[i]][0]<0:
            if (nreadsStart[i] - startPos[i]) < nreadsWidth[labels[i]]:
                start = nreadsStart[i] - nreadsWidth[labels[i]]
                s = startPos[i] - (nreadsStart[i] - nreadsWidth[labels[i]])
            else:
                start = startPos[i]
                s = 0
            d = 0
        else:
            start = preadsStart[i]
            s = startPos[i]-preadsStart[i]
            d = relDist[labels[i]][0]

        if relDist[labels[i]][1] < motifWidths[labels[i]]:
            if startPos[i]+motifWidths[labels[i]] - preadsStart[i] < preadsWidth[labels[i]]:
                end = preadsStart[i] + preadsWidth[labels[i]]
            else:
                end = start + d + motifWidths[labels[i]]
        else:
            end = start + d + relDist[labels[i]][1]
            
        if labels[i] not in posreadsPerMode.keys():
            posreadsPerMode[labels[i]]=[posreads[i][start:end]]
            negreadsPerMode[labels[i]]=[negreads[i][start:end]]
            totlength[labels[i]] = end-start
            plotMotStart[labels[i]] = s
            
        else:
            posreadsPerMode[labels[i]].append(posreads[i][start:end])
            if (len(posreads[i][start:end])!=totlength[labels[i]]):
                print 'PosreadsStart: ',preadsStart[i],'startpos: ',startPos[i],'negreadsStart: ',nreadsStart[i],'reldist ',relDist[labels[i]]
                exit(0)

            negreadsPerMode[labels[i]].append(negreads[i][start:end])

    #print "Reldist ",relDist
    for m in relDist.keys():
        if m == -1: continue
        readsPlot = outdir+'/reads_'+str(m)+'.png'
        try:
            colmeansPos = np.mean(np.array(posreadsPerMode[m]),axis=0)
            colmeansNeg = np.mean(np.array(negreadsPerMode[m]),axis=0)
        except (ValueError):
            print "ERROR: Reads window size problem in mode %d**************************************"%m 
            print "Positive reads start",[(i,x) for (i,x) in enumerate(preadsStart) if labels[i]==m]
            print relDist
            exit()
        ymax = max((max(colmeansPos),max(colmeansNeg)))
        fig,ax = plt.subplots(nrows=1,ncols=1)
        g = ax.plot(range(totlength[m]),colmeansPos,color='red',linewidth=1.5,label="Forward Strand")
        g = ax.plot(range(totlength[m]),colmeansNeg,color='blue',linewidth=1.5,label="Reverse Strand")
        g = ax.set_ylabel('Mean of read counts')
        g = ax.set_ylim(0,ymax)
        g = ax.set_xlim(0,totlength[m])
        
        g = ax.axvline(x=plotMotStart[m],color='m',linestyle='--') # motstart
        g = ax.axvline(x=plotMotStart[m]+int(motifWidths[m]),color='m',linestyle='--') #motend

        '''
        if relDist[m][0]>0:
            g = ax.axvline(x=relDist[m][0],color='m',linestyle='--') # motstart
            g = ax.axvline(x=relDist[m][0]+int(motifWidths[m]),color='m',linestyle='--') #motend
            #g = ax.axvline(x = 0,color='k',linestyle='--') #posreldist first
            #g = ax.axvline(x = preadsWidth[m],color='k',linestyle='--') #posreldist last
        else:
            g = ax.axvline(x=0,color='m',linestyle='--') # motstart
            g = ax.axvline(x=int(motifWidths[m]),color='m',linestyle='--') #motend
            #g = ax.axvline(x = abs(relDist[m][0]),color='k',linestyle='--') #posreldist first
            #g = ax.axvline(x = abs(relDist[m][0])+preadsWidth[m],color='k',linestyle='--') #posreldist last
        '''
        g = ax.legend(('Forward Strand','Reverse Strand'),loc = 'best')
        fig.savefig(readsPlot)
        plt.close(fig)
        del colmeansPos, colmeansNeg

def saveDetails(sequences,posreads,negreads,outdir,twobit,revflag,trainout):
    likefile = outdir + '/likelihood.txt'
    if os.path.isfile(likefile):
        imagefile = outdir + '/likelihood.png'
        if not(os.path.isfile(imagefile)):
            if (not(os.path.isfile(likefile))):
                pass
            else:
                os.system('Rscript plotLikelihood.r '+likefile+' '+imagefile)
                os.system('rm '+ likefile)

    n = len(sequences)
    relDist = dict()
    for i in range(n):
        if trainout['labels'][i]<0 or trainout['startPos'][i]<0:
            continue
        if trainout['labels'][i] not in relDist.keys():
            posrelDist = trainout['startPos'][i]-trainout['preadsStart'][i]
            negrelDist = trainout['nreadsStart'][i]+trainout['nreadsWidth'][trainout['labels'][i]]-trainout['startPos'][i]
            relDist[trainout['labels'][i]] = [posrelDist,negrelDist]
    createFiles(sequences,trainout['startPos'],trainout['labels'],trainout['motifWidth'],trainout['preadsStart'],trainout['nreadsStart'],trainout['preadsWidth'],trainout['nreadsWidth'],relDist,outdir,revflag)

    validLabels = set(trainout['labels'])

    bestModel = outdir +"/bestModel.out"
    bm = open(bestModel,'w')
    bm.write("Posterior: "+str(trainout['likelihood'])+"\n")
    bm.write("Labels:\n"+"\t".join([str(x) for x in trainout['labels']])+"\n")
    bm.write("StartPos:\n"+"\t".join([str(x) for x in trainout['startPos']])+"\n")
    bm.write("MotifWidth:\n"+"\n".join(["Mode_"+str(i)+": "+str(w) for (i,w) in enumerate(trainout['motifWidth']) if i in validLabels])+"\n")
    bm.write("Pos Reads Width:\n"+"\n".join(["Mode_"+str(i)+": "+str(w) for (i,w) in enumerate(trainout['preadsWidth']) if i in validLabels])+"\n")
    bm.write("Neg Reads Width:\n"+"\n".join(["Mode_"+str(i)+": "+str(w) for (i,w) in enumerate(trainout['nreadsWidth']) if i in validLabels])+"\n")
    bm.write("Relative Distances: \n")
    for m in sorted(relDist.keys()):
        if m in validLabels:
            bm.write("Mode"+str(m)+": ("+",".join(map(str,relDist[m]))+")"+"\n")
    bm.close()
    makeplots(posreads,negreads,trainout['preadsWidth'],trainout['nreadsWidth'],trainout['preadsStart'],trainout['nreadsStart'],trainout['labels'],trainout['startPos'],trainout['motifWidth'],relDist,outdir)
    

def createRPlots(outdir,minMode,maxMode,posreadsdir,negreadsdir,twobit):
    for m in range(minMode,maxMode+1):
        targetdir = outdir+"/"+str(m)+"modes/"
        os.system("python allplots.py "+targetdir+" "+posreadsdir+" "+negreadsdir+" "+twobit)


def createHTML(outdir,minMode,maxMode,twobit,fasta,bestmodel):
    bm = bestmodel.split("_")[1]
    hmflag = -1
    if twobit == '':
        hmflag = 0
    else:
        hmflag = 1
    for m in range(minMode,maxMode+1):
        gh.makehtml(outdir,m,hmflag,)
    gh.makelinksHTML(outdir,minMode,maxMode)
    gh.makeBestmodelHTML(outdir,bm,bestmodel,hmflag)

def seqsInProbSpace(alldata,modelfile,scoreOutfile,revFlag,maskReps,headerfile,outdir):
    hf = open(headerfile,'r')
    allheaders = []
    c = 0
    for h in hf:
        allheaders.append(h.strip()[1:])
        c+=1
    hf.close()
    model = loadmodel(modelfile)
    scoresFromModel = sc.scoreSequences(alldata,model,revFlag)
    nof = open(scoreOutfile,'w')
    for i in range(len(allheaders)):
        sumX = sum(scoresFromModel[i])
        finalscores = [x/sumX for x in scoresFromModel[i]]
        nof.write('\t'.join([allheaders[i]]+map(str,finalscores))+'\n')
    nof.close()
    infofile = outdir+'/info.txt'
    heatmap = outdir+'/seqsInProbSpace.png'
    os.system("Rscript plotSeqsProbHeatmap.r "+infofile+" "+scoreOutfile+" "+heatmap)
