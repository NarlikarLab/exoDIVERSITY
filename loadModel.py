import numpy as np
def loadmodel(modelfile):
    mf = open(modelfile)
    noOfmodes = int(mf.readline().strip().split(' ')[-1])
    features = int(mf.readline().strip().split(' ')[-1])
    totseqs = int(mf.readline().strip().split(' ')[-1])
    allpwms = {}
    l1 = mf.readline().strip()
    for i in range(noOfmodes):
        m = int(l1.split(' ')[-1])
        motif = []
        l1 = mf.readline().strip()
        while l1[0]!='M' and l1[0]!='P':
            l = map(int,l1.split('\t'))
            motif.append(l)
            l1 = mf.readline().strip()
        allpwms[m] = np.zeros((features,len(motif)),dtype=int)
        for j in range(len(motif)):
            allpwms[m][:,j] = motif[j]

    allpospwms = {}
    for i in range(noOfmodes):
        l1 = l1.split(' ')[-1]
        m = int(l1)
        pospwm = []
        l1 = mf.readline().strip()
        while l1[0]!='P' and l1[0]!='N':
            l = map(float,l1.split('\t'))
            pospwm.append(l)
            l1 = mf.readline().strip()
        allpospwms[m] = np.zeros((2,len(pospwm)),dtype=float)
        for j in range(len(pospwm)):
            allpospwms[m][:,j] = pospwm[j]
    
    allnegpwms = {}
    for i in range(noOfmodes):
        l1 = l1.split(' ')[-1]
        m = int(l1)
        negpwm = []
        l1 = mf.readline().strip()
        while l1[0]!='N' and l1[0]!='M':
            l = map(float,l1.split('\t'))
            negpwm.append(l)
            l1 = mf.readline().strip()
        allnegpwms[m] = np.zeros((2,len(negpwm)),dtype=float)
        for j in range(len(negpwm)):
            allnegpwms[m][:,j] = negpwm[j]
    modeSeqCounts = {}
    l1 = mf.readline().strip()
    while l1[0]!='M': # motif widths M is encountered
        l = map(int,l1.split('\t'))
        modeSeqCounts[l[0]] = l[1]
        l1 = mf.readline().strip()
    motifWidths= {}
    l1 = mf.readline().strip()
    while l1[0]!='P': # Positive strand widths is encountered
        l = map(int,l1.split('\t'))
        motifWidths[l[0]] = l[1]
        l1 = mf.readline().strip()

    posreadsWidths = {}
    l1 = mf.readline().strip()
    while l1[0]!='N': #Negative strand widths is encountered
        l = map(int,l1.split('\t'))
        posreadsWidths[l[0]] = l[1]
        l1 = mf.readline().strip()

    negreadsWidths = {}
    l1 = mf.readline().strip()
    while l1[0]!='R': # Reldist is encountered
        l = map(int,l1.split('\t'))
        negreadsWidths[l[0]] = l[1]
        l1 = mf.readline().strip()

    reldist = {}
    l1 = mf.readline().strip()
    while l1[0]!='A': # Alpha is encountered
        l = map(int,l1.split('\t'))
        reldist[l[0]] = (l[1],l[2]) # +ve and -ve reldist
        l1 = mf.readline().strip()

    alpha = float(l1.split(' ')[-1])
    l1 = mf.readline()
    pcReads = map(float,mf.readline().strip().split('\t'))
    mf.close()

    #model = [noOfmodes,features,totseqs,allpwms,allpospwms,allnegpwms,modeSeqCounts,motifWidths,posreadsWidths,negreadsWidths,reldist,alpha,pcReads]
    m = {}
    m['noOfmodes'] = noOfmodes
    m['features'] = features
    m['totseqs'] = totseqs
    m['allpwms'] = allpwms
    m['allpospwms'] = allpospwms
    m['allnegpwms'] = allnegpwms
    m['modeSeqCounts'] = modeSeqCounts
    m['motifWidths'] = motifWidths
    m['posreadsWidths'] = posreadsWidths
    m['negreadsWidths'] = negreadsWidths
    m['reldist'] = reldist
    m['alpha'] = alpha
    m['pcReads'] = pcReads
    return m
