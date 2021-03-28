from ctypes import *

# create ctypes for C data structures dataSet and trainOut

class DataSet(Structure):
    _fields_ = [("data",POINTER(POINTER(c_int))),("features",POINTER(c_int)),("posreads",POINTER(POINTER(c_double))),("negreads",POINTER(POINTER(c_double))),("lookahead",POINTER(POINTER(c_int))),("featureValues",c_int),("n",c_int)]

class TrainOut(Structure):
    _fields_ = [("labels",POINTER(c_int)),("startPos",POINTER(c_int)),("preadsStart",POINTER(c_int)),("nreadsStart",POINTER(c_int)),
                ("motifWidth",POINTER(c_int)),("preadsWidth",POINTER(c_int)),("nreadsWidth",POINTER(c_int)),("likelihood",c_double)]
datalib = cdll.LoadLibrary("/home/anushua/ReadsAndMotifs/scripts/binary-version/withReads/exoDiv_trials/exoDiv_finalVersion4/dataStructures.so")
def getCvar(l,datatype):
    n = len(l)
    #arr = (datatype*n)(*l)
    arrOfN = datatype*n
    num = arrOfN()
    for i in xrange(n):
        num[i] = l[i]
    return num

def makeC2darray(arr,datatype):
    D=datatype
    PD=POINTER(D)
    #PPD=POINTER(PD)
    rows = len(arr)
    ARR_rows = rows*PD
    ptr = ARR_rows()
    for i in xrange(rows):
        # fill out each pointer with an array of D
        cols = len(arr[i])
        ARR_cols = (cols*D)
        ptr[i]=ARR_cols()
        for j in xrange(cols):
            ptr[i][j] = arr[i][j]
    return ptr

def getTrainOut(to,mode,n):
    d={}
    d['labels']=[to.contents.labels[i] for i in xrange(n)]
    d['startPos']=[to.contents.startPos[i] for i in xrange(n)]
    d['preadsStart']=[to.contents.preadsStart[i] for i in xrange(n)]
    d['nreadsStart']=[to.contents.nreadsStart[i] for i in xrange(n)]
    
    d['motifWidth']=[to.contents.motifWidth[j] for j in xrange(mode)]
    d['preadsWidth']=[to.contents.preadsWidth[j] for j in xrange(mode)]
    d['nreadsWidth']=[to.contents.nreadsWidth[j] for j in xrange(mode)]
    d['likelihood'] = to.contents.likelihood
    return d

datalib.trainData.restype = POINTER(TrainOut)
datalib.freeTo.restype= None

