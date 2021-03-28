#ifndef _motifAndReadsFunctions_h
#define _motifAndReadsFunctions_h

void motifDecreaseRight(dataSet*, model*, double**, int*, int*, int, int*, double**);
void motifIncreaseRight(dataSet*, model*, double**, int*, int*, int, int*, double**);
void motifDecreaseLeft(dataSet*, model*, double**, int*, int*, int, int*, double**);
void motifIncreaseLeft(dataSet*, model*, double**, int*, int*, int, int*, double**);
motifStruct *delNodeEnd(motifStruct*);
motifStruct *delNodeFront(motifStruct*);
readParams *delNodeEndReads(readParams*);
readParams *delNodeFrontReads(readParams*);

motifStruct *insertNodeEnd(motifStruct*, motifStruct*);
motifStruct *insertNodeFront(motifStruct*, motifStruct*);
readParams *insertReadsFront(readParams*, readParams*);
readParams *insertReadsEnd(readParams*, readParams*);

double scoreNochangeRight(dataSet*, model*, double**, int*, int*, int, int*, double**);
double scoreDecreaseRight(dataSet*, model*, double**, int*, int*, int, int*, double**,int);
double scoreIncreaseRight(dataSet*, model*, double**, int*, int*, int, int*, double**, int);

double scoreDecreaseLeft(dataSet* ,model*, double**, int*, int*, int, int*, double**,int);
double scoreNochangeLeft(dataSet* ,model*, double**, int*, int*, int, int*, double**);
double scoreIncreaseLeft(dataSet* ,model*, double**, int*, int*, int, int*, double**,int);

void readsDecreaseRight(model*, int, int*, double**, int);
void readsIncreaseRight(dataSet*, model*, int*, int, int*, int);
void readsDecreaseLeft(dataSet*, model*, int*, int, int*, double**, int);
void readsIncreaseLeft(dataSet*, model*, int*, int, int *,int);

double readScoreFirstnLast(dataSet*, model*, int*, int*, int, double**, int, int);
double readScoreBeforeFirst(dataSet*, model*, int*, int*, int, double**, int*, int,int,int);
double readScoreAfterLast(dataSet*, model*, int*, int*, int, double**, int*, int,int,int);

double readScoreMotif(dataSet*, model*, int, int*, int*, double **, double**, int);
double readScoreMotifExtend(dataSet*, model*, int, int*, int*, double**, double**, int);
double getReadsScoreMode(int, int, int*, int, int*, double**);
motifStruct *getLastNode(motifStruct*);
readParams *getLastReadNode(readParams*);

void fillReadsPWM(model*, dataSet*, int, int*, int*, int*, readParams*, int, int);


int samplePositiveOffset(model*, dataSet*, int, int*, int*, int*, int, int, int, double**, unsigned int*, int,int);
int sampleNegativeOffset(model*, dataSet*, int, int*, int*, int*, int, int, int, double**, unsigned int*, int,int);
#endif
