#ifndef _motifAndReadsStructs_h
#define _motifAndReadsStructs_h

typedef struct motifstruct{
  int *modeMotifCount;
  struct motifstruct *next;
}motifStruct;

typedef struct motifCont{
  motifStruct *motif;
}motifContainer;

typedef struct readparams{
  double *modeReadsCount;
  struct readparams *next;
}readParams;

typedef struct readCont{
  readParams *readsWin;
}readsContainer;

motifStruct* createNode();
motifStruct* initializeMotifPerMode(int);
motifContainer* initializeMotifs(int *, int);
readParams* createReadsNode();
readParams* initializeReadParamsPerMode(int);
readsContainer* initializeReadParams(int*,int);

readsContainer *initializeMotifReadsParams(int*,int);
void freeMotifs(motifContainer*,int);
void freeMotif(motifStruct*);
void freeReadsParams(readsContainer*, int);
void freeReadsWin(readParams*);
#endif
