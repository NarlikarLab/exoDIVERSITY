#ifndef _modelStructures_h
#define _modelStructures_h

typedef struct dataset{
  int **data;
  int *features;
  double **posreads;
  double **negreads;
  int **lookahead;
  int featureValues;
  int n;
}dataSet;

typedef struct trainout{
  int *labels;
  int *startPos;
  int *preadsStart;
  int *nreadsStart;
  int *motifWidth;
  int *preadsWidth;
  int *nreadsWidth;
  double likelihood;
}trainOut;

typedef struct reldist{
  int preadsMotif;
  int nreadsMotif;
} relDist;

typedef struct modelstruct{
  int mode;
  int featureValues;
  int n;
  motifContainer *motifs;
  readsContainer *posreadsparams;
  readsContainer *negreadsparams;
  int *counts;  // count of sequences in each mode
  int *mWidth;
  int *preadsWidth;
  int *nreadsWidth;
  relDist *readMotifDist;
  float alpha;
  float *pcReads;
}model;

void initializeLabelStartPos(dataSet*, int*, int*, int, int*, int*, int*, int, unsigned int*);
void initializeReadWinStartPos(int*, int*, int*, int*, int, int*, int*, int*, int*, int*);
model* createModel(int, dataSet*, int*, int*, int*, int*, float, float*, int*,int*, int*, int*, int*);
int motifWithN(int*, int, int);
model* initializeModel(int, int*, int*, int*);
relDist *initializeRelDist(int);
void getMotifCount(model*, dataSet*, int*, int*);
void getReadParams(model*, dataSet*, int*, int*, int*, int*);
void getMotifReadParams(model*, dataSet*, int*, int*);
void getRelativeDistance(model*,int*,int*);
void getReadsback(dataSet*, float*, int, double*, double*);
void copyarr(int*, int*,int);
void copyOffsets(int*, int*, model*);
int arrayMax(double*, int);
int sample(double*, int, double);
void assignTrainOut(trainOut*, int, int, int*, int*, int*, int*, int*,int*, int*, double);
void freeData(dataSet*);
void freeModel(model*);
void freeTo(trainOut*);
#endif

