#ifndef _traindata_h
#define _traindata_h

extern const float modesPrior;
extern const int defaultReadWinsize;
double getBackScore(int, int, double*);
double getBackScoreReadsMotif(int, int, double*, double*);
double calculateLikelihoodMode(model*, dataSet*, int*, int*, int*, int*, double**, double**, double**, int, float*);
void addRemoveDataPoint(model*, dataSet*, int*, int*, int*, int*, int, int);
double calculateLikelihood(model*, dataSet*, int*, int*, int*, int*, double**, double**, double**, float*);
double scoreMotif(model*, dataSet*, int, int, int, double*);
double scoreReads(model*, dataSet*, int, int, int, int, int, double*, double*);
double scoreMotifRegionReads(model*, dataSet*, int, int, int, double*, double*);
int sampleStartPosn(model*, dataSet*, int, int, unsigned int*, int, int, int, double*, double*, double*);
int sampleLabel(model*, dataSet*, int*, int, unsigned int*, int, double*, double*, double*,int);
void updateReadsStartPos(model*, int*, int*, int, int, int);
void takeExp(double*, int*, int);
int getMax(double*, int*, int);

int sampleMotifWidthRight(dataSet*, model*, int*, int*, int*, int, int, int, double**, double **, double**, int, unsigned int*,int);
int sampleMotifWidthLeft(dataSet*, model*, int*, int*, int*, int, int, int, double**, double**, double**, int, unsigned int*,int);
int sampleReadsWidthRight(dataSet*, model*, int*, int*, int*, int, double**, int, unsigned int*, int,int);
int sampleReadsWidthLeft(dataSet*, model*, int*, int*, int*,int, double**, int, unsigned int*, int,int);

trainOut *trainData(dataSet*, int, float, float*, unsigned int, double**, double**, double**, int*, int,int, int, int, int, int,char*,char *);

double updateBestModel(model*, dataSet*, int*, int*, int*, int*, double**, double**, double**, double, int,int,int, FILE*, int);
void printModelParams(model*, int, char*);

void freeTo(trainOut*);
/// Diversity initialization
void initializeDivLabelStart(dataSet*, int*, int*, char*);
#endif
