#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include "motifAndReadsStructs.h"
#include "modelStructures.h"
#include "messages.h"

int motifWithN(int *seq, int start, int width){
  int i;
  for (i=start;i<(start+width);i++){
    if(seq[i]==4)
      return 1;
  } 
  return 0;
}

void initializeLabelStartPos(dataSet *ds, int* labels, int* startPos, int mode, int* mWidth, int *preadspart, int *nreadspart, int revFlag, unsigned int *seed){
  int i,count,ul,choice;
  for (i=0;i<ds->n;i++){
    labels[i]=rand_r(seed)%mode;
    if((preadspart[labels[i]]>0) && (nreadspart[labels[i]]>0)){
      if (((ds->features)[i] - (mWidth[labels[i]]+preadspart[labels[i]]+nreadspart[labels[i]]) + 1)<=0){
	startPos[i]=-1;
	labels[i]=-1;
	continue;
      }
    }
    else if(preadspart[labels[i]]<=0 && nreadspart[labels[i]]>0){
      if (((ds->features)[i] - (mWidth[labels[i]]+nreadspart[labels[i]]) + 1)<=0){
	startPos[i]=-1;
	labels[i]=-1;
	continue;
      }
    }
    else if(preadspart[labels[i]]>0 && nreadspart[labels[i]]<=0){
      if (((ds->features)[i] - (mWidth[labels[i]]+preadspart[labels[i]]) + 1)<=0){
	startPos[i]=-1;
	labels[i]=-1;
	continue;
      }
    }
    else if(preadspart[labels[i]]<=0 && nreadspart[labels[i]]<=0){
      if((ds->features)[i] - (mWidth)[labels[i]] + 1 <=0){
	startPos[i]=-1;
	labels[i]=-1;
	continue;
      }
    }
    
    count = 0;
    while(count<30){
      if(revFlag){
	choice = rand_r(seed)%2;
	if (preadspart[labels[i]]>0){
	  if(nreadspart[labels[i]]>0)
	    ul = ((ds->features)[i])/2 - (preadspart[labels[i]]+nreadspart[labels[i]]) - mWidth[labels[i]]+1;
	  else
	    ul = ((ds->features)[i])/2 - preadspart[labels[i]] - mWidth[labels[i]]+1;
	  if (choice){ //1: forward strand
	    startPos[i]=(rand_r(seed)%ul)+preadspart[labels[i]];
	  }
	  else{ //0: reverse strand
	    startPos[i]=(rand_r(seed)%ul)+((ds->features)[i]/2)+preadspart[labels[i]];
	  }
	}
	else {
	  if (nreadspart[labels[i]] > 0)
	    ul = ((ds->features)[i])/2 - nreadspart[labels[i]] - mWidth[labels[i]]+1;
	  else
	    ul = ((ds->features)[i])/2 - mWidth[labels[i]]+1;
	  if (choice){ //1: forward strand
	    startPos[i]=(rand_r(seed)%ul);
	  }
	  else{ //0: reverse strand
	    startPos[i]=(rand_r(seed)%ul)+((ds->features)[i]/2);
	  }
	}
      }
      else{
	if (preadspart[labels[i]]>0){
	  if (nreadspart[labels[i]]>0)
	    ul = (ds->features)[i]-(preadspart[labels[i]]+nreadspart[labels[i]])-mWidth[labels[i]]+1;
	  else
	    ul = (ds->features)[i]-preadspart[labels[i]]-mWidth[labels[i]]+1;
	  startPos[i]=(rand_r(seed)%ul)+preadspart[labels[i]];
	}
	else{
	  if (nreadspart[labels[i]]>0)
	    ul = (ds->features)[i]-nreadspart[labels[i]]-mWidth[labels[i]]+1;
	  else
	    ul = (ds->features)[i]-mWidth[labels[i]]+1;
	  startPos[i]=(rand_r(seed)%ul);
	}
      }
      if (motifWithN((ds->data)[i],startPos[i],mWidth[labels[i]]))
	count++;
      else
	break;
    } //while end
    if (count==30) {
      //printf("count became 30\n");
      startPos[i]=-1;
      labels[i]=-1;
    }
    
  }// for end
  /*Only Printing
  printf("Startpos\n");
  for (i=0;i<ds->n;i++)
    printf("%d\n",startPos[i]);
  printf("Labels\n");
  for (i=0;i<ds->n;i++)
    printf("%d\n",labels[i]);
  */
}

void initializeReadWinStartPos(int *posreadsStart, int *negreadsStart, int *startPos, int *labels, int n, int *mWidth, int *prWidth, int *nrWidth, int *preadspart, int *nreadspart){
  int i;
  for (i=0;i<n;i++){
    if(startPos[i]<0) {
      posreadsStart[i]=-1;
      negreadsStart[i]=-1;
    }
    else{
      posreadsStart[i] = startPos[i] - preadspart[labels[i]];
      negreadsStart[i] = startPos[i] + mWidth[labels[i]] + nreadspart[labels[i]] - nrWidth[labels[i]];
    }
  }
  //for(i=0;i<n;i++){
  //  printf("motstart: %d prstart: %d nrstart: %d\n",startPos[i],posreadsStart[i],negreadsStart[i]);
  //}
}

void getRelativeDistance(model *m, int *preadspart, int *nreadspart){
  int i;
  relDist *rl;
  rl = m->readMotifDist;
  for(i=0;i<(m->mode);i++){
    (rl[i].preadsMotif) = preadspart[i];
    (rl[i].nreadsMotif) = nreadspart[i] + (m->mWidth)[i];
  }
}

model* createModel(int mode, dataSet *ds, int *labels, int *startPos, int *posreadsStart, int *negreadsStart,float alpha, float *pcReads, int *mWidth, int *prWidth, int *nrWidth, int *preadspart, int *nreadspart){
  //Build PFM for each mode and gather the binary reads parameters
  model *m;
  int i;
  m = initializeModel(mode,mWidth,prWidth,nrWidth);
  m->mode = mode;
  m->featureValues = ds->featureValues;
  m->n = ds->n;
  m->alpha = alpha;
  for (i=0;i<2;i++){
    (m->pcReads)[i] = pcReads[i];
  }
  for (i=0;i<mode;i++){
    (m->mWidth)[i] = mWidth[i];
    (m->preadsWidth)[i] = prWidth[i];
    (m->nreadsWidth)[i] = nrWidth[i];
  }
  getMotifCount(m,ds,startPos,labels);
  getReadParams(m,ds,startPos,labels,posreadsStart,negreadsStart);
  getRelativeDistance(m,preadspart,nreadspart);
  return m;
}

relDist *initializeRelDist(int mode){
  int i;
  relDist *rmd;
  rmd = (relDist*)malloc(sizeof(relDist)*mode);
  for (i=0;i<mode;i++){
    rmd[i].preadsMotif = 0;
    rmd[i].nreadsMotif = 0;
  }
  return rmd;
}

model *initializeModel(int mode, int* mWidth, int *prWidth, int *nrWidth){

  model *m;
  m=(model*)malloc(sizeof(model));
  if(!m) printMessage(0);
  m -> motifs = initializeMotifs(mWidth,mode);
  m -> posreadsparams = initializeReadParams(prWidth,mode);
  m -> negreadsparams = initializeReadParams(nrWidth,mode);

  m -> readMotifDist = initializeRelDist(mode);
  m -> counts= (int*)malloc(sizeof(int)*mode);
  if (!(m->counts)) printMessage(0);
  m -> mWidth= (int*)malloc(sizeof(int)*mode);
  if(!(m->mWidth)) printMessage(0);
  m -> preadsWidth = (int*)malloc(sizeof(int)*mode);
  if(!(m->preadsWidth)) printMessage(0);
  m -> nreadsWidth = (int*)malloc(sizeof(int)*mode);
  if(!(m->nreadsWidth)) printMessage(0);
  m -> pcReads =(float*)malloc(sizeof(int)*2);
  if(!(m->pcReads)) printMessage(0);
  return m;
}

void getMotifCount(model *m,dataSet *ds,int *startPos, int *labels){
  motifStruct *m1;
  int i,j,k,mode;
  for (i=0;i<(m->mode);i++){
    m1 = (m->motifs)[i].motif;
    for(j=0;j<(m->mWidth)[i];j++){
      for(k=0;k<(ds->featureValues);k++){
	(m1->modeMotifCount)[k]=0;
      }
      m1=m1->next;
    }
    (m->counts)[i]=0;
  }
  m->n=0;
  for (i=0;i<(ds->n);i++){
    mode = labels[i];
    if (startPos[i] < 0 || labels[i] < 0) continue;
    (m->counts)[mode]++;
    m->n++;
    m1 = (m->motifs)[mode].motif;
    for (j=0;j<(m->mWidth)[mode];j++){
      (m1->modeMotifCount)[(ds->data)[i][startPos[i]+j]]++;
      m1 = m1->next;
    }
  }
  
  /* Printing */
  /*
  for(i=0;i< m->mode; i++){
    m1 = (m->motifs)[i].motif;
    printf("\nMotif PFM for mode %d\n",i);
    for(j=0;j<(m->mWidth)[i];j++){
      for(k=0;k<(ds->featureValues);k++){
	printf("%d\t",(m1->modeMotifCount)[k]);
      }
      m1=m1->next;
      printf("\n");
    }
    }*/
}

void getReadParams(model *m,dataSet *ds,int *startPos,int *labels, int *posreadsStart, int *negreadsStart){
  readParams *pr1;
  readParams *nr1;
  int i,j,k,mode;
  for (i=0;i<(m->mode);i++){
    pr1 = (m->posreadsparams)[i].readsWin;
    nr1 = (m->negreadsparams)[i].readsWin;
    for (j=0;j<(m->preadsWidth)[i];j++){
      for (k=0;k<2;k++){
	(pr1->modeReadsCount)[k]=0;
      }
      pr1=pr1->next;
    }
    for (j=0;j<(m->nreadsWidth)[i];j++){
      for(k=0;k<2;k++){
	(nr1->modeReadsCount)[k]=0;
      }
      nr1=nr1->next;
    }
  }
  for (i=0;i<(ds->n);i++){
    mode = labels[i];
    if (startPos[i]<0 || labels[i]<0) continue;
    pr1 = (m->posreadsparams)[mode].readsWin;
    nr1 = (m->negreadsparams)[mode].readsWin;
    for (j=0;j<(m->preadsWidth)[mode];j++){
      (pr1->modeReadsCount)[(int)((ds->posreads)[i][posreadsStart[i]+j])]++;
      pr1 = pr1->next;
    }
    for(j=0;j<(m->nreadsWidth)[mode];j++){
      (nr1->modeReadsCount)[(int)((ds->negreads)[i][negreadsStart[i]+j])]++;
      nr1 = nr1->next;
    }
  }
  /* Printing
  for(i=0;i< m->mode; i++){
    pr1 = (m->posreadsparams)[i].readsWin;
    printf("\nReads FM for mode with width: %d\n",i,(m->preadsWidth)[i]);
    for(j=0;j<(m->preadsWidth)[i];j++){
      for(k=0;k<2;k++){
	printf("%lf\t",(pr1->modeReadsCount)[k]);
      }
      pr1=pr1->next;
      printf("\n");
    }
  }
  */
}

void copyarr(int* source, int* dest, int n){
  int i;
  for (i=0;i<n;i++){
    dest[i]=source[i];
  }
}

void copyOffsets(int *preadspart, int *nreadspart, model *m){
  int i;
  relDist *rl;
  rl = m->readMotifDist;
  for (i=0;i<(m->mode);i++){
    preadspart[i] = (rl[i].preadsMotif);
    nreadspart[i] = (rl[i].nreadsMotif)-(m->mWidth)[i];
  }
}

int arrayMax(double *values,int n){
  int i,index=0;
  double tempmax;
  tempmax=values[0];
  for(i=1;i<n;i++){
    if(tempmax > values[i]){
      tempmax = values[i] ;
      index = i;
    }
  }
  return index;
}

int intervalBinsearch(int low, int high, double num, double *arr){
  int mid,count=0,orighigh,i;
  orighigh = high;
  if(num<arr[0])
    return 0;
  while(low<=high){
    if (count > orighigh) {
      printf("Stuck in loop with orighigh: %d\n",orighigh);
      for (i=0;i<orighigh;i++) printf("%lf \t",arr[i]);
      printf("\n");
      exit(0);
    }
    mid = (low+high)/2;
    count++;
    if(num>arr[mid-1] && num<=arr[mid])
      return mid;
    else if(num < arr[mid]){
      high=mid-1;
    }
    else if(num > arr[mid]){
      low = mid+1;
    }
  }
  return -1;
}
int sample(double *p, int n, double r){
  int i,index;
  double sum;
  sum=0;
  for (i=0;i<n;i++){
    sum=sum+p[i];
  }
  //Normalizing
  for (i=0;i<n;i++){
    p[i]=p[i]/sum;
  }
  //Cumulative  
  for (i=1;i<n;i++){
    p[i]=p[i]+p[i-1];
  }
  index = intervalBinsearch(1,(n-1),r,p);
  return index;
}

void getReadsback(dataSet *ds, float *pcReads, int mode, double *posreadsback, double *negreadsback){
  int i,j;
  double total=0.0;
  for (i=0;i<2;i++){
    posreadsback[i]=0.0;
    negreadsback[i]=0.0;
  }
  for (i=0;i<(ds->n);i++){
    for (j=0;j<(ds->features)[i];j++){
      posreadsback[(int)(ds->posreads)[i][j]]++;
      negreadsback[(int)(ds->negreads)[i][j]]++;
    }
    total+=(ds->features)[i];
  }
  total+= (mode*pcReads[0] + mode*pcReads[1]);
  for (i=0;i<2;i++){
    posreadsback[i] = (posreadsback[i]+mode*pcReads[i])/total;
    negreadsback[i] = (negreadsback[i]+mode*pcReads[i])/total;
  }
}

void assignTrainOut(trainOut *to, int n, int mode, int *labels, int *startPos, int *posreadsStart, int *negreadsStart, int *widths,int *preadsWidth, int *nreadsWidth, double maxLikelihood){
  int i;
  
  to->labels = (int*)malloc(sizeof(int)*n);
  if(!(to->labels)) printMessage(0);
  for (i=0;i<n;i++) (to->labels)[i] = labels[i];

  to->startPos = (int*)malloc(sizeof(int)*n);
  if(!(to->startPos)) printMessage(0);
  for(i=0;i<n;i++) (to->startPos)[i] = startPos[i];

  to->preadsStart = (int*)malloc(sizeof(int)*n);
  if(!(to->preadsStart)) printMessage(0);
  for(i=0;i<n;i++)  (to->preadsStart)[i] = posreadsStart[i];

  to->nreadsStart = (int*)malloc(sizeof(int)*n);
  if(!(to->nreadsStart)) printMessage(0);
  for(i=0;i<n;i++) (to->nreadsStart)[i] = negreadsStart[i];

  to->motifWidth = (int*)malloc(sizeof(int)*mode);
  if(!(to->motifWidth)) printMessage(0);
  for(i=0;i<mode;i++) (to->motifWidth)[i] = widths[i];

  to->preadsWidth = (int*)malloc(sizeof(int)*mode);
  if(!(to->preadsWidth)) printMessage(0);
  for(i=0;i<mode;i++) (to->preadsWidth)[i] = preadsWidth[i];

  to->nreadsWidth = (int*)malloc(sizeof(int)*mode);
  if(!(to->nreadsWidth)) printMessage(0);
  for(i=0;i<mode;i++) (to->nreadsWidth)[i] = nreadsWidth[i];
  to->likelihood = maxLikelihood;
}


void freeData(dataSet *ds){
  int i;
  for (i=0;i<(ds->n);i++){
    free((ds->data)[i]); (ds->data[i])=NULL;
    free((ds->lookahead)[i]); (ds->lookahead)[i]=NULL;
    free((ds->posreads)[i]); (ds->posreads)[i]=NULL;
    free((ds->negreads)[i]); (ds->negreads)[i]=NULL;
  }
  free(ds->features); ds->features = NULL;
  free(ds->data);   ds->data = NULL;
  free(ds->lookahead); ds->lookahead = NULL;
  free(ds->posreads);  ds->posreads = NULL;
  free(ds->negreads);  ds->negreads = NULL;
  free(ds);
}

void freeModel(model *m){
  freeMotifs(m->motifs,m->mode);
  freeReadsParams(m->posreadsparams,m->mode);
  freeReadsParams(m->negreadsparams,m->mode);
  free(m->readMotifDist); m->readMotifDist = NULL;
  free(m->counts); m->counts=NULL;
  free(m->mWidth); m->mWidth=NULL;
  free(m->preadsWidth); m->preadsWidth=NULL;
  free(m->nreadsWidth); m->nreadsWidth=NULL;
  free(m->pcReads); m->pcReads=NULL;
  free(m);
}

