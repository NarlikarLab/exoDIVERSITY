#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<limits.h>
#include "motifAndReadsStructs.h"
#include "modelStructures.h"
#include "motifAndReadsFunctions.h"
#include "traindata.h"
#include "messages.h"

motifStruct *createNode(){
  motifStruct *m = (motifStruct*)malloc(sizeof(motifStruct));
  if (!m) printMessage(0);
  m -> modeMotifCount = (int*)malloc(sizeof(int)*4);
  if(!(m->modeMotifCount)) printMessage(0);
  m->next=NULL;
  return m;
}

void freeNode(motifStruct *m){
  free(m->modeMotifCount);
  m->next = NULL;
  m->modeMotifCount = NULL;
  free(m);
}

void freeNodeReads(readParams *r){
  free(r->modeReadsCount);
  r->modeReadsCount = NULL;
  r->next = NULL;
  free(r);
}

motifStruct *delNodeEnd(motifStruct *m){
  motifStruct *mc, *mc1;
  mc = m;
  while (((mc->next)->next)!=NULL) mc = mc->next;
  mc1 = mc->next;
  freeNode(mc1);
  mc->next = NULL;
  return m;
}

motifStruct *delNodeFront(motifStruct *m){
  motifStruct *mc;
  mc = m->next;
  freeNode(m);
  return mc;
}

readParams *delNodeEndReads(readParams *r){
  readParams *rc, *rc1;
  rc = r;
  while (((rc->next)->next)!=NULL) rc = rc->next;
  rc1 = rc->next;
  freeNodeReads(rc1);
  rc->next=NULL;
  return r;
}

readParams *delNodeFrontReads(readParams *r){
  readParams *rc;
  rc = r->next;
  freeNodeReads(r);
  return rc;
}

readParams *insertReadsFront(readParams *f, readParams *r){
  f->next = r;
  return f;
}

readParams *insertReadsEnd(readParams *e, readParams*r){
  readParams *r1;
  r1 = r;
  while((r1->next)!=NULL) r1 = r1->next;
  r1->next = e;
  return r;
}

motifStruct *insertNodeEnd(motifStruct *e,motifStruct *m){
  motifStruct *mc;
  mc=m;
  while((mc->next)!=NULL) mc=mc->next;
  mc->next = e;
  return m;
}

motifStruct *insertNodeFront(motifStruct *f, motifStruct *m){
  f->next = m;
  return f;
}

motifStruct *initializeMotifPerMode(int width){
  int i;
  motifStruct *m, *mc;
  m = createNode();
  mc = m;
  for (i=1;i<width;i++){
    mc->next = createNode();
    mc = mc->next;
  }
  return m;
}

motifContainer *initializeMotifs(int *mWidth,int mode){
  motifContainer *mc;
  int i; 
  mc = (motifContainer*)malloc(sizeof(motifContainer)*mode);
  if (!mc) printMessage(0);
  for (i=0;i<mode;i++){
    mc[i].motif = initializeMotifPerMode(mWidth[i]);
  }
  return mc;
}

readParams *createReadsNode(){
  readParams *r=(readParams*)malloc(sizeof(readParams));
  if (!r) printMessage(0);
  r -> modeReadsCount = (double *)malloc(sizeof(double)*2);
  if(!(r->modeReadsCount)) printMessage(0);
  r->next=NULL;
  return r;
}

readParams *initializeReadParamsPerMode(int readWinsize){
  int i;
  readParams *r,*rc;
  r = createReadsNode();
  rc = r;
  for (i=1;i<readWinsize;i++){
    rc->next = createReadsNode();
    rc = rc->next;
  }
  return r;
}
readsContainer *initializeReadParams(int *readWinsize,int mode){
  readsContainer *rc;
  int i;
  rc = (readsContainer *)malloc(sizeof(readsContainer)*mode);
  if(!rc) printMessage(0);
  for(i=0; i<mode; i++){

    rc[i].readsWin = initializeReadParamsPerMode(readWinsize[i]);
  }
  return rc;
}

readsContainer *initializeMotifReadsParams(int *mWidth,int mode){
  readsContainer *rc;
  int i;
  rc = (readsContainer*)malloc(sizeof(readsContainer)*mode);
  if(!rc) printMessage(0);
  for(i=0;i<mode;i++){
    rc[i].readsWin = initializeReadParamsPerMode(mWidth[i]);
  }
  return rc;
}

void freeMotif(motifStruct* m){
  motifStruct *mt;
  while(m!=NULL){
    mt = m;
    m=m->next;
    freeNode(mt);
  }
}

void freeMotifs(motifContainer *mc,int mode){
  int i;
  for (i=0;i<mode;i++){
    freeMotif(mc[i].motif);
    mc[i].motif=NULL;
  }
  free(mc);
}

void freeReadsWin(readParams *r){
  readParams *rt;
  while (r!=NULL){
    rt=r;
    r=r->next;
    freeNodeReads(rt);
  }
}
/*Free the start pointer of the linkedlist of every mode*/
void freeReadsParams(readsContainer *rc, int mode){
  int i;
  for(i=0;i<mode;i++){    
    freeReadsWin(rc[i].readsWin);
    rc[i].readsWin=NULL;
  }
  free(rc);
}

motifStruct *getLastNode(motifStruct *m){
  motifStruct *mc;
  mc = m;
  while (mc->next!=NULL)
    mc = mc->next;
  return mc;
}

readParams *getLastReadNode(readParams *r){
  readParams *rc;
  rc = r;
  while(rc->next!=NULL)
    rc = rc->next;
  return rc;
}

double getReadsScoreMode(int n, int mode, int *labels, int posn, int *readsStart, double **readsback){
  int i;
  double rb=0.0;
  for (i=0;i<n;i++){
    if (labels[i]!=mode) continue;
    rb = rb + log(readsback[i][readsStart[i]+posn]);
  }
  return rb;
}

double readScoreFirstnLast(dataSet *ds, model *m, int *labels, int *readsStart, int mode, double **readsback,int fl,int pn){
  int i;
  double r,rb,rtot;
  readParams *rn;
  r=0.0;
  rb=0.0;
  rtot=0.0;

  if (pn){ //pos
    if((m->preadsWidth)[mode] == 0) return rtot;
    if(fl){
      //printf("Positive strand reads first\n");
      rn = (m->posreadsparams)[mode].readsWin;
      rb = getReadsScoreMode(ds->n,mode,labels,0,readsStart,readsback);
    }
    else{
      //printf("Positive strand reads last\n");
      rn = getLastReadNode((m->posreadsparams)[mode].readsWin);
      rb = getReadsScoreMode(ds->n,mode,labels,(m->preadsWidth)[mode]-1,readsStart,readsback);
    }
  }
  else{ //neg
    if((m->nreadsWidth)[mode] == 0) return rtot;
    if(fl){
      //printf("Negative strand reads first\n");
      rn = (m->negreadsparams)[mode].readsWin; 
      rb = getReadsScoreMode(ds->n,mode,labels,0,readsStart,readsback);
    }
    else{
      //printf("Negative strand reads last\n");
      rn = getLastReadNode((m->negreadsparams)[mode].readsWin);
      rb = getReadsScoreMode(ds->n,mode,labels,(m->nreadsWidth)[mode]-1,readsStart,readsback);
    }
  }
  for (i=0;i<2;i++){
    //printf("%lf\t",(rn->modeReadsCount)[i]);
    rtot += (rn->modeReadsCount)[i]+(m->pcReads)[i];
    r = r + ((rn->modeReadsCount)[i]+(m->pcReads)[i])*log((rn->modeReadsCount)[i]+(m->pcReads)[i]);
    }
  //printf("\nr: %lf \t rtot: %lf\trb:%lf\n",r,rtot,rb);
  r = r-rtot*log(rtot)-rb;
  return r;
}

double readScoreBeforeFirst(dataSet *ds, model *m,int *labels, int *startPos,int mode, double **readsback, int *readsStart,int revFlag,int pn,int mr){
  double rb,r,rtot;
  double *rr0;
  int i;
  r = 0.0;
  rr0 = (double*)malloc(sizeof(double)*2);
  rr0[0]=0.0;
  rr0[1]=0.0;
  if (pn){
    if ((m->preadsWidth)[mode]==0) {
      free(rr0);
      return r;
    }
    for (i=0;i<ds->n;i++){
      if (labels[i]!=mode) continue;
      if (revFlag && readsStart[i]==((ds->features)[i])/2) {
	free(rr0);
	return r;
      }
      if (readsStart[i]==0){
	free(rr0);
	return r;
      }
      else{
	rr0[(int)(ds->posreads)[i][readsStart[i]-1]]+=1;
	rb = rb + log(readsback[i][readsStart[i]-1]);
      }
    }
  }
  else{
    if ((m->nreadsWidth)[mode]==0){
      free(rr0);
      return r;
    }
    for (i=0;i<ds->n;i++){
      if (labels[i]!=mode) continue;
      if (!mr && (readsStart[i]==(startPos[i]+(m->mWidth)[mode]))){
	free(rr0);
	return r;
      }
      else{
	rr0[(int)(ds->negreads)[i][readsStart[i]-1]]+=1;
	rb = rb + log(readsback[i][readsStart[i]-1]);
      }
    }
  }
  rtot= 0.0;
  r = 0.0;
  for (i=0;i<2;i++){
    rtot += rr0[i]+(m->pcReads)[i];
    r = r+(rr0[i]+(m->pcReads)[i])*log(rr0[i]+(m->pcReads)[i]);
  }
  r = r-rtot*log(rtot)-rb;
  free(rr0);
  return r;
}
double readScoreAfterLast(dataSet *ds,model *m,int *labels, int *startPos,int mode,double **readsback,int *readsStart,int revFlag,int pn, int mr){
  int i;
  double *rrlastnext,r,rb,rtot;
  rrlastnext = (double *)malloc(sizeof(double)*2);
  rrlastnext[0]=0.0;
  rrlastnext[1]=0.0;
  r=0.0;
  rb=0.0;
  rtot=0.0;
  if (pn){
    if ((m->preadsWidth)[mode]==0){
      free(rrlastnext);
      return r;
    }
    for (i=0;i<ds->n;i++){
      if(labels[i]!=mode) continue;
      if(!mr && (readsStart[i]+(m->preadsWidth)[mode] >= startPos[i])){
	free(rrlastnext);
	return r;
      }
      rrlastnext[(int)(ds->posreads)[i][readsStart[i]+(m->preadsWidth)[mode]]]+=1;
      rb += log(readsback[i][readsStart[i]+(m->preadsWidth)[mode]]);
    }
  }
  else{
    if ((m->nreadsWidth)[mode]==0){
      free(rrlastnext);
      return r;
    }
    for (i=0;i<ds->n;i++){
      if (labels[i]!=mode) continue;
      if ((readsStart[i]+(m->nreadsWidth)[mode])>=((ds->features)[i])){
	free(rrlastnext);
	return r;
      }
      if (revFlag && readsStart[i]<((ds->features)[i])/2 && (readsStart[i]+(m->nreadsWidth)[mode])>=((ds->features)[i])/2){
	free(rrlastnext);
	return r;
      }
      rrlastnext[(int)(ds->negreads)[i][readsStart[i]+(m->nreadsWidth)[mode]]]+=1;
      rb += log(readsback[i][readsStart[i]+(m->nreadsWidth)[mode]]);
    }
  }
  rtot=0.0;
  r=0.0;
  for (i=0;i<2;i++){
    rtot += rrlastnext[i]+(m->pcReads)[i];
    r=r+(rrlastnext[i]+(m->pcReads)[i])*log(rrlastnext[i]+(m->pcReads)[i]);
  }
  r=r-rtot*log(rtot)-rb;
  free(rrlastnext);
  return r;
}

double scoreNochangeRight(dataSet *ds, model *m, double **background, int *labels, int *startPos, int mode, int *nreadsStart, double **negreadsback){
  double s,sb,tot;
  motifStruct *ml;
  int i;
  ml = getLastNode((m->motifs)[mode].motif);
  sb = 0.0;
  for (i=0;i<ds->n;i++){
    if (labels[i]!=mode) continue;
    if (startPos[i]<0) continue;
    if (background[i][startPos[i] + (m->mWidth)[mode]-1] < 0.00001) continue;
    sb = sb+log(background[i][startPos[i] + (m->mWidth)[mode] -1]);
  }
  s=0.0;
  tot = 0.0;
  for (i=0;i<(m->featureValues);i++){
    tot = tot + (ml->modeMotifCount)[i]+(m->alpha);
    s = s+((ml->modeMotifCount)[i]+(m->alpha))*log((ml->modeMotifCount)[i]+(m->alpha));
  }
  //printf("No change right: s: %lf \tsb: %lf \ttot: %lf \n",s,sb,tot);
  s = s - tot*log(tot);
  s = s - sb;
  //printf("No change Motif %lf \t rr1: %lf \t rrlast: %lf\n",s,rr1,rrlast);
  return s;
}

double scoreDecreaseRight(dataSet *ds, model *m, double **background, int *labels, int *startPos, int mode, int *nreadsStart, double **negreadsback,int revFlag){
  // Calculate the score based only on reads
  double rr1,rr0;
  rr1 = readScoreFirstnLast(ds,m,labels,nreadsStart,mode,negreadsback,1,0); //first and neg
  rr0 = readScoreBeforeFirst(ds,m,labels,startPos,mode,negreadsback,nreadsStart,revFlag,0,1); // 0 for neg 1 is for motif sampling
  //printf("For decrease right\nrr0: %lf \t rr1: %lf\n",rr0,rr1);
  return (rr0+rr1);
}

double scoreIncreaseRight(dataSet *ds, model *m, double **background, int *labels,int *startPos, int mode, int *readsStart, double **readsback,int revFlag){
  int i;
  double *values,s,sb,stot;
  motifStruct *ml;
  sb=0.0;
  ml = getLastNode((m->motifs)[mode].motif);
  //printf("Increase on Right\n");
  values = (double *)malloc(sizeof(double)*(m->featureValues));
  for (i=0;i<(m->featureValues);i++) values[i]=0.0;
  for(i=0;i<ds->n;i++){
    if (labels[i]!=mode) continue;
    if (startPos[i]<0) continue;
    if ((ds->data)[i][startPos[i]+(m->mWidth)[mode]] > 3) {
      //printf("Got an N\n");
      free(values);
      return 0;
    }
    if (startPos[i]+(m->mWidth)[mode]+1 > (ds->features)[i]){
      free(values);
      return 0;
    }
    if(revFlag && (startPos[i]<(ds->features)[i]/2) && (startPos[i]+(m->mWidth)[mode]+1 > (ds->features)[i]/2)){
      free(values);
      return 0;
    }
    //if (revFlag && (startPos[i]<((ds->features)[i]/2)) && (startPos[i] + ((m->readMotifDist)[mode].nreadsMotif)+1 > ((ds->features)[i]/2))){
      //printf("Going beyond half sequence length. start pos: %d\n",startPos[i]);
      //free(values);
      //return 0;
    //}
    values[(ds->data)[i][startPos[i]+(m->mWidth)[mode]]]++;
    sb = sb + (log(background[i][startPos[i]+(m->mWidth)[mode]])+(log(background[i][startPos[i]+(m->mWidth)[mode]-1])));
  }
  s=0.0;
  stot=0.0;
  for (i=0;i<(m->featureValues);i++){
    stot += (values[i]+0.5);
    s = s+(values[i]+0.5)*log(values[i]+0.5)+((ml->modeMotifCount)[i]+0.5)*log((ml->modeMotifCount)[i]+0.5); // lastnext and last positions
  }
  s = s - (2*stot*log(stot)) - sb;
  //RR(nrw-1) RR(nrw)
  //printf("For Increase right: motif: %lf \t rrlast: %lf \t rrlastnext: %lf\n",s,rrlast,rrlastnext);
  free(values);
  return s;
}

void motifDecreaseRight(dataSet *ds, model *m, double **background, int *labels, int *startPos, int mode, int *nreadsStart, double **negreadsback){
  (m->motifs)[mode].motif=delNodeEnd((m->motifs)[mode].motif);
  (m->mWidth)[mode]--;
}

void motifIncreaseRight(dataSet *ds, model *m,double **background, int *labels, int *startPos, int mode, int *nreadsStart, double **negreadsback){
  int i;
  motifStruct *me;
  me = createNode();
  for (i=0;i<(m->featureValues);i++) (me->modeMotifCount)[i]=0;
  for (i=0;i<ds->n;i++){
    if (labels[i]!=mode) continue;
    if (startPos[i]<0) continue;
    if ((ds->data)[i][startPos[i]+(m->mWidth)[mode]]>3) {
      freeNode(me);
      return;
    }
    (me->modeMotifCount)[(ds->data)[i][startPos[i]+(m->mWidth)[mode]]]++;
  }
  (m->motifs)[mode].motif = insertNodeEnd(me,(m->motifs)[mode].motif);
  (m->mWidth)[mode]++;
}

double scoreDecreaseLeft(dataSet *ds, model *m, double **background, int *labels, int *startPos, int mode, int *preadsStart, double **posreadsback,int revFlag){
  double rrlast, rrlastnext;
  rrlast = readScoreFirstnLast(ds,m,labels,preadsStart,mode,posreadsback,0,1);  //last and pos;
  rrlastnext = readScoreAfterLast(ds,m,labels,startPos,mode,posreadsback,preadsStart,revFlag,1,1); // 1 is for positive and 1 is for motif
  //printf("Decrease Left: rrlast: %lf\trrlastnext: %lf\n",rrlast,rrlastnext);
  return (rrlast+rrlastnext);
}

double scoreIncreaseLeft(dataSet *ds, model *m, double **background,int *labels, int *startPos, int mode, int *readsStart, double **readsback,int revFlag){
  //Based on motif and reads both: {ms-1,ms}, {rs, rs-1}
  int i;
  double *values,s,sb,stot;
  motifStruct *ml;
  sb=0.0;
  stot=0.0;
  ml = (m->motifs)[mode].motif;
  values = (double*)malloc(sizeof(double)*(m->featureValues));
  //printf("Motif increase left\n");
  for (i=0;i<(m->featureValues);i++) values[i]=0.0;
  for (i=0;i<ds->n;i++){
    if (labels[i]!=mode) continue;
    if (startPos[i] == 0) {
      free(values);
      return 0;
    }
    if (startPos[i] <0) continue;
    
    if ((ds->data)[i][startPos[i]-1] > 3){
      //printf("Got an N at pos: %d\n",startPos[i]-1);
      free(values);
      return 0;
    }
    if (revFlag && (startPos[i] == (ds->features)[i]/2)){
      free(values);
      return 0;
    }
    values[(ds->data)[i][startPos[i]-1]]++;
    sb = sb + (log(background[i][startPos[i]-1]))+ (log(background[i][startPos[i]]));
  }
  s = 0.0;
  stot = 0.0;
  for (i=0;i<(m->featureValues);i++){
    stot += (values[i]+0.5);
    s = s+ (values[i]+0.5)*log(values[i]+0.5)+((ml->modeMotifCount)[i] +0.5)*log((ml->modeMotifCount)[i]+0.5);
  }
  s = s - (2*stot*log(stot)) - sb;
  free(values);
  //printf("Increase left: s:%lf \t rr0:%lf \t rr1:%lf\n",s,rr0,rr1);
  return s;
}

double scoreNochangeLeft(dataSet *ds, model *m, double **background, int *labels,int *startPos, int mode, int *readsStart, double **readsback){
  double s,sb,tot;
  motifStruct *ml;
  int i;
  ml = (m->motifs)[mode].motif;
  sb = 0.0;
  for (i=0;i<(ds->n);i++){
    if (labels[i]!=mode) continue;
    if(startPos[i]<0) continue;
    if (background[i][startPos[i]] < 0.00001) continue;
    sb = sb+log(background[i][startPos[i]]);
  }
  s = 0.0;
  tot = 0.0;
  for (i=0;i<(m->featureValues);i++){
    tot = tot + (ml->modeMotifCount)[i]+(m->alpha);
    s = s+((ml->modeMotifCount)[i]+(m->alpha))*log((ml->modeMotifCount)[i]+(m->alpha));
  }
  //printf("No change left: s:%lf\t sb: %lf tot:%lf\n",s,sb,tot);
  s = s - tot*log(tot) - sb;
  //printf("For no change left:  motif:%lf \t rr1: %lf \t rrlast: %lf \n",s,rr1,rrlast);
  return s;
}

void motifDecreaseLeft(dataSet *ds, model *m, double **background, int *labels, int *startPos, int mode, int *preadsStart, double **posreadsback){
  int i;
  /* readParams *r; */
  /* r = createReadsNode(); */
  /* for (i=0;i<2;i++) (r->modeReadsCount)[i]=0.0; */
  for (i=0;i<(ds->n);i++){
    if (labels[i]!=mode) continue;
    if (startPos[i] < 0) continue;
    //(r->modeReadsCount)[(int)(ds->posreads)[i][preadsStart[i] + (m->preadsWidth)[mode]]]+=1;
    //preadsStart[i] += 1;
    startPos[i] += 1;
  }
  (m->motifs)[mode].motif = delNodeFront((m->motifs)[mode].motif);
  (m->mWidth)[mode]--;
  (m->readMotifDist)[mode].preadsMotif +=1;
  (m->readMotifDist)[mode].nreadsMotif -=1;
}

void motifIncreaseLeft(dataSet *ds, model *m, double **background, int *labels,int *startPos,int mode,int *preadsStart, double **posreadsback){
  int i;
  motifStruct *ms;
  ms = createNode();
  for (i=0;i<(m->featureValues);i++) (ms->modeMotifCount)[i]=0;
  for (i=0;i<ds->n;i++){
    if (labels[i]!=mode) continue;
    if (startPos[i]<=0) continue;
    (ms->modeMotifCount)[(ds->data)[i][startPos[i]-1]]++;
    startPos[i] = startPos[i]-1;
  }
  (m->motifs)[mode].motif = insertNodeFront(ms,(m->motifs)[mode].motif);
  (m->mWidth)[mode]++;
  (m->readMotifDist)[mode].preadsMotif -=1; 
  (m->readMotifDist)[mode].nreadsMotif +=1;
}

void fillReadsPWM(model *m, dataSet *ds, int mode,int *startPos, int *labels, int *readsStart, readParams *r, int width, int pn){
  int i,j;
  readParams *r1;
  r1 = r;
  for(i=0;i<width;i++){
    for (j=0;j<2;j++){
      (r1->modeReadsCount)[j]=0;
    }
    r1 = r1->next;
  }
  
  for (i=0;i<(ds->n);i++){
    if (labels[i]!=mode) continue;
    if (startPos[i]<0 || labels[i]<0) continue;
    r1 = r;
    for (j=0;j<width;j++){
      if (!pn){
	(r1->modeReadsCount)[(int)((ds->posreads)[i][readsStart[i]+j])]++;
	r1=r1->next;
      }
      else{
	(r1->modeReadsCount)[(int)((ds->negreads)[i][readsStart[i]+j])]++;
	r1=r1->next;
      }
    }
  }
}


int samplePositiveOffset(model *m, dataSet *ds, int mode, int *labels, int *startPos, int *preadsStart, int prWidth, int preadspart, int maxoffset, double **posreadsback,unsigned int *seed, int rev,int gobeyond){
  //preadsStart can change at the end
  int maxleftwidth,maxrightwidth,maxallowed;
  int i,j,k,flag=0,ind,posnBeforeMotif,start;
  readsContainer *tr;
  readParams *tv;
  int *tpreadsStart,*valid;
  double *poslikes,maxlike=0;
  double rp,rpt,s,totpc,sum=0;
  //printf("From pos offset\n");
  maxrightwidth = 0;
  maxleftwidth = maxoffset + prWidth;
  //if(2*prWidth < (m->mWidth)[mode]) start=2*prWidth-(m->mWidth)[mode];
  if (prWidth > (m->mWidth)[mode]) start = prWidth - (m->mWidth)[mode];
  else start = 0; 
  for (i=start;i<maxleftwidth;i++){
    flag=0;
    for (j=0;j<ds->n;j++){
      if(labels[j]!=mode) continue;
      if(startPos[j]<0) continue;
      if(startPos[j]-i<0) {
	flag=1;
	break;
      }
      if(rev){
	  // When motif start is in reverse strand and we cannot extend positive read window as it is getting into the forward strand
	if((startPos[j]>=(ds->features)[j]/2) && (startPos[j]-i <(ds->features)[j]/2)){
	  flag = 1;
	  break;
	}
      }
       
    }
    if (flag) break;
  }
  if (prWidth > (m->mWidth)[mode] && i==start) return (i-1);
  if(i >0) maxallowed = i-1;
  else maxallowed = 0;
  
  if (maxallowed < 0){
    printf("Positive maxallowed: %d\n",maxallowed);
    exit(0);
  }

  if (gobeyond){
    flag = 0;
    for (i=0;i<(m->mWidth)[mode];i++){
      for (j = 0;j<ds->n;j++){
	if(labels[j]!=mode) continue;
	if(startPos[j]<0) continue;
	if (startPos[j]+i+prWidth>(ds->features)[j]){
	  flag=1;
	  break;
	}
	if(rev && startPos[j]<(ds->features)[j]/2 && startPos[j]+i+prWidth>(ds->features)[j]/2){
	  flag=1;
	  break;
	}
      }
      if (flag) break;
    }
    maxrightwidth = i-1;
  }
  else{
    maxrightwidth = (m->mWidth)[mode]-prWidth;
  }
  //poslikes = (double*)malloc(sizeof(double)*(maxallowed+((m->mWidth)[mode]-prWidth)+1));
  //valid = (int*)malloc(sizeof(int)*(maxallowed+((m->mWidth)[mode]-prWidth)+1));
  poslikes = (double*)malloc(sizeof(double)*(maxallowed+maxrightwidth+1));
  valid = (int*)malloc(sizeof(int)*(maxallowed+maxrightwidth+1));
  tpreadsStart = (int*)malloc(sizeof(int)*(ds->n));
  totpc = (m->pcReads)[0]+(m->pcReads)[1];
  tr = (readsContainer*)malloc(sizeof(readsContainer)); 
  tr->readsWin = initializeReadParamsPerMode(prWidth);
  maxlike = -1.0*INT_MAX;

  for(i=maxallowed;i>=-maxrightwidth;i--){
    /* adjust read start positions as we slide*/
    for (j=0;j<ds->n;j++){
      if (labels[j]!=mode) continue;
      if (startPos[j]<0) continue;
      tpreadsStart[j] = startPos[j] - i; 
    }
    /* initialize the pos reads matrix  */
    tv = tr->readsWin;
    for (j=0;j<prWidth;j++){
      for (k=0;k<2;k++)
	(tv->modeReadsCount)[k]=0;
      tv = tv->next;
    }
    /* Fill the read parameters */
    rp = 0;
    for (j=0;j<ds->n;j++){
      if (labels[j]!=mode) continue;
      if (startPos[j]<0) continue;

      if (tpreadsStart[j]<0 || (tpreadsStart[j]+prWidth>(ds->features)[j])){
	printf("startpos: %d preadsstart: %d maxallowed: %d\n",startPos[j],tpreadsStart[j],maxallowed);
	continue;
      }
      rpt = getBackScore(prWidth,tpreadsStart[j],posreadsback[j]);
      //if (fabs(rpt)>0) rpt = log(rpt);
      rp += rpt;

      tv = tr->readsWin;
      for(k=0;k<prWidth;k++){
	(tv->modeReadsCount)[(int)((ds->posreads)[j][tpreadsStart[j]+k])]++;
	tv = tv->next;
      }
    }
    /* Find the log likelihood score of these params */
    tv = tr->readsWin;
    s = 0;
    for (j=0;j<prWidth;j++){
      for(k=0;k<2;k++){
	s += (tv->modeReadsCount)[k]*log((tv->modeReadsCount)[k]+(m->pcReads)[k]);
      }
      s = s-(m->counts)[mode]*log((m->counts)[mode]+totpc);
      tv = tv->next;
    }
    s = s-rp; // log likelihood when reads start from pos i
    
    poslikes[maxallowed-i] = s;
    //printf("Posn: %d\tposlikes: %.12lf\n",maxallowed-i,s);
    if (s>maxlike) maxlike = s;
    valid[maxallowed-i] = 1;
  }
  sum = 0;
  for (i=0;i<maxallowed+maxrightwidth+1;i++){
    if(fabs(maxlike-poslikes[i])<700) valid[i]=1;
    else{
      valid[i]=0;
      poslikes[i]=0;
    }
    sum += fabs(poslikes[i]);
  }
  if(isinf(log(sum))){
    printf("Sum of the elements for sampling positive offset is 0\n");
    printf("Mode: %d Maxallowed: %d and motifwidth: %d\n",mode,maxallowed,(m->mWidth)[mode]);
    for (i=0;i<maxallowed+maxrightwidth+1;i++){
      printf("poslikes[%d]: %lf valid: %d\n",i,poslikes[i],valid[i]);
    }
    exit(0);
  }
  takeExp(poslikes,valid,maxallowed+maxrightwidth+1);
  //for (i=0;i<(maxallowed+(m->mWidth)[mode]-prWidth+1);i++)
  //  printf("%.12lf\n",poslikes[i]);

  ind = sample(poslikes,maxallowed+maxrightwidth+1,(double)(rand_r(seed))/(RAND_MAX)); // ind: reads should start at maxallowed-ind number of positions before motifstart
  posnBeforeMotif = maxallowed-ind;
  //printf("posnBeforeMotif: %d\tmaxallowed: %d\tind: %d\n",posnBeforeMotif,maxallowed,ind);
  for (i=0;i<ds->n;i++){
    if(labels[i]!=mode) continue;
    if(startPos[i]<0) continue;
    preadsStart[i] = startPos[i]-posnBeforeMotif;
    if (preadsStart[i] < 0){
      printf("Positive reads went beyond sequence for mode: %d\n",mode);
      printf("startpos: %d   posnBeforeMotif: %d   preadsStart: %d   maxallowed: %d   ind: %d\n",startPos[i],posnBeforeMotif,preadsStart[i],maxallowed,ind);
      exit(0);
    }
    if (rev){
      if ((startPos[i] >= ((ds->features)[i])/2) && (preadsStart[i] < ((ds->features)[i])/2)){
	printf("Positive reads came into forward strand for mode: %d\n",mode);
	printf("startpos: %d   posnBeforeMotif: %d   preadsStart: %d   maxallowed: %d   ind: %d\n",startPos[i],posnBeforeMotif,preadsStart[i],maxallowed,ind);
	exit(0);
      }
    }
  }
  //fillReadsPWM(m,ds,mode,startPos,labels,preadsStart,(m->posreadsparams)[mode].readsWin,prWidth,0);
  free(valid);
  free(poslikes);
  free(tpreadsStart);
  freeReadsWin(tr->readsWin);
  free(tr);
  return posnBeforeMotif;
}

int sampleNegativeOffset(model *m, dataSet *ds, int mode, int *labels, int *startPos, int *nreadsStart, int nrWidth, int nreadspart, int maxoffset, double **negreadsback, unsigned int *seed, int rev,int gobeyond){
  //negreadsStart can change at the end
  int maxrightwidth, maxallowed, maxleftwidth=0,left=0;
  int i,j,k,flag=0,ind,posnAfterMotif;
  readsContainer *tr;
  readParams *tv;
  int *tnreadsStart, *valid;
  double *neglikes,maxlike=0;
  double rn,rnt,s,totpc,sum;
  //printf("From neg offset \n");
  maxrightwidth = maxoffset;
  for (i=-(m->mWidth)[mode];i<maxrightwidth;i++){
    flag=0;
    for (j=0;j < ds->n;j++){
      if (labels[j]!=mode) continue;
      if (startPos[j]<0) continue;
      if(startPos[j]+(m->mWidth)[mode]+nrWidth + i > (ds->features)[j]){
	//printf("Neg offset, starpos: %d i:%d mwidth: %d\n",startPos[i],i,(m->mWidth)[mode]);
	flag=1;
	break;
      }
      if(rev){
	if((startPos[j] < (ds->features)[j]/2) && ((startPos[j]+(m->mWidth)[mode]+nrWidth+i) > (ds->features)[j]/2)){
	  //printf("Neg offset, starpos: %d i:%d mwidth: %d\n",startPos[i],i,(m->mWidth)[mode]);
	  flag=1 ;
	  break;
	}
      }
    }
    if (flag) break;
  }
  maxallowed = i;
  //printf("Mode: %d Maxallowed: %d and motifwidth: %d\n",mode,maxallowed,(m->mWidth)[mode]);
  if (gobeyond){
    flag = 0;
    for (i=1;i<nrWidth;i++){
      for (j=0;j<(ds->n);j++){
	if (labels[j]!=mode) continue;
	if (startPos[j]<0) continue;
	if (startPos[j]-i<0){
	  flag = 1;
	  break;
	}
	if (rev && startPos[j] > (ds->features)[j]/2 && startPos[j]-i<(ds->features)[j]/2){
	  flag = 1;
	  break;
	}
      }
      if (flag) break;
    }
    left = i-1;
  }
  maxleftwidth = left + (m->mWidth)[mode];
  //printf("left obtained: %d\n",left);

  tnreadsStart = (int *)malloc(sizeof(int)*(ds->n));
  neglikes = (double*)malloc(sizeof(double)*(maxleftwidth+maxallowed));
  valid = (int*)malloc(sizeof(int)*(maxleftwidth+maxallowed));
  totpc = (m->pcReads)[0]+(m->pcReads)[1];
  tr = (readsContainer*)malloc(sizeof(readsContainer));
  tr->readsWin = initializeReadParamsPerMode(nrWidth);
  
  maxlike = -1.0*INT_MAX;

  for (i=-left;i<((m->mWidth)[mode]+maxallowed);i++){
    /*Adjust read positions as we slide*/
    for(j=0;j<ds->n;j++){
      if(labels[j]!=mode) continue;
      if(startPos[j]<0) continue;
      tnreadsStart[j] = startPos[j] + i;
    }
    /*Initialize the neg reads matrix*/
    tv = tr->readsWin;
    for (j=0;j<nrWidth;j++){
      for(k=0;k<2;k++)
	(tv->modeReadsCount)[k]=0;
      tv = tv->next;
    }
    /*Fill up the read parameters*/
    rn = 0;
    for(j=0;j<ds->n;j++){
      if(labels[j]!=mode) continue;
      if (labels[j]<0) continue;
      if (startPos[j]<0) continue;
      if (tnreadsStart[j]+nrWidth>(ds->features)[j]){
	printf("startpos: %d nreadsstart: %d maxallowed: %d\n",startPos[j],tnreadsStart[j],maxallowed);
	exit(0);
      }
      if (tnreadsStart[j]<0){
	printf("startpos: %d nreadsstart: %d maxallowed: %d\n",startPos[j],tnreadsStart[j],maxallowed);
	continue;
      }
      
      rnt = getBackScore(nrWidth,tnreadsStart[j],negreadsback[j]);
      //if (fabs(rnt)>0) rnt = log(rnt);
      rn += rnt;

      tv = tr->readsWin;
      for(k=0;k<nrWidth;k++){
	(tv->modeReadsCount)[(int)((ds->negreads)[j][tnreadsStart[j]+k])]++;
	tv = tv->next;
      }
    }
    /*Find the log likelihood score of these params*/
    tv = tr->readsWin;
    s=0;
    for(j=0;j<nrWidth;j++){
      for (k=0;k<2;k++){
	s += (tv->modeReadsCount)[k]*log((tv->modeReadsCount)[k]+(m->pcReads)[k]);
      }
      s = s - (m->counts)[mode]*log((m->counts)[mode]+totpc);
      tv = tv->next;
    }
    s = s-rn;
    neglikes[i+left] = s;
    //printf("Posn: %d neglikes: %lf\n",i,neglikes[i+left]);
    if (neglikes[i+left]>maxlike) maxlike = neglikes[i+left];
    valid[i+left] = 1;
  }
  sum=0;
  for(i=0;i<(maxleftwidth+maxallowed);i++){
    if (fabs(maxlike-neglikes[i])<700) valid[i] = 1;
    else{
      valid[i]=0;
      neglikes[i]=0;
    }
    sum+=fabs(neglikes[i]);
    //printf("%lf %d\n",neglikes[i],valid[i]);
  }
  if (isinf(log(sum))) {
    printf("Sum of the elements for sampling negative offset is 0\n");
    printf("maxallowed: %d and motif width: %d\n",maxallowed,(m->mWidth)[mode]);
    exit(0);
  }
  takeExp(neglikes,valid,maxallowed+maxleftwidth);
  //for (i=0;i<maxallowed;i++)
  //printf("%lf %d\n",neglikes[i],valid[i]);

  ind = sample(neglikes,maxleftwidth+maxallowed,(double)(rand_r(seed))/(RAND_MAX));
  posnAfterMotif = ind-maxleftwidth;
    
  //printf("For negative offset. Maxallowed: %d, posnAfterMotif: %d\n",maxallowed,posnAfterMotif);

  for (i=0;i<ds->n;i++){
    if(labels[i]!=mode) continue;
    nreadsStart[i] = startPos[i] + (m->mWidth)[mode] + posnAfterMotif;
     if (nreadsStart[i]+nrWidth > (ds->features)[i]){
      printf("Negreads window going beyond sequence for mode:%d\n",mode);
      exit(0);
    }
     if(rev && (startPos[i] < (ds->features)[i]/2) && (nreadsStart[i]+nrWidth > (ds->features)[i]/2)){
      printf("Negreads going into reverse strand for mode:%d\n",mode);
      printf("startpos: %d\tmotifwidth: %d\tnegreadsstart: %d\n",startPos[i],(m->mWidth)[mode],nreadsStart[i]);
      exit(0);
    }
  }
  //fillReadsPWM(m,ds,mode,startPos,labels,nreadsStart,(m->negreadsparams)[mode].readsWin,nrWidth,1);
  free(valid);
  free(neglikes);
  free(tnreadsStart);
  freeReadsWin(tr->readsWin);
  free(tr);
  return posnAfterMotif;
}
