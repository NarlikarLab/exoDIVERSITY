#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<math.h>
#include<string.h>
#include<limits.h>
#include<float.h>
#include "motifAndReadsStructs.h"
#include "modelStructures.h"
#include "motifAndReadsFunctions.h"
#include "messages.h"
#include "traindata.h"
#include "bestModelFunctions.h"

const float modesPrior = 1.0;
const int defaultReadWinsize = 20;
double getBackScore(int width, int start, double *backseq){
  double score=0;
  int i;
  if (start < 0) return score;
  for (i=start;i<(start+width);i++){
    if (backseq[i] < 0.000001){
      printf("Background is close to zero: %lf at position: %d\n",backseq[i],i);
      exit(0);
    }
    score = score+log(backseq[i]);
  }
  return score;
}

double calculateLikelihoodMode(model *m, dataSet *ds,int *labels ,int *startPos, int *posreadsStart, int *negreadsStart, double **posreadsback, double **negreadsback, double **background,int mode, float *pcReads){
  double s,s1,s2,totpseudocount;
  double rp,rn,rpt,rnt;
  motifStruct *mf;
  readParams *rf;
  int i,j;
  s=0;
  s = ((m->counts)[mode]+modesPrior)*(log((m->counts)[mode]+modesPrior) - log(m->n + (m->mode)*modesPrior));
  if (isnan(s)){
    printf("Nan is encountered in the modes component\n");
    printf("mode: %d \t seqs in mode: %d \ttotal seq: %d\n",mode,(m->counts)[mode],(m->n));
    exit(0);
  }
    
  
  s1=0;
  mf = (m->motifs)[mode].motif;
  for (i=0;i<(m->mWidth)[mode];i++){
    for(j=0;j<(m->featureValues);j++){
      s1 = s1 + ((mf->modeMotifCount)[j] + m->alpha)*(log((mf->modeMotifCount)[j]+m->alpha) - log((m->counts)[mode] + (m->featureValues)*(m->alpha))) ;
    }
    if (isnan(s1)){
      //printf("Nan got in mode: %d motif position: %d with seq count: %d\n",mode,i,(m->counts)[mode]);
	for (j=0;j<(m->featureValues);j++)
	  printf("%d ",(mf->modeMotifCount)[j]);
	printf("\n");
	exit(0);
    }
    mf = mf->next;
  }
  s=s+s1;
  s2=0;
  rp = 0.0;
  rn = 0.0;
  rpt = 0.0;
  rnt = 0.0;
  //to get the background of the motif part and reads parts
  for (i=0;i<(ds->n);i++){
    s1=0;
    if (labels[i]!=mode) continue;
    if (startPos[i] < 0) continue;
    else{
      if (startPos[i]>=0) {
	s1=getBackScore((m->mWidth)[mode],startPos[i],background[i]);
	rpt=getBackScore((m->preadsWidth)[mode],posreadsStart[i],posreadsback[i]); 
	rnt=getBackScore((m->nreadsWidth)[mode],negreadsStart[i],negreadsback[i]);
      }
      if (isnan(s1) || isnan(rpt) || isnan(rnt)){
	printf("Sequence no: %d, start position: %d posreadsstart: %d negreadsstart: %d\n",i,startPos[i],posreadsStart[i],negreadsStart[i]);
	exit(0);
      }
      rp+=rpt;
      rn+=rnt;
      s2+=s1;
    }
  }

  s=s-s2;
  s1=0;
  totpseudocount = pcReads[0]+pcReads[1];
  rf = (m->posreadsparams)[mode].readsWin;
  for (i=0;i<(m->preadsWidth)[mode];i++){
    for (j=0;j<2;j++){
      s1 = s1+((rf->modeReadsCount)[j]+pcReads[j])*(log((rf->modeReadsCount)[j]+pcReads[j]) - log((m->counts)[mode]+totpseudocount));     
    }

    if (isnan(s1)){
	printf("Nan got in mode: %d pos read window position: %d with seq count: %d\n",mode,i,(m->counts)[mode]);
	for (j=0;j<2;j++)
	  printf("%lf ",(rf->modeReadsCount)[j]);
	printf("\n");
	exit(0);
    }
    rf = rf->next;
  }
  s1 = s1 - rp;
  s = s+s1;
  s1 = 0.0;
  rf = (m->negreadsparams)[mode].readsWin;
  for(i=0;i<(m->nreadsWidth)[mode];i++){
    for (j=0;j<2;j++){
      s1 = s1+((rf->modeReadsCount)[j]+pcReads[j])*(log((rf->modeReadsCount)[j]+pcReads[j]) - log((m->counts)[mode]+totpseudocount));
    }

    if (isnan(s1)){
	printf("Nan got in mode: %d pos read window position: %d with seq count: %d\n",mode,i,(m->counts)[mode]);
	for (j=0;j<2;j++)
	  printf("%lf ",(rf->modeReadsCount)[j]);
	printf("\n");
	exit(0);
    }
    rf = rf->next;
  }
  s1 = s1 - rn;
  s = s+s1;  
  return s;
}

void addRemoveDataPoint(model *m, dataSet *ds, int *labels, int *startPos, int *posreadsStart, int *negreadsStart, int index, int ar){
  int i,mode;
  motifStruct *m1;
  readParams *r1,*r2;
  mode = labels[index];
  if (startPos[index]<0) return;
  (m->counts)[mode]=(m->counts)[mode]+ar;
  m->n=m->n +ar;
  m1=(m->motifs)[mode].motif;
  r1=(m->posreadsparams)[mode].readsWin;
  r2=(m->negreadsparams)[mode].readsWin;
  for(i=0;i<(m->mWidth)[mode];i++){
    if(i+startPos[index] >= (ds->features)[index]) {
      printf("Cannot add/remove sequence motif\n");
      exit(0);
    }
    (m1->modeMotifCount)[(ds->data)[index][startPos[index]+i]] = (m1->modeMotifCount)[(ds->data)[index][startPos[index]+i]] + ar;
    m1=m1->next;
  }
  for(i=0;i<(m->preadsWidth)[mode];i++){
    if(posreadsStart[index]+i >= (ds->features)[index]) {
      printf("Cannot add/remove pos read window\n");
      exit(0);
    }
    (r1->modeReadsCount)[(int)(ds->posreads)[index][posreadsStart[index]+i]] = (r1->modeReadsCount)[(int)(ds->posreads)[index][posreadsStart[index]+i]] + ar;
    r1 = r1->next;
  }
  for(i=0;i<(m->nreadsWidth)[mode];i++){
    if(negreadsStart[index]+i >= (ds->features)[index]){
      printf("Cannot add/remove neg reads window\n");
      exit(0);
    }
    (r2->modeReadsCount)[(int)(ds->negreads)[index][negreadsStart[index]+i]] = (r2->modeReadsCount)[(int)(ds->negreads)[index][negreadsStart[index]+i]] + ar;
    r2 = r2->next;
  }
}

double calculateLikelihood(model *m, dataSet *ds,int *labels, int *startPos, int *posreadsStart, int *negreadsStart, double **posreadsback,double **negreadsback, double **background, float *pcReads){
  double s,s1,s2,totpc;
  double rpt,rnt,rp,rn;
  int i, j,k;
  motifStruct *m1;
  readParams *r1,*r2;
  s=0;
  //Calculating the Posterior

  for (i=0; i < m->mode; i++){
    //if ((m->counts)[i]>0)
    s = s+((m->counts)[i]+modesPrior)*(log((m->counts)[i]+modesPrior)- log(m->n + (m->mode)*modesPrior));
  }
  if (isnan(s)){
    printf("Nan is encountered in the modes component\n");
    for (i=0; i < m->mode; i++){
      printf("seqs in mode: %d \ttotal seq: %d\n",(m->counts)[i],(m->n));
    }
    exit(0);
  }

  //printf("From modes: %lf \n",s);
  s1=0;
  for (i=0;i < (m->mode);i++){    
    m1 = (m->motifs)[i].motif;
    for (j=0;j< (m->mWidth)[i];j++){
      for (k=0;k< (m->featureValues); k++){
	s1 = s1 + ((m1->modeMotifCount)[k] + (m->alpha))*(log((m1->modeMotifCount)[k]+m->alpha) - log((m->counts)[i] + (m->featureValues)*(m->alpha)));
      }
      if (isnan(s1)){
	printf("Nan got in mode: %d motif position: %d with seq count: %d\n",i,j,(m->counts)[i]);
	for (k=0;k<(m->featureValues);k++)
	  printf("%d ",(m1->modeMotifCount)[k]);
	printf("\n");
	exit(0);
      }
      m1=m1->next;
    }

  }
  s =s+s1;
  s2=0;
  rp = 0.0;
  rn = 0.0;
  rpt = 0.0;
  rnt = 0.0;
  s1 = 0;
  // Background of the motif and the reads
  for (i=0;i<(ds->n);i++){
    if (labels[i]<0) continue;
    if (startPos[i]>=0) {
      s1=getBackScore((m->mWidth)[labels[i]],startPos[i],background[i]);
      
      rpt = getBackScore((m->preadsWidth)[labels[i]],posreadsStart[i],posreadsback[i]);
      rnt = getBackScore((m->nreadsWidth)[labels[i]],negreadsStart[i],negreadsback[i]);
      if (isnan(s1) || isnan(rpt) || isnan(rnt)){
	printf("Sequence no: %d, start position: %d posreadsstart: %d negreadsstart: %d\n",i,startPos[i],posreadsStart[i],negreadsStart[i]);
	exit(0);
      }
      rp+=rpt;
      rn+=rnt;
      s2 +=s1;
    }
  }
  
  s = s-s2;
  totpc = pcReads[0]+pcReads[1];
  for(i=0;i<(m->mode);i++){
    s1 = 0;
    //printf("Mode %d\n",i);
    r1 = (m->posreadsparams)[i].readsWin;
    for (j=0;j<(m->preadsWidth)[i];j++){
      for (k=0;k<2;k++){
	s1 = s1+((r1->modeReadsCount)[k] + pcReads[k])*(log((r1->modeReadsCount)[k]+pcReads[k]) - log((m->counts)[i]+totpc));
      }
      if (isnan(s1)){
	printf("Nan got in mode: %d pos read window position: %d with seq count: %d\n",i,j,(m->counts)[i]);
	for (k=0;k<2;k++)
	  printf("%lf ",(r1->modeReadsCount)[k]);
	printf("\n");
	exit(0);
      }
      r1=r1->next;
    }
    //s1 = s1 - (m->counts)[i]*(m->preadsWidth)[i]*log((m->counts)[i]+totpc);
    //printf("From positive reads: %lf\n",s1);
    s=s+s1;
    s1=0;
    r2 = (m->negreadsparams)[i].readsWin;
    for(j=0;j<(m->nreadsWidth)[i];j++){
      for(k=0;k<2;k++){
	s1 = s1+((r2->modeReadsCount)[k] + pcReads[k])*(log((r2->modeReadsCount)[k]+pcReads[k]) - log((m->counts)[i]+totpc));
      }
      if (isnan(s1)){
	printf("Nan got in mode: %d neg read window position: %d with seq count: %d\n",i,j,(m->counts)[i]);
	for (k=0;k<2;k++)
	  printf("%lf ",(r2->modeReadsCount)[k]);
	printf("\n");
	exit(0);
      }
      r2 = r2->next;
    }
    //s1 = s1 - (m->counts)[i]*(m->nreadsWidth)[i]*log((m->counts)[i]+totpc);
    //printf("From negative reads: %lf\n",s1);
    s=s+s1;
  }
  s=s-rp-rn;

  return s;
}

double scoreMotif(model *m,dataSet *ds,int mode,int start,int index,double *background){
  int i;
  double s;
  motifStruct *m1;
  s=1;
  m1=(m->motifs)[mode].motif;
  for(i=0;i<(m->mWidth)[mode];i++){
    if ((ds->data)[index][start+i] ==4) {
      printf("Motif contains N! seqno: %d start: %d\n",index,start);
      exit(0);
    }
    s = s * ((m1->modeMotifCount)[(ds->data)[index][start+i]] +0.5)/((m->counts)[mode] + (m->featureValues)*0.5);
    if (background[start+i] < 0.00001) continue;
    s = s / background[start+i];
    m1 = m1->next;
  }
  return s;
}

double scoreReads(model *m, dataSet *ds,int mode, int start, int posRelDist, int negRelDist, int index, double *posreadsback, double *negreadsback){
  int i,j,rc;
  readParams *r1,*r2;
  double s1,s2,s,tot;
  tot = (m->counts)[mode]+(m->pcReads)[0]+(m->pcReads)[1];
  s1=1;
  j = start-posRelDist;
  i=0;
  r1 = (m->posreadsparams)[mode].readsWin;
  while (i<(m->preadsWidth)[mode]){
    rc = (int)(ds->posreads)[index][j];
    s1=s1*(((r1->modeReadsCount)[rc]+(m->pcReads)[rc])/tot)/posreadsback[j];
    j=j+1;
    i=i+1;
    r1 = r1->next;
  }

  s2 = 1;
  r2 = (m->negreadsparams)[mode].readsWin;
  j = start+negRelDist-(m->nreadsWidth)[mode];
  i=0;
  while(i<(m->nreadsWidth)[mode]){
    rc = (int)(ds->negreads)[index][j];
    s2=s2*(((r2->modeReadsCount)[rc]+(m->pcReads)[rc])/tot)/negreadsback[j];
    i=i+1;
    j=j+1;
    r2 = r2->next;
  }
  s=s1*s2;
  return s;
}

void printModel(model *m,int mode){
  motifStruct *mc;
  readParams *r1, *r2;
  int i,j;
  mc = (m->motifs)[mode].motif;
  r1 = (m->posreadsparams)[mode].readsWin;
  r2 = (m->negreadsparams)[mode].readsWin;
  for (i=0;i<(m->mWidth)[mode];i++){
    for (j=0;j<(m->featureValues);j++){
      printf("%d  ",(mc->modeMotifCount)[j]);
    }
    mc = mc->next;
    printf("\n");
  }
  printf("Pos reads\n");
  for (i=0;i<(m->preadsWidth)[mode];i++){
    for(j=0;j<2;j++){
      printf("%lf  ",(r1->modeReadsCount)[j]);
    }
    r1=r1->next;
    printf("\n");
  }
  printf("Neg reads\n\n");
  for (i=0;i<(m->nreadsWidth)[mode];i++){
    for(j=0;j<2;j++){
      printf("%lf  ",(r2->modeReadsCount)[j]);
    }
    r2=r2->next;
    printf("\n");
  }
  printf("Motif Width: %d \t Pos reldist: %d \t Neg reldist: %d\n",(m->mWidth)[mode],(m->readMotifDist)[mode].preadsMotif,(m->readMotifDist)[mode].nreadsMotif);
}

int sampleStartPosn(model *m, dataSet *ds,int mode, int index, unsigned int *seed, int revFlag, int gobeyond, int posRelDist, int negRelDist, double *backgroundi, double *posreadsback, double *negreadsback){
  double *values,sum;
  int j,k,L,v,Ls,rhs,lhs=0;
  double motval,rval;
  int gbextraleft=0,gbextraright=0;
  if (gobeyond){
    gbextraleft = negRelDist-(m->nreadsWidth)[mode];
    gbextraright = ((m->preadsWidth)[mode]-posRelDist-(m->mWidth)[mode]);
  }
  if (gbextraright>0 && gbextraright+(m->mWidth)[mode]>negRelDist){
    rhs = (m->mWidth)[mode]+gbextraright;
  }
  else {
    if(negRelDist < (m->mWidth)[mode])
      rhs = (m->mWidth)[mode];
    else
      rhs = negRelDist;
  }
  if(posRelDist <= 0) {
    if (gbextraleft<0){
      values = (double *)malloc(sizeof(double)*((ds->features)[index] - rhs + gbextraleft+1));
      Ls=(ds->features)[index] - rhs + gbextraleft +1;
      j=abs(gbextraleft);
      lhs=j;
    }
    else{
      values = (double *)malloc(sizeof(double)*((ds->features)[index] - rhs+1));
      Ls=(ds->features)[index] - rhs +1;
      j=0;
      lhs=j;
    }
    if (!values) printMessage(0);
  }
  else{
    values = (double *)malloc(sizeof(double)*((ds->features)[index] - (posRelDist+rhs)+1));
    if (!values) printMessage(0);
    Ls= (ds->features)[index] - (posRelDist+rhs)+1;
    j=posRelDist;
    lhs=j;
  }
  
  //printf("Seq:%d\tposreldist:%d\tnegreldist:%d\tmotifwidth: %d\tmode: %d\n",index,posRelDist,negRelDist,(m->mWidth)[mode],mode);
  for (k=0;k<Ls;k++) values[k] = 0.0;
  sum=0;
  L=(ds->features)[index];
  k=0;
  while (j<L-rhs+1){
    if (revFlag){
      if ((posRelDist>0) && (j > (L/2)-rhs) && (j<(L/2)+posRelDist)){
	values[k]=0;
	j=j+1;
	k=k+1;
	continue;
      }
      else if ((posRelDist<=0) && (j > (L/2)-rhs) && (j<(L/2))){
	values[k]=0;
	j=j+1;
	k=k+1;
	continue;
      }
    }
    if ((ds->lookahead)[index][j]< (m->mWidth)[mode]){
      values[k]=0;
      v=(ds->lookahead)[index][j];
      j=j+v+1;
      k=k+v+1;
      continue;
    }
    rval = scoreReads(m,ds,mode,j,posRelDist,negRelDist,index,posreadsback,negreadsback);
    if (isinf(log(rval))){
      printf("Reads window score <=0. score: %lf\n",rval);
      exit(0);
    }
    motval = scoreMotif(m,ds,mode,j,index,backgroundi);
    values[k]=motval*rval;
    //printf("%d\t%lf \t %lf \t %lf \n",j,values[k],motval,rval);
    sum=sum+values[k];
    j=j+1;
    k=k+1;
  }
  if(isinf(log(sum))){
    printf("Sum is 0 for sampling motif start position for sequence: %d\n",index);
    exit(0);
  }
  v=0;
  v = sample(values,Ls,((double)rand_r(seed))/(RAND_MAX)); 
  //if (posRelDist > 0)
  //v = v+posRelDist;
  v = v+lhs;
  if (v==-1){
    printMessage(1);
    exit(1);
  }
  free(values);
  return v;
}

int sampleLabel(model *m, dataSet *ds, int *startPos, int index, unsigned int *seed, int revFlag, int gobeyond, double *background, double *posreadsback, double *negreadsback, int power){
  int i,v,motStart,L;
  double *values,sum=0,rval,motval;
  relDist *rd;
  int j;
  int gbextraleft=0,gbextraright=0;

  values = (double*)malloc(sizeof(double)*(m->mode));
  if (!values) printMessage(0);
  rd = m->readMotifDist;
  motStart = startPos[index];
  L = (ds->features)[index];
  for(i=0;i< m->mode; i++){
    if (gobeyond){
      gbextraleft = rd[i].nreadsMotif-(m->nreadsWidth)[i];
      gbextraright = ((m->preadsWidth)[i]-rd[i].preadsMotif-(m->mWidth)[i]);
    }
    values[i]=((m->counts)[i]+m->alpha)/(m->n+ (m->alpha)*(m->mode));
    if (rd[i].preadsMotif >0){
      if(motStart - rd[i].preadsMotif < 0){
	values[i]=0;
	continue;
      }
    }
    else{ //pos rel dist <= 0
      if (!gobeyond && (abs(rd[i].preadsMotif)+(m->preadsWidth)[i])>(m->mWidth)[i]){	
	values[i]=0;
	continue;
      }
      if (gbextraleft<0 && motStart+gbextraleft<0){
	values[i]=0;
	continue;
      }
      if (gbextraright>0 && motStart+(m->mWidth)[i]+gbextraright>L){
	values[i]=0;
	continue;
      }
    }
    if((rd[i].nreadsMotif <= (m->mWidth)[i]) && (motStart+(m->mWidth)[i] > L)){
      values[i] = 0;
      continue;
    }
    if(gbextraleft<0 && (rd[i].nreadsMotif <= (m->mWidth)[i]) && (motStart+gbextraleft < 0)){
      values[i] = 0;
      continue;
    }
    if((rd[i].nreadsMotif > (m->mWidth)[i]) && (motStart+rd[i].nreadsMotif > L)) {
      values[i]=0;
      continue;
    }
    else if(revFlag){
      if ((motStart < L/2) && (motStart + (rd[i].nreadsMotif)>(L/2))){
	values[i]=0;
	continue;
      } 
      if (gbextraright>0 && motStart<L/2 && (motStart+gbextraright)>L/2){
	values[i]=0;
	continue;
      }
      if ((motStart > L/2) && (rd[i].preadsMotif > 0) && (motStart - (rd[i].preadsMotif) < (L/2))){
	values[i]=0;
	continue;
      }
      if (gbextraleft<0 && (motStart > L/2) && (motStart+gbextraleft<L/2)){
	values[i]=0;
	continue;
      }
      if((rd[i].nreadsMotif <= (m->mWidth)[i]) && (motStart < L/2) && (motStart + (m->mWidth)[i] > L/2)){
	values[i]=0;
	continue;
      }
    }
   if(motifWithN((ds->data)[index],motStart,(m->mWidth)[i])){
     values[i]=0;
     continue;
   }
   //score the window motStart-posRelDist to motStart+ negRelDist
   rval=scoreReads(m,ds,i,motStart,(rd[i].preadsMotif), (rd[i].nreadsMotif),index, posreadsback, negreadsback);
   if(isinf(log(rval))){
     printf("Score for read window <=0. score: %lf\n",rval);
     exit(0);     
   }
   motval = scoreMotif(m,ds,i,motStart,index,background);

   if(power == 0){
     values[i] = values[i] *rval * motval;
   }
   else{
     values[i] = values[i] * pow(rval,1.0/power) * motval;
   }
   //if (motval > 100000000000) printf("Motifval is too big\n");
   if (motval < 0){
     printf("Motifval is too small\n");
     for (j=0;j<(m->mWidth)[i];j++){
       printf("(%lf %d)\t",background[motStart+j],(ds->data)[index][motStart+j]);
     }
     printf("\n");
     printModel(m,i);
     exit(0);
   }
   sum +=values[i];
   //printf("values at i: %d is %lf and sum: %lf\n",i,values[i],sum);
  }
  if (isinf(log(sum))){
    for (i=0;i<m->mode;i++) values[i]=0.0;
    printf("Sampling of labels have sum 0, seq no: %d with motstart: %d\n",index,motStart);
    for(i=0;i< m->mode; i++){
      printf("mode: %d\tposreldist: %d\tnegreldist: %d\tmotifwidth: %d\n",i,rd[i].preadsMotif,rd[i].nreadsMotif,(m->mWidth)[i]);
      printf("Values[%d]: %lf, #Seqs: %d, totalseqs: %d\n",i,values[i],(m->counts)[i],m->n);
      
    }
    exit(1);
  }
  rval = ((double)rand_r(seed))/(RAND_MAX);

  v=sample(values,m->mode,rval);
  if (v==-1){
    printMessage(1);
    exit(1);
  }
  free(values);
  return v;
}

void updateReadsStartPos(model *m, int *posreadsStart, int *negreadsStart,int newlabel, int newmotStart,int index){
  relDist rl;
  rl = (m->readMotifDist)[newlabel];
  posreadsStart[index] = newmotStart - (rl.preadsMotif);
  negreadsStart[index] = newmotStart + (rl.nreadsMotif)-(m->nreadsWidth)[newlabel];
}

void takeExp(double *values, int *valid,int n){
  int i,v=0;
  double avg=0.0;
  for (i=0;i<n;i++){
    if (!valid[i]) continue;
    else{
      avg+=values[i];
      v+=1;
    }
  }
  avg = avg/v;
  for (i=0;i<n;i++){
    if (!valid[i]) continue;
    else {
      values[i]=pow(exp(1),values[i]-avg);
      if (isinf(values[i])) values[i]=(double)(FLT_MAX);
    }
  }
}

int getMax(double *values,int *valid, int n){
  int i, maxindex;
  double *v,maxscore;
  v = (double *)malloc(sizeof(double)*n);
  for(i=0;i<n;i++) v[i] = values[i]*valid[i];
  maxscore = v[0];
  maxindex = 0;
  for (i=1;i<n;i++){
    if(v[i]>maxscore){
      maxscore = v[i];
      maxindex = i;
    }
  }
  return maxindex;
}

int sampleMotifWidthRight(dataSet *ds, model *m, int *nreadsStart, int *labels, int *startPos, int mode, int maxWidth, int minWidth, double **posreadsback, double **negreadsback, double **background,int revFlag, int gobeyond, unsigned int *seed,int hcflag){
  double *values;
  int ind,i,*valid;
  values = (double *)malloc(sizeof(double)*3);
  if (!values) printMessage(0);
  valid = (int*)malloc(sizeof(int)*3);
  for(i=0;i<3;i++) valid[i]=1;
  if ((m->mWidth)[mode]-1 < minWidth) {
    values[0]=0.0;
    valid[0]=0;
  }
  if ((m->readMotifDist)[mode].preadsMotif<0){
    if (!gobeyond && abs((m->readMotifDist)[mode].preadsMotif)+(m->preadsWidth)[mode] == (m->mWidth)[mode]){
      values[0]=0.0;
      valid[0]=0;
    }
    if (gobeyond && abs((m->readMotifDist)[mode].preadsMotif)==(m->mWidth)[mode]){
      values[0]=0.0;
      valid[0] = 0;
    }
  }
  if ((m->preadsWidth)[mode] >= (m->mWidth)[mode]){
    for (i=0;i<ds->n;i++){
      if (gobeyond && startPos[i]+(m->mWidth)[mode]==(m->nreadsWidth)[mode]){
	valid[0]=0;
	break;
      }
      if(startPos[i]-(m->readMotifDist)[mode].preadsMotif==0){
	valid[0]=0;
	break;
      }
      if (revFlag && startPos[i]>(ds->features)[i]/2 && startPos[i]-(m->readMotifDist)[mode].preadsMotif==(ds->features)[i]/2){
	valid[0]=0;
	break;
      }
    }
  }
  values[0] = 0.0;

  values[1] = scoreNochangeRight(ds,m,background,labels,startPos,mode,nreadsStart,negreadsback);
  if (((m->mWidth)[mode]+1)>maxWidth){
    values[2]=0.0;
    valid[2]=0;
  }
  else if(!gobeyond && (m->readMotifDist)[mode].nreadsMotif <= (m->mWidth)[mode]+(m->nreadsWidth)[mode]){
    values[2]=0.0;
    valid[2]=0;
  }// motif cannot increase inside nrwindow 
  else values[2] = scoreIncreaseRight(ds,m,background,labels,startPos,mode,nreadsStart,negreadsback,revFlag);
  if (fabs(values[2])<0.000000000000001){
    valid[2]=0;
  }  
  takeExp(values,valid,3);

  if (hcflag) ind = getMax(values,valid,3);
  else ind = sample(values,3,(double)(rand_r(seed))/(RAND_MAX));

  if (ind == 0) motifDecreaseRight(ds,m,background,labels,startPos,mode,nreadsStart,negreadsback);
  else if (ind == 2) motifIncreaseRight(ds,m,background,labels,startPos,mode,nreadsStart,negreadsback);
  //printf("Index chosen: %d\n",ind);
  free(values);
  free(valid);
  return ind;
}

int sampleMotifWidthLeft(dataSet *ds, model *m, int *preadsStart, int *labels, int *startPos, int mode, int maxWidth, int minWidth, double **posreadsback, double **negreadsback, double **background,int revFlag, int gobeyond, unsigned int *seed, int hcflag){
  double *values;
  int ind,i,*valid;
  values = (double *)malloc(sizeof(double)*3);
  if (!values) printMessage(0);
  valid = (int*)malloc(sizeof(int)*3);
  for (i=0;i<3;i++) valid[i]=1;
  if ((m->mWidth)[mode]-1 < minWidth) {
    values[0] = 0.0;
    valid[0] = 0;
  }
  if ((m->nreadsWidth)[mode] >= (m->mWidth)[mode]){
    for (i=0;i<ds->n;i++){
      if (startPos[i]+(m->readMotifDist)[mode].nreadsMotif+1 > (ds->features)[i]){
	valid[0]=0;
	break;
      }
      if (revFlag && startPos[i]<(ds->features)[i]/2 && (startPos[i]+(m->readMotifDist)[mode].nreadsMotif+1)>(ds->features)[i]/2){
	valid[0]=0;
	break;
      }
      if(gobeyond && startPos[i]+(m->readMotifDist)[mode].preadsMotif == (ds->features)[i]){
	valid[0]=0;
	break;
      }
    }
  }
  values[0] = 0.0;
  values[1] = scoreNochangeLeft(ds,m,background,labels,startPos,mode,preadsStart,posreadsback);
  if (((m->mWidth)[mode]+1)>maxWidth){
    values[2] = 0.0;
    valid[2] = 0;
  }
  /* else if((m->readMotifDist)[mode].preadsMotif <= (m->preadsWidth)[mode]){ */
  /*   values[2] = 0.0; */
  /*   valid[2] = 0; */
  /* } */
  else values[2] = scoreIncreaseLeft(ds,m,background,labels,startPos,mode,preadsStart,posreadsback,revFlag);
  if(fabs(values[2])<0.0000000000000001) valid[2]=0;

  takeExp(values,valid,3);

  if (hcflag) ind = getMax(values,valid,3);
  else ind = sample(values,3,(double)(rand_r(seed))/(RAND_MAX));

  if (ind == 0) motifDecreaseLeft(ds,m,background,labels,startPos,mode,preadsStart,posreadsback);
  else if (ind == 2) motifIncreaseLeft(ds,m,background,labels,startPos,mode,preadsStart,posreadsback);
  free(values);
  free(valid);
  //printf("Index: %d\n",ind);
  return ind;
}

double updateBestModel(model *m, dataSet *ds, int *labels, int *startPos, int *posreadsStart,int *negreadsStart, double **posreadsback, double **negreadsback, double **background, double maxLikelihood, int revFlag, int gobeyond, int maxWidth, int minWidth, FILE *fp, int iterNo){
  int i;
  //int j;
  double tmpLikelihood = 0.0;
  for (i=0;i<(ds->n);i++){
    tmpLikelihood = getBestParams(m,ds,i,background,posreadsback,negreadsback,revFlag,gobeyond,startPos,labels,posreadsStart,negreadsStart,maxLikelihood);
    if (maxLikelihood < tmpLikelihood){
      maxLikelihood = tmpLikelihood;
      if (fp!= NULL) fprintf(fp,"%lf\n",tmpLikelihood);
    }
    if (tmpLikelihood < maxLikelihood){
      printf("%d\t%lf\t%lf\n",iterNo,maxLikelihood,tmpLikelihood);
      exit(0);
    }
  }
  
  /* if(iterNo%10 == 0){ */
  /*   tmpLikelihood = 0; */
  /*   for(i=0;i<(m->mode);i++){ */
  /*     tmpLikelihood += getBestMotifWidth(m,ds,i,background,posreadsback,negreadsback,revFlag,gobeyond,maxWidth,minWidth,startPos,labels,posreadsStart,negreadsStart); */
  /*   } */
  /*   if(maxLikelihood < tmpLikelihood){ */
  /*     maxLikelihood = tmpLikelihood; */
  /*     if (fp!=NULL) fprintf(fp,"%lf\n",tmpLikelihood); */
  /*   } */
  /*   //printf("Maxlikelihood: %lf\n",maxLikelihood); */
  /* } */
  return maxLikelihood;
}

void printModelParams(model *m, int revflag, char *bestmodelfile){
  motifStruct *mc;
  readParams *r1,*r2;
  int i,j,k;
  FILE *fp;
  fp=NULL;
  if (bestmodelfile[0]!='\0') fp = fopen(bestmodelfile,"w");
  if (fp!=NULL) {
    //printf("Total modes: %d\n",(m->mode));
    fprintf(fp,"Modes: %d\n",(m->mode));
    //printf("Total features: %d\n",(m->featureValues));
    fprintf(fp,"Features: %d\n",(m->featureValues));
    //printf("Sequences: %d\n",(m->n));
    fprintf(fp,"Sequences: %d\n",(m->n));
    for (k=0;k<(m->mode);k++){
      fprintf(fp,"Motif: %d\n",k);
      //printf("Mode %d \n",k);
      mc = (m->motifs)[k].motif;
      for(i=0;i<(m->mWidth)[k];i++){
	for(j=0;j<(m->featureValues -1);j++){
	  fprintf(fp,"%d\t",(mc->modeMotifCount)[j]);
	}
	fprintf(fp,"%d\n",(mc->modeMotifCount)[j]);
	mc = mc->next;
      } 
      // printf("Motif printed\n");
    } // motifs written
    for (k=0;k<(m->mode);k++){
      fprintf(fp,"Posreads: %d\n",k);
      r1 = (m->posreadsparams)[k].readsWin;
      for(i=0;i<(m->preadsWidth)[k];i++){
	fprintf(fp,"%lf\t%lf\n",(r1->modeReadsCount)[0],(r1->modeReadsCount)[1]);
	r1 = r1->next;
      }
      //printf("Posreads for mode %d printed\n",k);
    } //posreads written
    for (k=0;k<(m->mode);k++){
      fprintf(fp,"Negreads: %d\n",k);
      r2 = (m->negreadsparams)[k].readsWin;
      for(i=0;i<(m->nreadsWidth)[k];i++){
	fprintf(fp,"%lf\t%lf\n",(r2->modeReadsCount)[0],(r2->modeReadsCount)[1]);
	r2 = r2->next;
      }
      //printf("Negreads for mode %d printed\n",k);
    } //negreads written
    fprintf(fp,"%s","Mode sequence counts\n");
    for (k=0;k<(m->mode);k++){
      fprintf(fp,"%d\t%d\n",k,(m->counts)[k]);
    } //seq counts written
    fprintf(fp,"%s","Motif widths \n");
    for (k=0;k<(m->mode);k++){
      fprintf(fp,"%d\t%d\n",k,(m->mWidth)[k]);
    } //motif widths written
    fprintf(fp,"%s","Positive strand read widths \n");
    for (k=0;k<(m->mode);k++){
      fprintf(fp,"%d\t%d\n",k,(m->preadsWidth)[k]);
    } //positive reads widths written
    fprintf(fp,"%s","Negative strand read widths \n");
    for (k=0;k<(m->mode);k++){
      fprintf(fp,"%d\t%d\n",k,(m->nreadsWidth)[k]);
    } //negative reads widths written
    fprintf(fp,"%s","Relative distances\n");
    for (k=0;k<(m->mode);k++){
      fprintf(fp,"%d\t%d\t%d\n",k,(m->readMotifDist)[k].preadsMotif,(m->readMotifDist)[k].nreadsMotif);
    } //relative distances written
    fprintf(fp,"Alpha: %lf\n",(m->alpha));
    fprintf(fp,"%s","Pseudo count for reads\n");
    fprintf(fp,"%lf\t%lf\n",(m->pcReads)[0],(m->pcReads)[1]);
    fclose(fp);
  }
  //printf("All other details printed\n");

}


trainOut* trainData(dataSet *ds,int mode, float alpha, float *pcReads, unsigned int seed, double **background, double **posreadsback, double **negreadsback, int *mWidth, int rWidth,int minWidth, int maxWidth, int posOffset, int negOffset, int revFlag, int gobeyond, char *filename,char *bestmodelfile){
  trainOut *to=NULL;
  int i,j,k,iterations,oldLabel,oldStart,flag=0,maxoffset,power,countIter;
  int *labels, *startPos, *posreadsStart, *negreadsStart,*preadsStart,*nreadsStart;
  int *preadspart, *nreadspart, *bmpreadspart, *bmnreadspart, *prWidth, *nrWidth;
  int *lpc, *spc, *widths;
  double *modeWiseLikelihood;
  double maxLikelihood,tmpLikelihood;
  model *m;
  FILE *fp;
  int mw1,mw2,pr,nr;
  //int rw1,rw2;
  fp=NULL;
  if (filename[0]!='\0') fp = fopen(filename,"w");

  to=(trainOut*)malloc(sizeof(trainOut));
  if (!to) printMessage(0);
  labels = (int*)malloc(sizeof(int)*ds->n);
  if (!labels) printMessage(0);
  startPos = (int*)malloc(sizeof(int)*ds->n);
  if (!startPos) printMessage(0);
  lpc = (int*)malloc(sizeof(int)*(ds->n));
  if(!lpc) printMessage(0);
  spc = (int*)malloc(sizeof(int)*(ds->n));
  if(!spc) printMessage(0);
  widths = (int*)malloc(sizeof(int)*mode);
  if(!widths) printMessage(0);
  
  posreadsStart = (int*)malloc(sizeof(int)*(ds->n));
  if(!posreadsStart) printMessage(0);
  negreadsStart = (int*)malloc(sizeof(int)*(ds->n));
  if(!negreadsStart) printMessage(0);
  preadsStart = (int*)malloc(sizeof(int)*(ds->n));
  if (!preadsStart) printMessage(0);
  nreadsStart = (int*)malloc(sizeof(int)*(ds->n));
  if(!nreadsStart) printMessage(0);
  
  preadspart = (int*)malloc(sizeof(int)*mode);
  nreadspart = (int*)malloc(sizeof(int)*mode);
  bmpreadspart = (int*)malloc(sizeof(int)*mode);
  bmnreadspart = (int*)malloc(sizeof(int)*mode);
  prWidth = (int*)malloc(sizeof(int)*mode);
  nrWidth = (int*)malloc(sizeof(int)*mode);
  for (i=0;i<mode;i++){    
    prWidth[i]=rWidth;
    nrWidth[i]=rWidth;
  }
  for(i=0;i<mode;i++) {
    if (posOffset > 0)
      preadspart[i] = posOffset+prWidth[i];
    else if (!gobeyond && posOffset<=(mWidth[i]-prWidth[i]))
      preadspart[i]= posOffset;
    else if (gobeyond && posOffset<=mWidth[i])
      preadspart[i]= posOffset;
    else if (abs(posOffset) > mWidth[i]){
      posOffset= 0;
      preadspart[i] = posOffset+prWidth[i];
    }
    if (negOffset < 0){
      if (gobeyond && abs(negOffset)>=(mWidth[i]+nrWidth[i]))
	negOffset = 0;
      else if(!gobeyond && abs(negOffset)>mWidth[i])
	negOffset = 0;
    }
    nreadspart[i]= negOffset+nrWidth[i];
  }
  maxoffset = 10;
  countIter = 0;
  initializeLabelStartPos(ds,labels,startPos,mode,mWidth,prWidth,nrWidth,preadspart,nreadspart,revFlag,gobeyond,&seed);

  initializeReadWinStartPos(posreadsStart,negreadsStart,startPos,labels,ds->n,mWidth, prWidth,nrWidth,preadspart,nreadspart);

  m=createModel(mode,ds,labels,startPos,posreadsStart,negreadsStart,alpha,pcReads,mWidth,prWidth,nrWidth,preadspart,nreadspart);
  copyarr(labels,lpc,ds->n);
  copyarr(startPos,spc,ds->n);
  copyarr(m->mWidth,widths,m->mode);
  copyarr(m->preadsWidth,prWidth,m->mode);
  copyarr(m->nreadsWidth,nrWidth,m->mode);
  copyarr(posreadsStart,preadsStart,ds->n);
  copyarr(negreadsStart,nreadsStart,ds->n);
  copyarr(preadspart, bmpreadspart, m->mode);
  copyarr(nreadspart, bmnreadspart, m->mode);

  modeWiseLikelihood = (double*)malloc(sizeof(double)*mode);
  if (!modeWiseLikelihood) printMessage(0);
  for (i=0;i<mode;i++){  
    modeWiseLikelihood[i] = calculateLikelihoodMode(m,ds,lpc,spc,preadsStart,nreadsStart,posreadsback,negreadsback,background,i,pcReads);
    //printf("Posterior for mode:%d is %lf\n",i,modeWiseLikelihood[i]);
  }
  maxLikelihood = calculateLikelihood(m,ds,lpc,spc,preadsStart,nreadsStart,posreadsback,negreadsback,background,pcReads);
  tmpLikelihood = maxLikelihood;
  //printf("Likelihood: %lf\n",maxLikelihood);
  iterations = (ds->n);
  if (iterations<=2000) iterations *= 2;
  if (iterations >= 5000) iterations = 2000;
  power = iterations*10;
  //power = 0;
  i=0;
  j=0;
  
  while(j!=iterations){
    j++;
    //printf("Iteration: %d\n",j);
    for (i=0;i<(ds->n);i++){
    oldLabel = lpc[i];
    oldStart = spc[i];
    addRemoveDataPoint(m,ds,lpc,spc,preadsStart,nreadsStart,i,-1); 
    if (lpc[i]>=0) {
      spc[i] = sampleStartPosn(m,ds,lpc[i],i,&seed,revFlag,gobeyond,(m->readMotifDist)[lpc[i]].preadsMotif,(m->readMotifDist)[lpc[i]].nreadsMotif,background[i],posreadsback[i],negreadsback[i]);
      if (spc[i]<0){
	printf("sampling start pos makes it negative. start:%d label:%d\n",spc[i],lpc[i]);
	exit(0);
      }
    }
    else spc[i]=-1;    

    if(spc[i]>=0) {
      lpc[i] = sampleLabel(m, ds, spc, i, &seed, revFlag,gobeyond,background[i],posreadsback[i],negreadsback[i],power);
    }
    else lpc[i]=-1;
    if (spc[i]>=0) {
      updateReadsStartPos(m,preadsStart,nreadsStart,lpc[i],spc[i],i);
      addRemoveDataPoint(m,ds,lpc,spc,preadsStart,nreadsStart,i,1);
    }
    if (oldLabel == lpc[i] && oldStart == spc[i]){
      if (fp!=NULL) fprintf(fp,"%lf\n",tmpLikelihood);
      continue;
    }
    modeWiseLikelihood[oldLabel] = calculateLikelihoodMode(m,ds,lpc,spc,preadsStart,nreadsStart,posreadsback,negreadsback, background,oldLabel,pcReads);
    if(oldLabel!=lpc[i]) modeWiseLikelihood[lpc[i]] = calculateLikelihoodMode(m,ds,lpc,spc,preadsStart,nreadsStart,posreadsback,negreadsback, background,lpc[i],pcReads);
    tmpLikelihood = 0;
    for (k=0;k< m->mode;k++) { 
      tmpLikelihood += modeWiseLikelihood[k];
    }
    
    if (fp!=NULL) fprintf(fp,"%lf\n",tmpLikelihood);

    /* Update bestModel if current likelihood is greater then next */
    if(tmpLikelihood > maxLikelihood){
      copyarr(lpc, labels, ds->n);
      copyarr(spc, startPos, ds->n);
      copyarr(m->mWidth, widths, m->mode);
      copyarr(m->preadsWidth, prWidth, m->mode);
      copyarr(m->nreadsWidth, nrWidth, m->mode);      
      copyarr(preadsStart, posreadsStart, ds->n);
      copyarr(nreadsStart, negreadsStart, ds->n);
      copyarr(preadspart, bmpreadspart, m->mode);
      copyarr(nreadspart, bmnreadspart, m->mode);
      maxLikelihood = tmpLikelihood;
      countIter = 0;
    }
    else countIter+=1;
    }
    
    if (j%10==0){
      for (k=0; k < m->mode; k++){
    	if ((m->counts)[k] == 0) continue;
    	//printf("Width sampling. Mode: %d\n",k);
    	mw1 = sampleMotifWidthRight(ds,m,nreadsStart,lpc,spc,k,maxWidth,minWidth,posreadsback,negreadsback,background,revFlag,gobeyond,&seed,0);
    	//printf("Mode: %d right mw1: %d motif width: %d\n",k,mw1,(m->mWidth)[k]);
    	mw2 = sampleMotifWidthLeft(ds,m,preadsStart,lpc,spc,k,maxWidth,minWidth,posreadsback,negreadsback,background,revFlag,gobeyond,&seed,0);
    	//printf("Mode: %d left mw2: %d motif width: %d\n",k,mw2,(m->mWidth)[k]);
    	if (mw1 !=1 || mw2!= 1) countIter = 0;
    	if (mw1 == 0){
    	  nreadspart[k] += 1;
    	  }
    	else if (mw1 == 2){
    	  nreadspart[k] -= 1;
    	}
    	if (mw2 == 0){
    	  preadspart[k] += 1 ;
    	}
    	else if(mw2 == 2){
    	  preadspart[k] -=1 ;
    	}
      }
    }
    // Sample positive and negative offset
    for (k=0;k<m->mode;k++){
      if (prWidth[k] == 0 && nrWidth[k] == 0) continue;
      if ((m->counts)[k] == 0) continue;
      if(prWidth[k] > 0){
	//printf("Mode : %d ************ Sample positive offset (Current: %d)\n",k,preadspart[k]);
	pr = samplePositiveOffset(m,ds,k,lpc,spc,preadsStart,prWidth[k],preadspart[k],maxoffset,posreadsback,&seed,revFlag,gobeyond);
	preadspart[k] = pr;
	(m->readMotifDist)[k].preadsMotif = preadspart[k];
	//printf("New positive offset: %d for mode: %d\n",preadspart[k],k);
      }
      if (nrWidth[k] > 0){
	//printf("Sample negative offse: (Current: %d)t\n",nreadspart[k]);
	nr = sampleNegativeOffset(m,ds,k,lpc,spc,nreadsStart,nrWidth[k],nreadspart[k],maxoffset,negreadsback,&seed,revFlag,gobeyond);
	nreadspart[k] = nr + nrWidth[k];
	(m->readMotifDist)[k].nreadsMotif = nreadspart[k]+(m->mWidth)[k];
	//printf("New negative offset: %d for mode: %d\n",nreadspart[k],k);
      }
      //printf("Mode: %d poff : %d noff:%d\n",k,pr,nr);
      // Recompute all the read params
      getReadParams(m, ds, spc, lpc, preadsStart, nreadsStart);
    }
    
    tmpLikelihood = calculateLikelihood(m,ds,lpc,spc,preadsStart,nreadsStart,posreadsback,negreadsback,background,pcReads);
    power = power/2;
    if (fp!=NULL) fprintf(fp,"%lf\n",tmpLikelihood);
    if(tmpLikelihood > maxLikelihood){
      copyarr(lpc, labels, ds->n);
      copyarr(spc, startPos, ds->n);
      copyarr(m->mWidth, widths, m->mode);
      copyarr(m->preadsWidth, prWidth, m->mode);
      copyarr(m->nreadsWidth, nrWidth, m->mode);      
      copyarr(preadsStart, posreadsStart, ds->n);
      copyarr(nreadsStart, negreadsStart, ds->n);
      copyarr(preadspart, bmpreadspart, m->mode);
      copyarr(nreadspart, bmnreadspart, m->mode);
      maxLikelihood = tmpLikelihood;
      countIter = 0;
    }
    if (countIter >= 5*(ds->n)) break;
    // There was no update to maxlikelihood for 5n iterations
  }

  //printf("Sampler ran for %d iterations and countIter value: %d\n",j,countIter);
  //Hill Climbing
  //Build new model according to the parameters corresponding to the maxLikelihood
  freeModel(m);

  m=createModel(mode,ds,labels,startPos,posreadsStart,negreadsStart,alpha,pcReads,widths,prWidth,nrWidth,bmpreadspart,bmnreadspart);
  tmpLikelihood=0;
  for (k=0;k< m->mode;k++) {
    modeWiseLikelihood[k] = calculateLikelihoodMode(m,ds,labels,startPos,posreadsStart,negreadsStart,posreadsback,negreadsback,background,k,pcReads);
    tmpLikelihood += modeWiseLikelihood[k];
  }
  j = 0;
  flag = 1;
  while((flag == 1) && (j < 100)){
    tmpLikelihood = updateBestModel(m,ds,labels,startPos,posreadsStart,negreadsStart,posreadsback,negreadsback,background,maxLikelihood,revFlag,gobeyond,maxWidth,minWidth,fp,j);
    if (tmpLikelihood > maxLikelihood) {
      maxLikelihood = tmpLikelihood;
      flag = 1;
    }
    else flag = 0;
    j++;
  }
  maxLikelihood = calculateLikelihood(m,ds,labels,startPos,posreadsStart,negreadsStart,posreadsback,negreadsback,background,pcReads);
  
  copyarr(m->mWidth,widths,m->mode);
  copyarr(m->preadsWidth,prWidth,m->mode);
  copyarr(m->nreadsWidth,nrWidth,m->mode);

  printModelParams(m,revFlag,bestmodelfile);

  if (fp!=NULL) fclose(fp);
  assignTrainOut(to,ds->n,mode,labels,startPos,posreadsStart,negreadsStart,widths,prWidth,nrWidth,maxLikelihood); 

  freeModel(m);
  free(modeWiseLikelihood);
  free(lpc);
  free(spc);
  free(widths);
  free(prWidth);
  free(nrWidth);
  free(preadspart);
  free(nreadspart);
  free(preadsStart);
  free(nreadsStart);
  free(bmpreadspart);
  free(bmnreadspart);
  free(posreadsStart);
  free(negreadsStart);

  return to;
}

void freeTo(trainOut *to){
  free(to->labels); to->labels=NULL;
  free(to->startPos); to->startPos=NULL;
  free(to->preadsStart); to->preadsStart=NULL;
  free(to->nreadsStart); to->nreadsStart=NULL;
  free(to->motifWidth); to->motifWidth=NULL;
  free(to->preadsWidth); to->preadsWidth=NULL;
  free(to->nreadsWidth); to->nreadsWidth=NULL;
  free(to);
}
