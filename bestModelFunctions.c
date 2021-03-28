#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<limits.h>
#include "motifAndReadsStructs.h"
#include "modelStructures.h"
#include "messages.h"
#include "motifAndReadsFunctions.h"
#include "traindata.h"


double getBestParams(model *tmpmodel, dataSet *ds, int index, double **background, double **posreadsback, double **negreadsback, int revFlag, int gobeyond, int *startPos, int *labels, int *posreadsStart, int *negreadsStart, double prevLikelihood){
  double currscore,rval,motval,maxscore,like;

  int j,L,v,i,maxlabel,maxindex,preldist,nreldist,oldlabel,oldstart,oldprStart,oldnrStart,flag=0;
  int rhs;
  int gbextraleft=0, gbextraright=0;
  maxscore = -1;
  maxlabel = -1;
  maxindex = -1;
  oldlabel = labels[index];
  oldstart = startPos[index];
  
  //printf("Index: %d\n",index);
  //Score each window and check if the score exceed the currscore for all modes
  //printf("Seq: %d Mode: %d motifstart: %d\tposreadstart: %d negreadstart: %d\n",index,labels[index],startPos[index],posreadsStart[index],negreadsStart[index]);
  L = (ds->features)[index];
  for(i=0;i<(tmpmodel->mode);i++){
    preldist = (tmpmodel->readMotifDist)[i].preadsMotif;
    nreldist = (tmpmodel->readMotifDist)[i].nreadsMotif;
    if (gobeyond){
      gbextraleft = nreldist-(tmpmodel->nreadsWidth)[i];
      gbextraright = ((tmpmodel->preadsWidth)[i]-preldist-(tmpmodel->mWidth)[i]-1);
    }
    if (gbextraright > 0 && gbextraright+(tmpmodel->mWidth)[i]>nreldist){
      rhs = (tmpmodel->mWidth)[i]+gbextraright;
    }
    else{
      if (nreldist < (tmpmodel->mWidth)[i]) 
	rhs = (tmpmodel->mWidth)[i];
      else
	rhs = nreldist;
    }
    // If posreldist > motifwidth do not check in this mode
    if (!gobeyond){
      if ((preldist<0) && (abs(preldist) + (tmpmodel->preadsWidth)[i])>(tmpmodel->mWidth)[i])
	continue;
    }
    else if (gobeyond){ //positive read window is beyond the right edge of motif
      if ((preldist<0) && (abs(preldist)>(tmpmodel->mWidth)[i]))
	continue;
      // if nreldist == nreadswidth :neg read window is at the left edge of motif 
    }
    if (preldist<=0){
      if (gbextraleft<0)
	j = abs(gbextraleft);
      else
	j=0;
    }
    else j = preldist;

    while(j<(L-rhs+1)){
      if (revFlag){
	if ((preldist>0) && (j>(L/2)-rhs) && (j<(L/2)+preldist)){
	  j=j+1;
	  continue;
	}
	else if ((preldist<=0) && (j>(L/2)-rhs) && (j<L/2)){
	  j=j+1;
	  continue;
	}
      }
      if ((ds->lookahead)[index][j] < (tmpmodel->mWidth)[i]){
	v = (ds->lookahead)[index][j];
	j= j+v+1;
	continue;
      }
      rval = scoreReads(tmpmodel,ds,i,j,preldist,nreldist,index,posreadsback[index],negreadsback[index]);
      if (isinf(log(rval))){
	printf("Reads window score <=0 i hill climbing. score: %lf\n",rval);
      exit(0);
      }
      motval = scoreMotif(tmpmodel,ds,i,j,index,background[index]);
      currscore = rval*motval;
      if (currscore > maxscore){
	maxscore = currscore;
	maxindex = j;
	maxlabel = i;
      }
      j=j+1;
    }
  }
  like = 0;
  oldprStart = posreadsStart[index];
  oldnrStart = negreadsStart[index];
  if(maxlabel !=labels[index]){ 
    if(maxlabel == -1 && labels[index]!= -1){
      //printf("label becomes -1\nSeq: %d\tstartPos: %d\tlabel: %d preldist: %d\tnreldist: %d\n",index,startPos[index],labels[index],(tmpmodel->readMotifDist)[labels[index]].preadsMotif,(tmpmodel->readMotifDist)[labels[index]].nreadsMotif);
      addRemoveDataPoint(tmpmodel,ds,labels,startPos,posreadsStart,negreadsStart,index,-1);
      labels[index] = -1; 
      posreadsStart[index] = -1;
      negreadsStart[index] = -1;
      flag = 1;
    }
    else if(maxlabel != -1 && labels[index] == -1){
      labels[index] = maxlabel; 
      startPos[index] = maxindex;
      posreadsStart[index] = maxindex - (tmpmodel->readMotifDist)[maxlabel].preadsMotif;
      negreadsStart[index] = maxindex + (tmpmodel->readMotifDist)[maxlabel].nreadsMotif - (tmpmodel->nreadsWidth)[maxlabel];
      addRemoveDataPoint(tmpmodel,ds,labels,startPos,posreadsStart,negreadsStart,index,1);
      flag = 2;
    }
    else{ // maxlabel and previous label differ and both are not -1   
      addRemoveDataPoint(tmpmodel,ds,labels,startPos,posreadsStart,negreadsStart,index,-1);
      labels[index] = maxlabel;
      startPos[index] = maxindex;
      posreadsStart[index] = maxindex - (tmpmodel->readMotifDist)[maxlabel].preadsMotif;
      negreadsStart[index] = maxindex + (tmpmodel->readMotifDist)[maxlabel].nreadsMotif - (tmpmodel->nreadsWidth)[maxlabel];
      addRemoveDataPoint(tmpmodel,ds,labels,startPos,posreadsStart,negreadsStart,index,1);
      flag = 3;
    }
    like = calculateLikelihood(tmpmodel,ds,labels,startPos,posreadsStart,negreadsStart,posreadsback,negreadsback,background,(tmpmodel->pcReads));
    if(like < prevLikelihood){
      switch(flag){
      case 1:
	labels[index] = oldlabel;
	posreadsStart[index] = oldprStart;
	negreadsStart[index] = oldnrStart;
	addRemoveDataPoint(tmpmodel,ds,labels,startPos,posreadsStart,negreadsStart,index,1);
	break;
      case 2:
	addRemoveDataPoint(tmpmodel,ds,labels,startPos,posreadsStart,negreadsStart,index,-1);
	labels[index] = -1;
	startPos[i] = -1;
	posreadsStart[index] = -1;
	negreadsStart[index] = -1;
	break;
      case 3:
	addRemoveDataPoint(tmpmodel,ds,labels,startPos,posreadsStart,negreadsStart,index,-1);
	labels[index] = oldlabel;
	startPos[index] = oldstart;
	posreadsStart[index] = oldprStart;
	negreadsStart[index] = oldnrStart;
	addRemoveDataPoint(tmpmodel,ds,labels,startPos,posreadsStart,negreadsStart,index,1);
	break;
      default:
	break;
      }
    }
    else prevLikelihood = like;
  }   
  else if (startPos[index] != maxindex){
    /* if (maxindex==-1){ */
    /*   printf("label becomes -1\nSeq: %d\tstartPos: %d\tlabel: %d preldist: %d\tnreldist: %d\n",index,startPos[index],labels[index],(tmpmodel->readMotifDist)[labels[index]].preadsMotif,(tmpmodel->readMotifDist)[labels[index]].nreadsMotif); */
    /* } */
    addRemoveDataPoint(tmpmodel,ds,labels,startPos,posreadsStart,negreadsStart,index,-1);
    startPos[index]=maxindex;
    posreadsStart[index] = startPos[index] - (tmpmodel->readMotifDist)[maxlabel].preadsMotif;
    negreadsStart[index] = startPos[index] + (tmpmodel->readMotifDist)[maxlabel].nreadsMotif - (tmpmodel->nreadsWidth)[maxlabel];
    addRemoveDataPoint(tmpmodel,ds,labels,startPos,posreadsStart,negreadsStart,index,1);
    like = calculateLikelihood(tmpmodel,ds,labels,startPos,posreadsStart,negreadsStart,posreadsback,negreadsback,background,(tmpmodel->pcReads));
    if (like < prevLikelihood){
      addRemoveDataPoint(tmpmodel,ds,labels,startPos,posreadsStart,negreadsStart,index,-1);
      startPos[index] = oldstart;
      //labels[index] = oldlabel;
      posreadsStart[index] = oldprStart;
      negreadsStart[index] = oldnrStart;
      addRemoveDataPoint(tmpmodel,ds,labels,startPos,posreadsStart,negreadsStart,index,1);
    }
    else prevLikelihood = like;
  }
  //printf("Seq: %d Mode: %d motstart: %d\tposreadstart: %d negreadstart: %d\n",index,labels[index],startPos[index],posreadsStart[index],negreadsStart[index]);

  return prevLikelihood;
}

double getBestMotifWidth(model *m, dataSet *ds, int mode, double **background, double **posreadsback, double **negreadsback, int revFlag, int gobeyond,int maxWidth, int minWidth, int *startPos, int *labels, int *posreadsStart, int *negreadsStart){
  // first best width on right and then on left
  int indL, indR;
  unsigned int dummyseed=0;
  double like,newlike;
  like = calculateLikelihoodMode(m,ds,labels,startPos,posreadsStart,negreadsStart,posreadsback,negreadsback,background,mode,(m->pcReads));
  //printf("For mode %d  ********** Like: %lf\n",mode,like);
  indR = sampleMotifWidthRight(ds,m,negreadsStart,labels,startPos,mode,maxWidth,minWidth,posreadsback,negreadsback,background,revFlag,gobeyond,&dummyseed,1);

  
  //printf("motwidth: %d\tposreldist: %d\tnegreldist: %d\n",(m->mWidth)[mode],(m->readMotifDist)[mode].preadsMotif,(m->readMotifDist)[mode].nreadsMotif);
  //printf("IndR: %d\n",indR);
  
  if (indR != 1){
    newlike = calculateLikelihoodMode(m,ds,labels,startPos,posreadsStart,negreadsStart,posreadsback,negreadsback,background,mode,(m->pcReads));
    //printf("Newlike on right: %lf \t indR: %d\n",newlike,indR);
    if (newlike < like){
      if (indR==0){ //Add back
	motifIncreaseRight(ds,m,background,labels,startPos,mode,negreadsStart,negreadsback);
      }
      else if (indR==2){ //Subtract
	motifDecreaseRight(ds,m,background,labels,startPos,mode,negreadsStart,negreadsback);
      }
    }
    else {
      like = newlike;
    }
  }
  //printf("motwidth: %d\tposreldist: %d\tnegreldist: %d\n",(m->mWidth)[mode],(m->readMotifDist)[mode].preadsMotif,(m->readMotifDist)[mode].nreadsMotif);
  
  //printf("After right width sampling Like: %lf\n",like);
  indL = sampleMotifWidthLeft(ds,m,posreadsStart,labels,startPos,mode,maxWidth,minWidth,posreadsback,negreadsback,background,revFlag,gobeyond,&dummyseed,1);
  
  //printf("motwidth: %d\tposreldist: %d\tnegreldist: %d\n",(m->mWidth)[mode],(m->readMotifDist)[mode].preadsMotif,(m->readMotifDist)[mode].nreadsMotif);
  //printf("IndL: %d\n",indL);
  
  if(indL != 1){
    newlike = calculateLikelihoodMode(m,ds,labels,startPos,posreadsStart,negreadsStart,posreadsback,negreadsback,background,mode,(m->pcReads));
    if(newlike < like){
      if(indL == 0){ //Add back
	motifIncreaseLeft(ds,m,background,labels,startPos,mode,posreadsStart,posreadsback);
      }
      else if (indL==2){ //subtract
	motifDecreaseLeft(ds,m,background,labels,startPos,mode,posreadsStart,posreadsback);
      }
    }
    else{
      like = newlike;  
    }
  }
  
  //printf("motwidth: %d\tposreldist: %d\tnegreldist: %d\n",(m->mWidth)[mode],(m->readMotifDist)[mode].preadsMotif,(m->readMotifDist)[mode].nreadsMotif);
  
  //printf("Before exiting mode: %d likelihood: %lf\n",mode,calculateLikelihoodMode(m,ds,labels,startPos,posreadsStart,negreadsStart,posreadsback,negreadsback,background,mode,(m->pcReads)));
  return like;
}

