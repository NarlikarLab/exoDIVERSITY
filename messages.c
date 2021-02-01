#include<stdlib.h>
#include<stdio.h>
void printMessage(int msgno){
  switch(msgno){
  case 0:
    printf("Error in memory allocation.\n");
    break;
  case 1:
    printf("Error in sampling.\n");
    break;
  default:
    printf("Error\n");
  }
}
