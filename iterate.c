#include <math.h>
#include <stdio.h>
#include "lattice2d.h"
#include "mt64.h"

//simulate one gillespie step, return denote the type of event: 0:tumble, 1:hop
int iterate(int* lattice, ecoli* ecoliList, status* pstatus)
{
  double TumblingRate=pstatus->TumblingRate;
  double HoppingRate=pstatus->HoppingRate;
  double diagHoppingRate=pstatus->diagHoppingRate;
  int* movingEcoliDiag=pstatus->movingEcoliDiag;
  int* movingEcoliEdge=pstatus->movingEcoliEdge;
  int NumOfCells=pstatus->NumOfCells;
  int NumOfMovingCellsEdge=pstatus->NumOfMovingCellsEdge;
  int NumOfMovingCellsDiag=pstatus->NumOfMovingCellsDiag;

  double tumblingProb=NumOfCells * TumblingRate;
  double diagProb=NumOfMovingCellsDiag * diagHoppingRate + tumblingProb;
  double totalProb=NumOfMovingCellsEdge * HoppingRate + diagProb;
  double r1=genrand64_real3();
  double tau=-log(r1)/totalProb;  //time increment
  pstatus->time+=tau;

  //finding the event which is going to happen
  double threshold=genrand64_real3()*totalProb;

  if (threshold<tumblingProb)  //tumbling event
  {
    int idx=floor(threshold/TumblingRate);
    sampleDirection(&ecoliList[idx]);
    checkIfBlockedSingleCell(pstatus, lattice, ecoliList, idx);
    return 0;
  }
  else   //hopping event
  {
    int idx, i;
    if (threshold<diagProb) //hopping along diagonal
    {
      idx=floor((threshold-tumblingProb)/diagHoppingRate);
      i=movingEcoliDiag[idx];
    }
    else //hopping along edge
    {
      idx=floor((threshold-diagProb)/HoppingRate);
      i=movingEcoliEdge[idx];
    }

    int LatticeDim=pstatus->LatticeDim;
    int posX=ecoliList[i].posX;
    int posY=ecoliList[i].posY;
    int index=posY*LatticeDim+posX;

    //the mod operation is to make periodic boundary condition
    //plus with LatticeDim helps eliminate the possible of getting a negative value
    int newX=(posX+ecoliList[i].directionX+LatticeDim)%LatticeDim;
    int newY=(posY+ecoliList[i].directionY+LatticeDim)%LatticeDim;
    int newIndex=newY*LatticeDim+newX;

    //lattice[index]=0;     //now the previous site is empty
    //lattice[newIndex]=i;  //and the new site is occupied
    ecoliList[i].posX=newX;
    ecoliList[i].posY=newY;
    ecoliList[i].deltaX+=ecoliList[i].directionX;
    ecoliList[i].deltaY+=ecoliList[i].directionY;
    checkIfBlockedSingleCell(pstatus, lattice, ecoliList, i);
    return 1;
  }
  return -1;
}
