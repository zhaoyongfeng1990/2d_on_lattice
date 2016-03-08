#include <math.h>

#include "lattice2d.h"
#include "mt64.h"

//simulate one gillespie step, return denote the type of event: 0:tumble, 1:hop, 2:blocked
int iterate(int* lattice, ecoli* ecoliList, status* pstatus)
{
  double TumblingRate=pstatus->TumblingRate;
  double HoppingRate=pstatus->HoppingRate;
  double diagHoppingRate=pstatus->diagHoppingRate;
  int NumOfCells=pstatus->NumOfCells;

  double totalProb=NumOfCells*TumblingRate;
  for (int i = 0; i < NumOfCells; ++i)
  {
    if (ecoliList[i].ifDiagnal)
    {
      totalProb+=diagHoppingRate;
    }
    else
    {
      totalProb+=HoppingRate;
    }
  }

  double r1=genrand64_real3();
  double tau=-log(r1)/totalProb;  //time increment
  pstatus->time+=tau;

  //finding the event which is going to happen
  double threshold=genrand64_real3()*totalProb;
  double partialProb=0;
  for (int i = 0; i < NumOfCells; ++i)
  {
    partialProb+=TumblingRate;
    if (partialProb>threshold)   //the i-th cell goes to tumble
    {
      int dx=ecoliList[i].directionX;
      int dy=ecoliList[i].directionY;
      sampleDirection(&ecoliList[i]);
      return 0;
    }

    if (ecoliList[i].ifDiagnal)
    {
      partialProb+=diagHoppingRate;
    }
    else
    {
      partialProb+=HoppingRate;
    }

    if (partialProb>threshold) //the i-th cell runs to the next site
    {
      int LatticeDim=pstatus->LatticeDim;
      int posX=ecoliList[i].posX;
      int posY=ecoliList[i].posY;
      int index=posY*LatticeDim+posX;

      //the mod operation is to make periodic boundary condition
      //plus with LatticeDim helps eliminate the possible of getting a negative value
      int newX=(posX+ecoliList[i].directionX+LatticeDim)%LatticeDim;
      int newY=(posY+ecoliList[i].directionY+LatticeDim)%LatticeDim;
      int newIndex=newY*LatticeDim+newX;
      if (0==lattice[newIndex]) //check if the new site is empty
      {
        lattice[index]=0;     //now the previous site is empty
        lattice[newIndex]=2;  //and the new site is occupied
        ecoliList[i].posX=newX;
        ecoliList[i].posY=newY;
        ecoliList[i].deltaX+=ecoliList[i].directionX;
        ecoliList[i].deltaY+=ecoliList[i].directionY;
        return 1;
      }
      return 2;
    }
  }
  return -1;
}
