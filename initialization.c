#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "lattice2d.h"
#include "mt64.h"

//generate a lattices with obstruction. 0 denotes empty site, and 1 denotes blocked site.
void generateLattice(int* lattice, const status* pstatus)
{
  double GelConcentration=pstatus->GelConcentration;
  int NumSite=pstatus->NumSite;

  for (int i = 0; i < NumSite; ++i)
  {
    if (genrand64_real1()<GelConcentration)   //each site has probability GelConcentration to be blocked.
    {
      lattice[i]=-1;
    }
    else
    {
      lattice[i]=0;
    }
  }
}

//uniformly sample a new direction of bacteria.
void sampleDirection(ecoli* cell)
{
  int randomDirection=(int)(genrand64_int64()%8); //a random number from 0 to 7
  switch (randomDirection)
  {
    case 0:
      cell->directionX=1;
      cell->directionY=0;
      cell->ifDiagnal=0;
      break;
    case 1:
      cell->directionX=1;
      cell->directionY=1;
      cell->ifDiagnal=1;
      break;
    case 2:
      cell->directionX=0;
      cell->directionY=1;
      cell->ifDiagnal=0;
      break;
    case 3:
      cell->directionX=-1;
      cell->directionY=1;
      cell->ifDiagnal=1;
      break;
    case 4:
      cell->directionX=-1;
      cell->directionY=0;
      cell->ifDiagnal=0;
      break;
    case 5:
      cell->directionX=-1;
      cell->directionY=-1;
      cell->ifDiagnal=1;
      break;
    case 6:
      cell->directionX=0;
      cell->directionY=-1;
      cell->ifDiagnal=0;
      break;
    case 7:
      cell->directionX=1;
      cell->directionY=-1;
      cell->ifDiagnal=1;
      break;
  }
}

//put particles on lattices. lattice[idx]==2 means the site is occupied by a cell.
void putParticles(int* lattice, ecoli* ecoliList, const status* pstatus)
{
  int NumOfCells=pstatus->NumOfCells;
  int LatticeDim=pstatus->LatticeDim;

  for (int i = 0; i < NumOfCells; ++i)
  {
    int randomX, randomY, idx;  //random coordinates and index in lattice array
    do
    {
      randomX=floor(LatticeDim*genrand64_real2());
      randomY=floor(LatticeDim*genrand64_real2());
      idx=randomY*LatticeDim+randomX;   //index in lattice array
    } while (!lattice[idx]);
    lattice[idx]=i;   //the site is occupied
    ecoliList[i].posX=randomX;
    ecoliList[i].posY=randomY;
    ecoliList[i].deltaX=0;
    ecoliList[i].deltaY=0;
    sampleDirection(&ecoliList[i]);
  }
}

//initialize pstatus with the pre-setting values
void setDefaultStatus(status* pstatus)
{
  FILE* parafile=fopen("parameters.txt", "r");

  int setLatticeDim;
  int setNumOfCells;
  double setGelConcentration;
  double setTumblingRate;
  double setHoppingRate;

  fscanf(parafile, "%d/n", &setLatticeDim);
  fscanf(parafile, "%d/n", &setNumOfCells);
  fscanf(parafile, "%lf/n", &setGelConcentration);
  fscanf(parafile, "%lf/n", &setTumblingRate);
  fscanf(parafile, "%lf/n", &setHoppingRate);

  fclose(parafile);
  int setNumSite=setLatticeDim*setLatticeDim;
  double setdiagHoppingRate=setHoppingRate/1.414213562373095;

  pstatus->time=0;
  pstatus->GelConcentration=setGelConcentration;
  pstatus->TumblingRate=setTumblingRate;
  pstatus->HoppingRate=setHoppingRate;
  pstatus->diagHoppingRate=setdiagHoppingRate;
  pstatus->LatticeDim=setLatticeDim;
  pstatus->NumOfCells=setNumOfCells;
  pstatus->NumSite=setNumSite;

  pstatus->movingEcoliEdge=(int*) malloc(setNumOfCells*sizeof(int));
  pstatus->movingEcoliDiag=(int*) malloc(setNumOfCells*sizeof(int));
  pstatus->NumOfMovingCellsEdge=0;
  pstatus->NumOfMovingCellsDiag=0;
}

void destructeStatus(status* pstatus)
{
  free(pstatus->movingEcoliEdge);
  free(pstatus->movingEcoliDiag);
}

void checkingIfBlocked(status* pstatus, const int* lattice, ecoli* ecoliList)
{
  int NumOfCells=pstatus->NumOfCells;
  int LatticeDim=pstatus->LatticeDim;
  int NumOfMovingCellsEdge=0;
  int NumOfMovingCellsDiag=0;
  int* movingEcoliDiag=pstatus->movingEcoliDiag;
  int* movingEcoliEdge=pstatus->movingEcoliEdge;

  for (int i = 0; i < NumOfCells; ++i)  //loop over all cells
  {
    int nextPosX=(ecoliList[i].posX+ecoliList[i].directionX+LatticeDim)%LatticeDim;
    int nextPosY=(ecoliList[i].posY+ecoliList[i].directionY+LatticeDim)%LatticeDim;
    int idx=nextPosY*LatticeDim+nextPosX;
    if (0==lattice[idx])   //cell can move
    { 
      ecoliList[i].inEdgeOrDiag=ecoliList[i].ifDiagnal;
      ecoliList[i].ifMoving=1;
      if (ecoliList[i].ifDiagnal)
      {
        movingEcoliDiag[NumOfMovingCellsDiag]=i;
        ecoliList[i].whereInMove=NumOfMovingCellsDiag;
        ++NumOfMovingCellsDiag;
      }
      else
      {
        movingEcoliEdge[NumOfMovingCellsEdge]=i;
        ecoliList[i].whereInMove=NumOfMovingCellsEdge;
        ++NumOfMovingCellsEdge;
      }
    }
    else  //cell cannot move
    {
      ecoliList[i].ifMoving=0;
    }
  }
  pstatus->NumOfMovingCellsEdge=NumOfMovingCellsEdge;
  pstatus->NumOfMovingCellsDiag=NumOfMovingCellsDiag;
}

void checkingIfBlockedSingleCell(status* pstatus, const int* lattice, ecoli* ecoliList, const int idx)
{
  int LatticeDim=pstatus->LatticeDim;
  int NumOfMovingCellsEdge=pstatus->NumOfMovingCellsEdge;
  int NumOfMovingCellsDiag=pstatus->NumOfMovingCellsDiag;
  int* movingEcoliDiag=pstatus->movingEcoliDiag;
  int* movingEcoliEdge=pstatus->movingEcoliEdge;

  int nextPosX=(ecoliList[idx].posX+ecoliList[idx].directionX+LatticeDim)%LatticeDim;
  int nextPosY=(ecoliList[idx].posY+ecoliList[idx].directionY+LatticeDim)%LatticeDim;
  int newIdx=nextPosY*LatticeDim+nextPosX;

  if (ecoliList[idx].ifMoving)  //cell can move before the event, it must be in movingEcoli list
  {
    if (0==lattice[newIdx])   //now cell can move
    {
      if (ecoliList[idx].ifDiagnal != ecoliList[idx].inEdgeOrDiag)  //need to swap between movingEcoliEdge and movingEcoliDiag
      {
        if (ecoliList[idx].inEdgeOrDiag)
        {
          //deleting from movingEcoliDiag
          movingEcoliDiag[ecoliList[idx].whereInMove]=movingEcoliDiag[NumOfMovingCellsDiag-1];
          --NumOfMovingCellsDiag;

          //adding to movingEcoliEdge
          movingEcoliEdge[NumOfMovingCellsEdge]=idx;
          ++NumOfMovingCellsEdge;

          ecoliList[idx].whereInMove=NumOfMovingCellsEdge;
          ecoliList[idx].inEdgeOrDiag=0;

          pstatus->NumOfMovingCellsEdge=NumOfMovingCellsEdge;
          pstatus->NumOfMovingCellsDiag=NumOfMovingCellsDiag;
        }
        else
        {
          //deleting from movingEcoliEdge
          movingEcoliEdge[ecoliList[idx].whereInMove]=movingEcoliEdge[NumOfMovingCellsDiag-1];
          --NumOfMovingCellsEdge;

          //adding to movingEcoliDiag
          movingEcoliDiag[NumOfMovingCellsDiag]=idx;
          ++NumOfMovingCellsDiag;

          ecoliList[idx].whereInMove=NumOfMovingCellsDiag;
          ecoliList[idx].inEdgeOrDiag=1;

          pstatus->NumOfMovingCellsEdge=NumOfMovingCellsEdge;
          pstatus->NumOfMovingCellsDiag=NumOfMovingCellsDiag;
        }
      }
    }
    else  //now cell cannot move, deleting from corresponding movingEcoli array
    {
      ecoliList[idx].ifMoving=0;
      if (ecoliList[idx].inEdgeOrDiag)
      {
        movingEcoliDiag[ecoliList[idx].whereInMove]=movingEcoliDiag[NumOfMovingCellsDiag-1];
        --NumOfMovingCellsDiag;

        pstatus->NumOfMovingCellsDiag=NumOfMovingCellsDiag;
      }
      else
      {
        movingEcoliEdge[ecoliList[idx].whereInMove]=movingEcoliEdge[NumOfMovingCellsDiag-1];
        --NumOfMovingCellsEdge;

        pstatus->NumOfMovingCellsEdge=NumOfMovingCellsEdge;
      }
    }
  }
  else  //cell cannot move before the event, it's not in any movingEcoli array
  {
    if (0==lattice[newIdx])  //now cell can move
    {
      ecoliList[idx].ifMoving=1;
      if (ecoliList[idx].ifDiagnal)
      {
        movingEcoliDiag[NumOfMovingCellsDiag]=idx;
        ++NumOfMovingCellsDiag;

        ecoliList[idx].whereInMove=NumOfMovingCellsDiag;
        ecoliList[idx].inEdgeOrDiag=1;

        pstatus->NumOfMovingCellsDiag=NumOfMovingCellsDiag;
      }
      else
      {
        movingEcoliEdge[NumOfMovingCellsEdge]=idx;
        ++NumOfMovingCellsEdge;

        ecoliList[idx].whereInMove=NumOfMovingCellsEdge;
        ecoliList[idx].inEdgeOrDiag=0;

        pstatus->NumOfMovingCellsEdge=NumOfMovingCellsEdge;
      }
    }
  }
}