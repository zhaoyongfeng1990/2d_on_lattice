#include <math.h>

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
      lattice[i]=1;
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
    } while (lattice[idx]>0);
    lattice[idx]=2;   //the site is occupied
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
  pstatus->time=0;
  pstatus->Delta=setDelta;
  pstatus->GelConcentration=setGelConcentration;
  pstatus->TumblingRate=setTumblingRate;
  pstatus->HoppingRate=setHoppingRate;
  pstatus->diagHoppingRate=setdiagHoppingRate;
  pstatus->LatticeDim=setLatticeDim;
  pstatus->NumOfCells=setNumOfCells;
  pstatus->NumSite=setNumSite;
}
