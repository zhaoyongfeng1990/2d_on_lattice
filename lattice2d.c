#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "lattice2d.h"


int main(int argc, const char * argv[])
{
  time_t it, ft;
  it=time(0);

  status cstatus;   //the structure for status of current system
  int* lattice=(int*) malloc(cstatus.NumSite*sizeof(int));  //0: empty; -1: blocked; n: occupied by n-th cell
  ecoli* ecoliList=(ecoli*) malloc(cstatus.NumOfCells*sizeof(ecoli));
  setDefaultStatus(&cstatus);

  generateLattice(lattice, &cstatus);
  putParticles(lattice, ecoliList, &cstatus);
  checkingIfBlocked(&cstatus, lattice, ecoliList);



  //snapshotLattices(lattice, ecoliList, &cstatus, 0);
  double totalTime=100;
  int fileIdx=0;
  FILE* outfile=fopen("result.txt","w");
  while (cstatus.time<totalTime)
  {
    int ifHopping=iterate(lattice, ecoliList, &cstatus);
    if(1==ifHopping)
    {
      fprintf(outfile, "%f ", cstatus.time);
      for (size_t i = 0; i < cstatus.NumOfCells; ++i)
      {
        int dx=ecoliList[i].deltaX;
        int dy=ecoliList[i].deltaY;
        fprintf(outfile, "%d ", dx*dx+dy*dy);
      }
      fprintf(outfile, "\n");
    }
    //printf("%f, %d, %d \n", cstatus.time, ecoliList[0].deltaX, ecoliList[0].deltaY);
    //snapshotParticles(lattice, ecoliList, &cstatus, fileIdx);
    //++fileIdx;
  }
  fclose(outfile);
  /////////////////////

  /////////////////////
  destructeStatus(&cstatus);
  free(lattice);
  free(ecoliList);

  ft=time(0);
  printf("%ld\n",ft-it);
  return 0;
}
