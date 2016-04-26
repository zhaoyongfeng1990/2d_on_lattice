#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "lattice2d.h"


int main(int argc, const char * argv[])
{
  time_t it, ft;
  it=time(0);

  status cstatus;   //the structure for status of current system
  setDefaultStatus(&cstatus);
  int* lattice=(int*) malloc(cstatus.NumSite*sizeof(int));  //0: empty; -1: blocked; n: occupied by n-th cell
  ecoli* ecoliList=(ecoli*) malloc(cstatus.NumOfCells*sizeof(ecoli));
  //printf("test\n");
  generateLattice(lattice, &cstatus);
  putParticles(lattice, ecoliList, &cstatus);
  checkIfBlocked(&cstatus, lattice, ecoliList);

  //snapshotLattices(lattice, ecoliList, &cstatus, 0);
  double totalTime=cstatus.totalTime;
  int fileIdx=0;
  FILE* outfile=fopen("result.txt","w");
  int step=0;
  while (cstatus.time<totalTime)
  {
    ++step;
    int ifHopping=iterate(lattice, ecoliList, &cstatus);
    if(1==ifHopping)
    {
      fprintf(outfile, "%f ", cstatus.time);
      double sum=0;
      for (int i = 0; i < cstatus.NumOfCells; ++i)
      {
        int dx=ecoliList[i].deltaX;
        int dy=ecoliList[i].deltaY;
        sum+=dx*dx+dy*dy;
        //fprintf(outfile, "%d ", dx*dx+dy*dy);
      }
      sum=sum/cstatus.NumOfCells;
      fprintf(outfile, "%f \n", sum);
    }
    //for (int i = 0; i < cstatus.NumOfCells; ++i)
    //printf("%f, %d, %d \n", cstatus.time, ecoliList[i].deltaX, ecoliList[i].deltaY);
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
