#include <stdio.h>
#include <string.h>

#include "lattice2d.h"

//save a binary snapshot, containing all the information of the system at current time
void snapshot(const int* lattice, const ecoli* ecoliList, const status* pstatus, const int fileIndex)
{
  char filename[10];  //file name will be like s0000.bin
  char buffer[5];
  sprintf(buffer,"%04d", fileIndex);
  strcpy(filename, "s");
  strcat(filename, buffer);
  strcat(filename, ".bin");
  FILE* pOutputFile=fopen(filename,"wb");

  fwrite(pstatus, sizeof(status), 1, pOutputFile);

  fwrite(lattice, sizeof(int), pstatus->NumSite, pOutputFile);
  fwrite(ecoliList, sizeof(ecoli), pstatus->NumOfCells, pOutputFile);
  fclose(pOutputFile);
}

void snapshotLattices(const int* lattice, const ecoli* ecoliList, const status* pstatus, const int fileIndex)
{
  char filename[10];  //file name will be like s0000.bin
  char buffer[5];
  sprintf(buffer,"%04d", fileIndex);
  strcpy(filename, "l");
  strcat(filename, buffer);
  strcat(filename, ".bin");
  FILE* pOutputFile=fopen(filename,"wb");

  fwrite(pstatus, sizeof(status), 1, pOutputFile);

  fwrite(lattice, sizeof(int), pstatus->NumSite, pOutputFile);
  fclose(pOutputFile);
}

void snapshotParticles(const int* lattice, const ecoli* ecoliList, const status* pstatus, const int fileIndex)
{
  char filename[10];  //file name will be like s0000.bin
  char buffer[5];
  sprintf(buffer,"%04d", fileIndex);
  strcpy(filename, "p");
  strcat(filename, buffer);
  strcat(filename, ".bin");
  FILE* pOutputFile=fopen(filename,"wb");

  fwrite(pstatus, sizeof(status), 1, pOutputFile);

  fwrite(ecoliList, sizeof(ecoli), pstatus->NumOfCells, pOutputFile);
  fclose(pOutputFile);
}

//read the time and parameters of the system from a binary snapshot
void readSnapshotHead(status* pstatus, const int fileIndex)
{
  char filename[10];
  char buffer[5];
  sprintf(buffer,"%04d", fileIndex);
  strcpy(filename, "s");
  strcat(filename, buffer);
  strcat(filename, ".bin");
  FILE* pInputFile=fopen(filename,"rb");

  fread(pstatus, sizeof(status), 1, pInputFile);

  fclose(pInputFile);
}

//read the data of lattices and cells from a binary snapshot
void readSnapshotData(int* lattice, ecoli* ecoliList, const status* pstatus, const int fileIndex)
{
  char filename[10];
  char buffer[5];
  sprintf(buffer,"%04d", fileIndex);
  strcpy(filename, "s");
  strcat(filename, buffer);
  strcat(filename, ".bin");
  FILE* pInputFile=fopen(filename,"rb");

  fseek(pInputFile, sizeof(status), SEEK_SET);
  fread(lattice, sizeof(int), pstatus->NumSite, pInputFile);
  fread(ecoliList, sizeof(ecoli), pstatus->NumOfCells, pInputFile);

  fclose(pInputFile);
}
