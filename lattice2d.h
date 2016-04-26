//the default settings
// #define setLatticeDim 2048
// #define setDelta 1  //unit: um, useless for now.
// #define setNumOfCells 1000
// #define setGelConcentration 0.3
// #define setTumblingRate 0
// #define setHoppingRate 10

// static const int setLatticeDim=512*4;
// static const double setDelta=1;  //unit: um, useless for now.
// static const int setNumOfCells=1000;
// static const double setGelConcentration=0.3;
// static const double setTumblingRate=10;
// static const double setHoppingRate=10;
//
// //two constant for speed
// static const int setNumSite=setLatticeDim*setLatticeDim;
// static const double setdiagHoppingRate=setHoppingRate/1.414213562373095;

typedef struct
{
  int posX;     //position
  int posY;
  int deltaX;   //displacement
  int deltaY;
  int directionX;   //direction
  int directionY;
  int ifDiagnal;    //bool value, true if direction is along diagonal
  int ifMoving;     //bool value, true if cell is moving
  int whereInMove;  //the index of this cell in movingEcoli array
  int inEdgeOrDiag; //0: in movingEcoliEdge, 1: in movingEcoliDiag
} ecoli;

typedef struct
{
  double time;
  double totalTime;
  double GelConcentration;  //volume concentration of the gel
  double TumblingRate;
  double HoppingRate;
  double diagHoppingRate;   //hopping rate/sqrt(2)
  int* movingEcoliEdge;     //list of index of cells moving along edge
  int* movingEcoliDiag;     //list of index of cells moving along diagonal
  int LatticeDim;           //number of sites in a row
  int NumOfCells;           //number of cells
  int NumSite;              //total number of sites
  int NumOfMovingCellsEdge; //number of cells moving along edge
  int NumOfMovingCellsDiag; //number of cells moving along diagonal
} status;  //status and parameters of the system

//generate a lattices with obstruction. 0 denotes empty site, and 1 denotes blocked site.
void generateLattice(int* lattice, const status* pstatus);

//put particles on lattices. lattice[idx]==2 means the site is occupied by a cell.
void putParticles(int* lattice, ecoli* ecoliList, const status* pstatus);

//generate movingEcoliEdge and movingEcoliDiag
void checkIfBlocked(status* pstatus, const int* lattice, ecoli* ecoliList);

//modifying movingEcoliEdge and movingEcoliDiag after one event
void checkIfBlockedSingleCell(status* pstatus, const int* lattice, ecoli* ecoliList, const int idx);

//initialize pstatus with the pre-setting values
void setDefaultStatus(status* pstatus);

//destruction function
void destructeStatus(status* pstatus);

//uniformly sample a new direction of bacteria.
void sampleDirection(ecoli* cell);

//save a binary snapshot, containing all the information of the system at current time
void snapshot(const int* lattice, const ecoli* ecoliList, const status* pstatus, const int fileIndex);


void snapshotLattices(const int* lattice, const ecoli* ecoliList, const status* pstatus, const int fileIndex);
void snapshotParticles(const int* lattice, const ecoli* ecoliList, const status* pstatus, const int fileIndex);

//read the time and parameters of the system from a binary snapshot
void readSnapshotHead(status* pstatus, const int fileIndex);

//read the data of lattices and cells from a binary snapshot
void readSnapshotData(int* lattice, ecoli* ecoliList, const status* pstatus, const int fileIndex);

//simulate one gillespie step, return denote the type of event: 0:tumble, 1:hop, 2:blocked
int iterate(int* lattice, ecoli* ecoliList, status* pstatus);
