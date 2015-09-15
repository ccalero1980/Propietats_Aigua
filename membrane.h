#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>    
// #include <iostream>


// using namespace std;
#define REAL float
// #define FORMAT "%f"

class Water {

 public:
  REAL xW[20][20000], yW[20][20000], zW[20][20000], zVW[20][20000];
  REAL Lx, Ly, Lz;
  int NP, NW, NmolecsW, Nframes, numframe;     
  
  Water();
  REAL cpc(REAL, REAL);
  REAL InteractionLJ(int, int);
  REAL InteractionC(int, int);
  REAL area();
  void InterfacialLJ(int, int, REAL, int, REAL);
  void InterfacialC(int, REAL, int, REAL);
  void stats(REAL E1[1000], int ncasos, REAL* AvE,REAL* sigmaE);
  void readcoord(int, int);
  



};
class Heads {

 public:
  REAL xH[100][1000], yH[100][1000], zH[100][1000];  
  REAL Lx, Ly, Lz;
  int NP, NW, Nframes, numframe;     
  Heads();
  void readcoord(int, int);
  void AreaVoronoi(int);
  void Density(int);
  void Diffusion(int);
  REAL cpc(REAL, REAL);
  
};