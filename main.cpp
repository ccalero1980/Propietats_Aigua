
#include "membrane.h"

// comment

int main(){
// -----------------------
//  Declaration variables
// -----------------------

 REAL inter, zV0, zVMax, width, surf;
 int Numframes, Nsteps, elec, Nin;
 Water water;
 Heads heads;
 printf("%d\t%d\n", heads.NP, heads.NW);
 printf("%d\t%d\n", water.NP, water.NW);
// -----------------------
//  Parameters of the calculation
// -----------------------

 printf("Please input Numframes: ");
 scanf("%d", &Numframes);
 printf("Please input Nin: ");
 scanf("%d", &Nin); 
 printf("Please input zV0: ");
 scanf("%f", &zV0);
 printf("Please input zVMax: ");
 scanf("%f", &zVMax);
 printf("Please input width: ");
 scanf("%f", &width);
 
 printf("LJ (0) or Electrostatics (1)?: ");
 scanf("%d", &elec);


 Nsteps = int((zVMax-zV0)/width);

 
// -----------------------
//  Execution
// ----------------------- 

//  LIPID HEADS:
 
//  heads.AreaVoronoi(Numframes);
//  heads.Density(Numframes);
//  heads.Diffusion(Numframes);
 
 
//  INTERFACIAL WATER:
 

   if(elec==0){
    water.InterfacialLJ(Nin, Numframes, zV0, Nsteps, width);
   }
   else{
     water.InterfacialC(Numframes, zV0, Nsteps, width);
    }

 return 0;
  
}
