// #include <iostream>
#include "membrane.h"
// using namespace std;


 Water::Water() { 
  FILE *infile=0;
  infile = fopen("dump_water_DMPC27.xyz","r");
  fscanf(infile, "%f%f%f" , &Lx, &Ly, &Lz);
  fscanf(infile, "%d", &NP);
  fscanf(infile, "%d", &NW);   
  fclose(infile);   
}   


void Water::readcoord(int nframe, int numframes)
{
  REAL x, y, z, zV;
  int n, i, num;
  
  FILE *infile=0;
  FILE *logfile=0;
  logfile = fopen("readcoord.log","w");
  infile = fopen("dump_water_DMPC27.xyz","r");
  fscanf(infile, "%f%f%f" , &Lx, &Ly, &Lz);
   fprintf(logfile,"%f\t%f\t%f\n", Lx, Ly, Lz);
//    printf("%f\t%f\t%f\n", Lx, Ly, Lz);
  fscanf(infile, "%d", &NP);
  fprintf(logfile,"%d\n", NP);
//   printf("%d\n", NP);
  fscanf(infile, "%d", &NW);
  fprintf(logfile,"%d\n", NW);
  NmolecsW = NW/3;
//   printf("%d\t%d\n", NW, NmolecsW);
  
  for(n=0; n<numframes; n++){
    for(i=0; i<NW; i++){
      xW[n][i] = 0.0;
      yW[n][i] = 0.0;
      zW[n][i] = 0.0;
      zVW[n][i] = 0.0;

    }
  }

  for(num=0; num<nframe+numframes; num++){   

      fscanf(infile, "%d", &numframe);
      if(num>=nframe){printf("%d\n", num);}
 
      for(i=0; i<NW; i++){
	  fscanf(infile, "%f%f%f%f", &x, &y, &z, &zV);
	  if(num>=nframe){
	    
	    xW[num-nframe][i] = x;
	    yW[num-nframe][i] = y;
	    zW[num-nframe][i] = z;
	    zVW[num-nframe][i] = zV;
	    
	  }
      }
  }  
  fclose(infile); 
  
}

REAL Water::cpc(REAL x, REAL L){
  if(x>L/2.0){
    x = x-L;
  }
  else if(x<-L/2.0){
   x = x+L; 
  }
  return x;
}

 REAL Water::area() {return (Lx*Ly);}