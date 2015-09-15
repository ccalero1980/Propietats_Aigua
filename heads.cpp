// #include <iostream>
#include "membrane.h"
// using namespace std;


 Heads::Heads() { 
  FILE *infile=0;
  infile = fopen("dump_heads_DMPC27.xyz","r");
  fscanf(infile, "%f%f%f" , &Lx, &Ly, &Lz);
  fscanf(infile, "%d", &NP);
  fscanf(infile, "%d", &NW);   
  fclose(infile);   
}




void Heads::readcoord(int nframe, int numframes)
{
  REAL x, y, z, zV;
  int n, i, num;
  
  FILE *infile=0;
  FILE *logfile=0;
  logfile = fopen("readcoord.log","w");
  infile = fopen("dump_heads_DMPC27.xyz","r");
  fscanf(infile, "%f%f%f" , &Lx, &Ly, &Lz);
  fscanf(infile, "%d", &NP);
  fscanf(infile, "%d", &NW);

  printf("Reading Coordinates\n");
  
  for(n=0; n<numframes; n++){
    for(i=0; i<NP; i++){
      xH[n][i] = 0.0;
      yH[n][i] = 0.0;
      zH[n][i] = 0.0;

    }
  }
  
  for(num=0; num<nframe+numframes; num++){   

      fscanf(infile, "%d", &numframe);
      fprintf(logfile, "%d\n", numframe);
      if(num>=nframe){printf("%d\n", num);}
 
      for(i=0; i<NP; i++){
	  fscanf(infile, "%f%f%f", &x, &y, &z);
	  if(num>=nframe){
	    
	    xH[num-nframe][i] = x;
	    yH[num-nframe][i] = y;
	    zH[num-nframe][i] = z;
	    fprintf(logfile, "%f\t%f\t%f\n", x, y, z);
	    
	  }
      }
  }  
  fclose(infile); 
  fclose(logfile); 
  
}


// ----------------------------------------------------------
//  Calculation of the distribution of the area of the Voronoi
//  regions. 
//   Monte Carlo method to compute the areas
// ----------------------------------------------------------

  void Heads::AreaVoronoi(int numframes){
  REAL rx, ry, Area[1000],Areaold[1000], histogram[1000], Areatot, AreaMax, x0, y0, dx, dy, d2, d2min, AvA1, sigmaA1, tol, tolmax, TOL, dA;
  int n, i, imin, npoints, ibin, num, Npoints, Nread, N0, N1, Nmax;
  FILE *logfile=0;
  FILE *outfile=0;
  logfile = fopen("area.log","w");
  outfile = fopen("Area_Voronoi_Histogram.dat","w");
  
  for(i=0; i<1000; i++){histogram[i] = 0.0;}
  Npoints = 10000000; 
  Nmax = 100;

  for(num=0;num<numframes;num++){
    if(num%Nmax ==0){
      N0 = num;
      N1 = numframes-num;
      if(N1>=Nmax){Nread = Nmax;}
      else if(N1<Nmax){Nread = N1;}
      readcoord(N0,Nread);
      printf("Area Voronoi Regions\n");
    }
    printf("%d\n",num);
  // ------- Initialization -------
    
    TOL = 0.01;
    srand ((unsigned int) time(NULL));
    n = 1;
    tolmax = 1.0;  
    AreaMax = 3.0*Lx*Ly/NP;
    for(i=0; i<1000; i++){
      Area[i]=0.0;
      Areaold[i]=0.0;   
      
    }
  // ------------------------------

  // ------- Monte Carlo Integration -------  

    while((tolmax>TOL)and(n<Npoints)){
    //   Generation of random point (x,y) in the Lx*Ly area:
      rx = ((REAL)rand() / (RAND_MAX));
      ry = ((REAL)rand() / (RAND_MAX));
      x0 = Lx*rx-Lx/2.0;
      y0 = Ly*ry-Ly/2.0;
      
    
      d2min = 1000000.0;
      for(i=0; i<NP; i++){
	dx = cpc(xH[num-N0][i]-x0, Lx);
	dy = cpc(yH[num-N0][i]-y0, Ly);
	d2 = dx*dx + dy*dy;
	
	if(d2<d2min){
	  imin = i;
	  d2min = d2;
	}
      }

      Area[imin] = Area[imin] + 1.0;
  //  ------------------------------------------    
  //  --- Checking quality of current result ---
      if(n%1000 == 1){
	tolmax = 0.0;
	for(i=0; i<NP; i++){
	  tol = (Area[i]-Areaold[i])/Area[i];
	  if(tol > tolmax){tolmax = tol;}
	  Areaold[i] = Area[i];
	}
      }
  //  ------------------------------------------ 
    n = n + 1;
    }
    
    
  // ---------- Writing results ----------
    npoints = n;
    
    Areatot = 0.0;
    dA = AreaMax/100.;
    for(i=0; i<NP; i++){
      Area[i]=(Area[i]/npoints)*Lx*Ly;
      Areatot = Areatot + Area[i];
      
      if(Area[i]<AreaMax){
	ibin = (int)(Area[i]/dA);
	histogram[ibin] = histogram[ibin]+1.;
      }
    }
    
    if(tolmax > 0.01){
      fprintf(logfile, "WARNING: tolerance not reached, increase Npoints\n");
      printf("WARNING: tolerance not reached, increase Npoints\n");
    }
  }
  
  // histogram
  for(n=0; n<100;n++){
      histogram[n] = histogram[n]/(numframes*NP);
      fprintf(outfile, "%8f\t%8f\n", n*dA, histogram[n]);
  }
  
  
  
// -----------------------------------
  fclose(logfile); 
  fclose(outfile);   
}


// ----------------------------------------------------------
// 
// Calculation of the density profile in z of phospholipid heads
// 
// ----------------------------------------------------------

void Heads::Density(int numframes){
  int Nmax, N0, N1, Nread, num, i, ibin;
  REAL avzH, densz[1000], dz, z;
  FILE *outfile=0;
  FILE *logfile=0;
  logfile = fopen("density.log","w");
  outfile = fopen("Heads_Density.dat","w");
    
  dz = 2.0*Lz/1000.0;
  Nmax = 100;  
  for(i=0; i<1000; i++){densz[i] = 0.0;}
  
  for(num=0;num<numframes;num++){
    if(num%Nmax ==0){
      N0 = num;
      N1 = numframes-num;
      if(N1>=Nmax){Nread = Nmax;}
      else if(N1<Nmax){Nread = N1;}
//       printf("%d\t%d\n",N0, Nread);  
      readcoord(N0,Nread);
      printf("Calculating Density\n");
    }
    
    avzH = 0.0;
    for(i=0; i<NP; i++){
      avzH = avzH + zH[num-N0][i]/NP;
    }
    
    for(i=0; i<NP; i++){
      z = ((zH[num-N0][i]-avzH)+Lz/2.0)/dz;
      if(z>0){
 	ibin = (int) z;
	densz[ibin] = densz[ibin] + 1.0;
      }
    }

  }

  for(i=0;i<1000;i++){
    densz[i] = densz[i]/(dz*Lx*Lz*numframes);
    fprintf(outfile,"%8f\t%8f\n", i*dz-Lz/2.0, densz[i]);
  }  
  fclose(outfile);    
  fclose(logfile);      
}

// ----------------------------

void Heads::Diffusion(int numframes){
  int Nmax, N0, N1, Nread, num, i, ibin;
  REAL dz, z, msd, desplz[1000], z0[1000], zold[1000], avdesplz;
  FILE *outfile=0;
  FILE *logfile=0;
  logfile = fopen("diffusion.log","w");
  outfile = fopen("Heads_Diffusion.dat","w");
    

  Nmax = 100;  
  msd = 0.0;
  for(num=0;num<numframes;num++){
    if(num%Nmax ==0){
      N0 = num;
      N1 = numframes-num;
      if(N1>=Nmax){Nread = Nmax;}
      else if(N1<Nmax){Nread = N1;}
//       printf("%d\t%d\n",N0, Nread);  
      readcoord(N0,Nread);
      dz = 2.0*Lz/1000.0;
      printf("Calculating Diffusion\n");
    }
    if(num == 0){
      for(i=0; i<NP; i++){
	desplz[i] = zH[0][i];
	z0[i] = zH[0][i];    
	zold[i] = zH[0][i];
      }
    }
    avdesplz = 0.0;
    for(i=0; i<NP; i++){
      dz = cpc(zH[num][i] - zold[i], Lz);
      desplz[i] = desplz[i] + dz;
      avdesplz = avdesplz + dz/NP;  
    }
    for(i=0; i<NP; i++){
      desplz[i] = desplz[i] -avdesplz;
      msd = msd + (desplz[i]-z0[i])*(desplz[i]-z0[i])/NP;    
      zold[i] = zH[num][i];      
    }    
    
    fprintf(outfile, "%8d\t%8f\n", num, msd);

  }

  fclose(outfile);    
  fclose(logfile);      
}

// ----------------------------

REAL Heads::cpc(REAL x, REAL L){
  if(x>L/2.0){
    x = x-L;
  }
  else if(x<-L/2.0){
   x = x+L; 
  }
  return x;
}

