#include "membrane.h"

void Water::InterfacialLJ(int Nin, int Nframes, REAL zV0, int Nsteps, REAL width){
  int n, m, i,k, ncasos, N0, N1, Nmax;
  REAL e1,zV, E1[100000], Etotal[1000], AvE, sigmaE;
  char strNframes[15], strNin[15], noutfile1[20], noutfile2[20];
  FILE *outfile1=0;
  FILE *outfile2=0;
  FILE *logfile=0;
  
  sprintf(strNframes, "%d", Nframes);
  sprintf(strNin, "%d", Nin);
  sprintf(noutfile1, "Energies_LJ_%s_%s.dat", strNin, strNframes);
  sprintf(noutfile2, "EnergiesAv_LJ_%s_%s.dat", strNin, strNframes);


  outfile1 = fopen(noutfile1,"w");
  outfile2 = fopen(noutfile2,"w");
  logfile = fopen("logfile.log","w");

  Nmax = 20;
  for(m=0;m<Nframes;m++){Etotal[m] = 0.0;}
    
  for(n=Nin;n<Nframes+Nin; n++){
      if((n-Nin)%Nmax==0){
	N0 = n;
	N1 = Nframes+Nin-n;
	if(N1>=Nmax){readcoord(n, Nmax);}
	else if(N1<Nmax){readcoord(n, N1);}
      }
    Etotal[n-Nin] = InteractionLJ(n-N0,-1);
    printf("%8f\n" , Etotal[n-Nin]);
  }
  

  

  
  for(k=0; k<Nsteps; k++){
    for(m=0;m<100000;m++){E1[m] = 0.0;}
    ncasos = 0;
    zV = zV0 + k*width;
    
    fprintf(outfile1,"\n%8f\n" , zV);
    fprintf(outfile2,"%8f\t" , zV);
    printf("%8f\n" , zV);
    
    for(n=Nin;n<Nframes+Nin; n++){
      if((n-Nin)%Nmax==0){
	N0 = n;
	N1 = Nframes+Nin-N0;
	if(N1>=Nmax){readcoord(n, Nmax);}
	else if(N1<Nmax){readcoord(n, N1);}
      }
      for(i=0; i<NmolecsW; i++){

	if((zVW[n-N0][3*i]>zV)&&(zVW[n-N0][3*i] < zV+width)){
	  e1 = Etotal[n-Nin] - InteractionLJ(n-N0,3*i);
	  
	  fprintf(outfile1,"%8f\n" , e1);
	  
	  E1[ncasos] = e1;
	  ncasos = ncasos + 1;
	  
	}
      }
    }
    AvE = 0.0;
    sigmaE = 0.0;
    stats(E1, ncasos, &AvE, &sigmaE);
    fprintf(outfile2,"%8f\t%8f\n" , AvE, sigmaE);
    printf("%8f\t%8f\n" , AvE, sigmaE);
    
  }
 fclose(outfile1); 
 fclose(outfile2);
 fclose(logfile);
}

REAL Water::InteractionLJ(int nframe, int index){
//    REAL cpc(REAL x, REAL L);
  int i, j;
  REAL dx, dy, dz, r2, aux2, aux6, aux12, ULJ, coeff1, coeff2, Ucorr,R2cutoff;
  coeff1 = 25.80524056;
  coeff2 = 25246.059;
  R2cutoff = 144.2401; //12.01*12.01
  ULJ = 0.0;
  Ucorr = -0.114486*2*3.14159*NmolecsW*NmolecsW/(Lx*Ly*Lz);
  
//   printf("%d\t%f\n", nframe, Ucorr);
  
    for(i=0; i<NmolecsW-1; i++){
      for(j=i+1; j<NmolecsW; j++){
	if((3*i!=index)&&(3*j!=index)){
	  dx = cpc(xW[nframe][3*i] - xW[nframe][3*j],Lx);
	  dy = cpc(yW[nframe][3*i] - yW[nframe][3*j],Ly);
	  dz = cpc(zW[nframe][3*i] - zW[nframe][3*j],Lz);
	  
	  r2 = dx*dx+dy*dy+dz*dz;
//  	  printf("%f\t", r2);
	  
	  if(r2<R2cutoff){
// 	    if(r2<1.0){printf("%f\t", r2);}
	    aux2 = 1.0/r2;
	    aux6 = aux2*aux2*aux2;
	    aux12 = aux6*aux6;
	    ULJ = ULJ  + coeff2*aux12 - coeff1*aux6;
	  }
	}
      }
    }
  ULJ = ULJ + Ucorr;
  return ULJ;
}




void Water::InterfacialC(int Nframes, REAL zV0, int Nsteps, REAL width){
  int n, m, i,k, ncasos, N0, N1,Nin, Nmax;
  REAL e1,zV, E1[100000], Etotal[1000], AvE, sigmaE;
  FILE *outfile2=0;
  FILE *outfile3=0;
  FILE *logfile=0;
  outfile2 = fopen("Energies_Electrostatics.dat","w");
  outfile3 = fopen("EnergiesAv_Electrostatics.dat","w");
  logfile = fopen("logfile.log","w");


  Nin = 0;
  Nmax = 20;
  for(m=0;m<Nframes;m++){Etotal[m] = 0.0;}
  for(n=Nin;n<Nframes+Nin; n++){
      if((n-Nin)%Nmax==0){
	N0 = n;
	N1 = Nframes+Nin-n;
	if(N1>=Nmax){readcoord(n, Nmax);}
	else if(N1<Nmax){readcoord(n, N1);}
      }
    Etotal[n-Nin] = InteractionC(n-N0,-1);
    printf("%8f\n" , Etotal[n-Nin]);
  }
  

  
  for(k=0; k<Nsteps; k++){
    for(m=0;m<100000;m++){E1[m] = 0.0;}
    ncasos = 0;
    zV = zV0 + k*width;
    
    fprintf(outfile2,"\n%8f\n" , zV);
    fprintf(outfile3,"%8f\t" , zV);
    printf("%8f\n" , zV);
    
    for(n=Nin;n<Nframes+Nin; n++){
      if((n-Nin)%Nmax==0){
	N0 = n;
	N1 = Nframes+Nin-N0;
	if(N1>=Nmax){readcoord(n, Nmax);}
	else if(N1<Nmax){readcoord(n, N1);}
      }
      for(i=0; i<NmolecsW; i++){

	if((zVW[n-N0][3*i]>zV)&&(zVW[n-N0][3*i] < zV+width)){
	  e1 = Etotal[n-Nin] - InteractionC(n-N0,3*i);
	  
	  fprintf(outfile2,"%8f\n" , e1);
	  
	  E1[ncasos] = e1;
	  ncasos = ncasos + 1;
	  
	}
      }
    }
    AvE = 0.0;
    sigmaE = 0.0;
    stats(E1, ncasos, &AvE, &sigmaE);
    fprintf(outfile3,"%8f\t%8f\n" , AvE, sigmaE);
    printf("%8f\t%8f\n" , AvE, sigmaE);
    
  }
 fclose(outfile2); 
 fclose(outfile3);
 fclose(logfile);
}




REAL Water::InteractionC(int nframe, int index){
  int i, j, ir, jr, nmax, nsqmax, ntot, nx, ny, nz, nsq;
  REAL dx, dy, dz, r2, R2cutoff, ke, eps0, V, alpha, alphasq, UC, UCreal, qi, qj,r;
  REAL UCself, qOw, qH, Uk, kx, ky, kz, ksq, sk, skcos, sksin, UCk, pi;
  
  pi = 3.141592654;  
  R2cutoff = 144.2401; //12.01*12.01
  ke = 332.0636; // kcal*A/(mol*e^2)
  eps0 = 1./(ke*4.*pi);
  V = Lx*Ly*Lz;
  qOw = -0.834;
  qH = 0.417;

  
  alpha=0.22643442364279220;
  alphasq = alpha*alpha;
  nmax = 10;
  nsqmax = 3*nmax*nmax;
  
  UC = 0.0;
  
  
// REAL SPACE TERM  
  UCreal = 0.0;
  
    for(i=0; i<NmolecsW-1; i++){
      for(j=i+1; j<NmolecsW; j++){
	if((3*i!=index)&&(3*j!=index)){
	  for(ir=0;ir<3; ir++){
	    for(jr=0;jr<3; jr++){
	      dx = cpc(xW[nframe][3*i+ir] - xW[nframe][3*j+jr],Lx);
	      dy = cpc(yW[nframe][3*i+ir] - yW[nframe][3*j+jr],Ly);
	      dz = cpc(zW[nframe][3*i+ir] - zW[nframe][3*j+jr],Lz);
	      r2 = dx*dx+dy*dy+dz*dz;
	      if(r2<R2cutoff){
		if(ir==0){qi = qOw;}
		else{qi = qH;}
		if(jr==0){qj = qOw;}
		else{qj = qH;}
		r = sqrtf(r2);
		UCreal = UCreal + qi*qj*(1-erf(alpha*r))/r;
	      }
	    }
	  }
	}
      }
    }
  UCreal = ke*UCreal;

// SELF INTERACTION TERM
  UCself = 0.0;
  for(i=0; i<NmolecsW; i++){
    if(3*i!=index){
	UCself = UCself + qOw*qOw*alpha/sqrtf(pi); 
	UCself = UCself + 2*qH*qH*alpha/sqrtf(pi);
	for(ir=0;ir<2; ir++){
	  for(jr=ir+1; jr<3; jr++){
	      dx = cpc(xW[nframe][3*i+ir] - xW[nframe][3*i+jr],Lx);
	      dy = cpc(yW[nframe][3*i+ir] - yW[nframe][3*i+jr],Ly);
	      dz = cpc(zW[nframe][3*i+ir] - zW[nframe][3*i+jr],Lz);
	      r2 = dx*dx+dy*dy+dz*dz;
	      r = sqrtf(r2);
	      if(ir==0){qi = qOw;}
	      else{qi = qH;}
	      UCself = UCself +qi*qH*erf(alpha*r)/r; 
	  }
	  
	}
    }
  }
  UCself = ke*UCself;
  
  
// RECIPROCAL TERM
  UCk = 0.0;
  ntot = 0;
  for(nx=0;nx<nmax+1; nx++){
      for(ny=-nmax;ny<nmax+1;ny++){
	  for(nz=-nmax;nz<nmax+1;nz++){
	      nsq = nx*nx + ny*ny + nz*nz;
	      if((nsq<nsqmax)&&(nsq!=0)){
		  ntot = ntot + 1;
		  kx = nx*2.0*pi/Lx;
		  ky = ny*2.0*pi/Ly;
		  kz = nz*2.0*pi/Lz;
		  ksq = kx*kx + ky*ky + kz*kz;
		  
		  skcos = 0.0;
		  sksin = 0.0;
		  for(i=0; i<NmolecsW; i++){
		    if(3*i!=index){
			for(ir=0;ir<3; ir++){
			    if(ir==0){qi = qOw;}
			    else{qi = qH;}
			    skcos = skcos + qi*cos(kx*xW[nframe][3*i+ir] +  ky*yW[nframe][3*i+ir] + kz*zW[nframe][3*i+ir]);
			    sksin = sksin + qi*sin(kx*xW[nframe][3*i+ir] +  ky*yW[nframe][3*i+ir] + kz*zW[nframe][3*i+ir]);
		
			}
		    }
		  }
		  sk = skcos*skcos + sksin*sksin;
		  if(nx==0){UCk = UCk + sk*exp(-ksq/(4*alphasq))/ksq; }
		  else{UCk = UCk + 2.0*sk*exp(-ksq/(4*alphasq))/ksq;}
		  
	      }
	  }
      }
  }
  UCk = UCk/(2.0*V*eps0);
  
//   TOTAL ELECTROSTATIC ENERGY:
  printf("%8f\t%8f\t%8f\n" , UCreal, UCself, UCk);
  UC = UCk + UCreal - UCself;
  return UC;  
}


void Water::stats(REAL E1[1000], int ncasos, REAL* AvE,REAL* sigmaE){
  REAL AvE1, sigma2E1, sigmaE1;
  int n;
  AvE1 = 0.0;
  sigma2E1 = 0.0;
  for(n=0;n<ncasos;n++){
    AvE1 = AvE1 + E1[n]/ncasos;
  }
  for(n=0;n<ncasos;n++){
    sigma2E1 = sigma2E1 + (E1[n]-AvE1)*(E1[n]-AvE1)/ncasos;
  }
  sigmaE1 = sqrtf(sigma2E1);
  *AvE = AvE1;
  *sigmaE = sigmaE1;
}

