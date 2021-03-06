// Filename: forceBondedGPU.cu
//
// Copyright (c) 2010-2015, Florencio Balboa Usabiaga
//
// This file is part of Fluam
//
// Fluam is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Fluam is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Fluam. If not, see <http://www.gnu.org/licenses/>.

__device__ void harmonic_spring(double r0, double kSpring, double rx, double ry, double rz, double x,
			   double y, double z, double &fx, double &fy, double &fz){
  double r = sqrt( (x-rx)*(x-rx) + (y-ry)*(y-ry) + (z-rz)*(z-rz) );
  if(r>0){//If r=0 -> f=0
    fx += -kSpring * (1 - r0/r) * (rx - x);
    fy += -kSpring * (1 - r0/r) * (ry - y);
    fz += -kSpring * (1 - r0/r) * (rz - z);
  }
}

__device__ void fene_spring(double rinf, double kSpring, double rx, double ry, double rz, double x, double y,
		       double z,double &fx, double &fy, double &fz){
  double r2 = ( (x-rx)*(x-rx) + (y-ry)*(y-ry) + (z-rz)*(z-rz) );
  if(r2>0){//If r=0 -> f=0
    double rinf2 = rinf*rinf;
    if(r2>=rinf2) printf("ERROR IN FENE FORCE, R > RINF!!!!, r=%.3f   rinf= %.3f\n", sqrt(r2), rinf);
    double A = -kSpring/(1 - r2/rinf2); 
    fx += A * (rx - x);
    fy += A * (ry - y);
    fz += A * (rz - z);
  }
  

}


__device__ void forceBondedParticleParticleGPU(const int i,
					       double& fx, //Pass by reference
					       double& fy,
					       double& fz,
					       const double rx,
					       const double ry,
					       const double rz,
					       const bondedForcesVariables* bFV){


 double x, y, z;
  double r, r0;
  double kSpring;
  int index;

  int stype;

  //Particle-Particle Force
  int nBonds = bFV->bondsParticleParticleGPU[i];
  int offset = bFV->bondsParticleParticleOffsetGPU[i];
  
  
  for(int j=0;j<nBonds;j++){

    index = bFV->bondsIndexParticleParticleGPU[offset + j];

    //if(i==0) index=1;
    //if(i==1) index=0;


    //Particle bonded coordinates
    x = fetch_double(texrxboundaryGPU,nboundaryGPU+index);
    y = fetch_double(texryboundaryGPU,nboundaryGPU+index);
    z = fetch_double(texrzboundaryGPU,nboundaryGPU+index);
    
    //Spring type CAUTION TEMPORAL WORK AROUND 
    stype = (int) bFV->r0ParticleParticleGPU[offset+j];
	
    // printf(" %d %d %d\n", i, index, stype);   
 
    //Spring constant
    kSpring = bFV->kSpringParticleParticleGPU[offset+j];

    if(stype==0){  // harmonic potential  TODO
      r0 = 3.200  ;  //TODO
      harmonic_spring(r0,kSpring,rx,ry,rz,x,y,z,fx,fy,fz);
    }
    else if(stype==1){  // FENE potential  TODO
      r0 = 1.200 ;   //TODO
      fene_spring(r0,kSpring,rx,ry,rz,x,y,z,fx,fy,fz);

    }

  }











  //Particle-FixedPoint Force
  nBonds = bFV->bondsParticleFixedPointGPU[i];
  offset = bFV->bondsParticleFixedPointOffsetGPU[i];
  
  
  for(int j=0;j<nBonds;j++){
    

    //Fixed point coordinates
    x = bFV->rxFixedPointGPU[offset+j];
    y = bFV->ryFixedPointGPU[offset+j];
    z = bFV->rzFixedPointGPU[offset+j];
    
    //Equilibrium distance 
    r0 = bFV->r0ParticleFixedPointGPU[offset+j];

    //Spring constant
    kSpring = bFV->kSpringParticleFixedPointGPU[offset+j];
    
    if(r0==0){
      fx += -kSpring * (rx - x);
      fy += -kSpring * (ry - y);
      fz += -kSpring * (rz - z);
    }  
    else{     //If r0!=0 calculate particle particle distance
      r = sqrt( (x-rx)*(x-rx) + (y-ry)*(y-ry) + (z-rz)*(z-rz) );
      if(r>0){//If r=0 -> f=0
	fx += -kSpring * (1 - r0/r) * (rx - x);
	fy += -kSpring * (1 - r0/r) * (ry - y);
	fz += -kSpring * (1 - r0/r) * (rz - z);
      }
    }
  }

    
}



__device__ void forceBondedThreeParticleGPU(const int i,
					    double& fx, //Pass by reference
					    double& fy,
					    double& fz,
					    const double rx,
					    const double ry,
					    const double rz,
					    const threeParticleBondsVariables* tPBV){	


  double ri[3], rj[3], rl[3];
  double rij[3], ril[3];
  double rij2, ril2, rij1, ril1;
  double a1, a2, a3;
  double r, r0;
  double kSpring;
  double ampli;
  int p1, p2, p3;

  //Particle-Particle Force
  int nBonds = tPBV->Nbonds[i];
  int offset = tPBV->cumulative_index[i];
  int bond;
  for(int j=0;j<nBonds;j++){
    bond = tPBV->isinbonds[offset+j];
    
    p1 = tPBV->bondList[3*bond];
    p2 = tPBV->bondList[3*bond+1];
    p3 = tPBV->bondList[3*bond+2];

    //Particle bonded coordinates
    rj[0] = fetch_double(texrxboundaryGPU,nboundaryGPU+p1);
    rj[1] = fetch_double(texryboundaryGPU,nboundaryGPU+p1);
    rj[2] = fetch_double(texrzboundaryGPU,nboundaryGPU+p1);
    //Particle bonded coordinates
    ri[0] = fetch_double(texrxboundaryGPU,nboundaryGPU+p2);
    ri[1] = fetch_double(texryboundaryGPU,nboundaryGPU+p2);
    ri[2] = fetch_double(texrzboundaryGPU,nboundaryGPU+p2);
    //Particle bonded coordinates
    rl[0] = fetch_double(texrxboundaryGPU,nboundaryGPU+p3);
    rl[1] = fetch_double(texryboundaryGPU,nboundaryGPU+p3);
    rl[2] = fetch_double(texrzboundaryGPU,nboundaryGPU+p3);
    
    /*
    printf("i= %d     current bond: %d %d %d\n   pos\n p1: %.3f %.3f %.3f\n p2: %.3f %.3f %.3f\n p3: %.3f %.3f %.3f\n     ",
	   i, p1, p2, p3, rj[0], rj[1], rj[2], ri[0], ri[1], ri[2], rl[0], rl[1], rl[2]);
    */
    rij[0] = ri[0]-rj[0];    rij[1] = ri[1]-rj[1];    rij[2] = ri[2]-rj[2];
    rij2 = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];
    rij1 = sqrt(rij2);

    ril[0] = ri[0]-rl[0];    ril[1] = ri[1]-rl[1];    ril[2] = ri[2]-rl[2];
    ril2 = ril[0]*ril[0] + ril[1]*ril[1] + ril[2]*ril[2];
    ril1 = sqrt(ril2);

    a1 = rij[0]*ril[0] + rij[1]*ril[1] + rij[2]*ril[2];
    a2 = rij1*ril1;
    a3 = a1/a2;             //a3 = cos (teta) = rij*ril / mod(rij)*mod(ril)                           


    //Spring constant
    kSpring = tPBV->kSprings[bond];
    //Equilibrium distance 
    r0 = tPBV->r0Springs[bond];  
    /*    
    if(a3<=-0.9999){
      ampli = -kSpring*sqrt(2.0/(1-a3));
    }
    */
     ampli = kSpring * (acos(a3)-3.1415)/sqrt(1-a3*a3);
    

    //p1 is j, p2 is i, p3 is l
    
    if(i==p1){
      fx += ampli * (a3*rij[0]/rij2 - ril[0]/a2);
      fy += ampli * (a3*rij[1]/rij2 - ril[1]/a2);
      fz += ampli * (a3*rij[2]/rij2 - ril[2]/a2);
      
      r = rij1;
      fx += -kSpring * (1 - r0/r) * (rj[0] - ri[0]);
      fy += -kSpring * (1 - r0/r) * (rj[1] - ri[1]);
      fz += -kSpring * (1 - r0/r) * (rj[2] - ri[2]);
      
    }
    else if(i==p2){
      
      fx += ampli * (-a3*(rij[0]/rij2 + ril[0]/ril2) + (1/a2)*(rij[0] + ril[0]));
      fy += ampli * (-a3*(rij[1]/rij2 + ril[1]/ril2) + (1/a2)*(rij[1] + ril[1]));
      fz += ampli * (-a3*(rij[2]/rij2 + ril[2]/ril2) + (1/a2)*(rij[2] + ril[2]));
      
      //First spring
      r = rij1;

      fx += -kSpring * (1 - r0/r) * (ri[0] - rj[0]);
      fy += -kSpring * (1 - r0/r) * (ri[1] - rj[1]);
      fz += -kSpring * (1 - r0/r) * (ri[2] - rj[2]);
      
      //Second spring
      r = ril1;

      fx += -kSpring * (1 - r0/r) * (ri[0] - rl[0]);
      fy += -kSpring * (1 - r0/r) * (ri[1] - rl[1]);
      fz += -kSpring * (1 - r0/r) * (ri[2] - rl[2]);
      
    }
    else if(i==p3){
      fx += ampli * (a3*ril[0]/ril2 - rij[0]/a2);
      fy += ampli * (a3*ril[1]/ril2 - rij[1]/a2);
      fz += ampli * (a3*ril[2]/ril2 - rij[2]/a2);
      
      r = ril1;

      fx += -kSpring * (1 - r0/r) * (rl[0] - ri[0]);
      fy += -kSpring * (1 - r0/r) * (rl[1] - ri[1]);
      fz += -kSpring * (1 - r0/r) * (rl[2] - ri[2]);
            
    }

    //printf("i=%d\n %.3f %.3f %.3f %.3f %.3f %.3f %.3f\n",i, a1, a2, a3, ampli, fx, fy, fz);
    
  }
}

