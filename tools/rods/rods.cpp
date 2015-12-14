/*
  Raul P. Pelaez 2015, initial configuration of rods generator.
  compile with g++ -std=c++11 -Ofast -ffast-math -march=native -funroll-loops rods.cp
  read.in contains:
  Nrods Np_per_rod Lrod
  Lbox
  Kspring
 */

#include<iostream>
#include<stdio.h>
#include<cmath>
#include<vector>
#include<fstream>

#define MIN_DIST 4.0
#define RANDESP (rand()/(float)RAND_MAX)
#define fori(x,y) for(int i=x; i<y;i++)
#define forj(x,y) for(int j=x; j<y;j++)
using namespace std;

float vnorm2(const float *a, const float *b, int n){
  float res = 0.0;
  fori(0,n) res += pow(a[i]-b[i],2);
  return res;
}
float vnorm2(const float *a, int n){
  float res = 0.0;
  fori(0,n) res += pow(a[i],2);
  return res;
}

class Rod{
public: 
  Rod(int Np, float Lrod, const float *pos, const float *dir){
    this->Lrod = Lrod;
    rod.resize(3*Np,0);
    fori(0,3)this->pos[i] = pos[i];
    fori(0,Np)
      forj(0,3) rod[3*i+j] = pos[j]+(2.0*i/Np-1)*0.5*Lrod*dir[j];
  }

  bool collides_with_walls(float Lbox){
    //    printf("%.3f %.3f %.3f\t %.3f %.3f %.3f\n",rod[0], rod[1], rod[2], rod[rod.size()-3], rod[rod.size()-3+1], rod[rod.size()-3+2]);
    fori(0,3) if(rod[i]>0.5*Lbox || rod[i]<-0.5*Lbox) return true;
    fori(0,3) if(rod[rod.size()-3+i]>0.5*Lbox || rod[rod.size()-3+i]<-0.5*Lbox) return true;
    return false;
  }
  bool collides(const Rod &other_rod){
    int Np = rod.size()/3;
    float dist2= 0.0;
    
    dist2 = vnorm2(&pos[0], &other_rod.pos[0], 3);
    
    if(dist2>pow(Lrod+other_rod.Lrod,2)) return false;
    
    fori(0,Np)forj(0,Np){
      dist2 = vnorm2(&rod[3*i], &other_rod.rod[3*j], 3);
      if(dist2<MIN_DIST) return true;
    }
      return false;
  }
  void plist(ofstream &out){
    int Np = rod.size()/3;
    fori(0,Np){
      forj(0,3) out<<rod[3*i+j]<<" ";
      out<<"\n";
    }
  }
  void bondlist(ofstream &out, int irod, float k, float r0){
    int Np = rod.size()/3;
    fori(1,Np-1){
      forj(-1,2) out<<Np*irod+i-j<<" ";
      out<<k<<" "<<r0;
      out<<"\n";
    }

  }
  
  vector<float> rod;
  float pos[3];
  float Lrod;
};

Rod random_rod(int Np, float Lrod, float Lbox){
  float dir[3]; fori(0,3) dir[i]=RANDESP-0.5;
  float dirnorm = sqrt(vnorm2(dir,3));
  fori(0,3) dir[i] /= dirnorm;
  
  float pos[3]; fori(0,3) pos[i]=(RANDESP-0.5)*Lbox;
  Rod a_rod = Rod(Np, Lrod, pos, dir);
  if(a_rod.collides_with_walls(Lbox-1.0)) a_rod = random_rod(Np, Lrod, Lbox);
  return a_rod;
}

int main(int argc, char *argv[]){
  int Nrods, Np;
  float Lrod, L, k, r0;

  ifstream readin("read.in");
  readin >> Nrods >> Np >> Lrod;
  readin >> L;
  readin >> k;
  readin.close();
  r0 = Lrod/(float)Np;

  srand(time(NULL));
  
  vector<Rod> rods;
  rods.push_back(random_rod(Np, Lrod, L));
  int trycount = 0, count = 0;
  bool accepted = true;
  int maxcount = 100000;
  while(rods.size()<Nrods && trycount<maxcount){
    count++;
    if(count%1000==0 || trycount == 0){
      printf("\r\t Introducing rod %d...    ", (int)rods.size());
      printf("Attempt %d of %d", trycount, maxcount);
      fflush(stdout);
    }
    accepted = true;
    Rod a_rod = random_rod(Np, Lrod, L);
    fori(0, rods.size()){
      if(rods[i].collides(a_rod)){
	trycount++;
	accepted = false;
	break;
      }
    }
    if(accepted){
      rods.push_back(a_rod);
      trycount = 0;
    }

  }
  printf("\r\t Introducing rod %d...    ", (int)rods.size());
  printf("Attempt %d of %d", trycount+1, maxcount);
  fflush(stdout);
  printf("\n");
  Nrods = rods.size();
  if(trycount > 0)
    printf("ERROR!: Could not achieve desired density, only %d were created\n", Nrods);

  printf("Approximate volume fraction: %.3f\n", Nrods*M_PI*(Lrod+3.0/4.0)/pow(L,3));
  
  ofstream pout("rods.pos");
  pout<<Np*Nrods<<"\n";
  fori(0,Nrods) rods[i].plist(pout);

  ofstream bout("rods.3bonds");
  bout<<Nrods*(Np-2)<<"\n";
  fori(0,Nrods) rods[i].bondlist(bout, i, k, r0);
  bout.close();
  
  return 0;
}
