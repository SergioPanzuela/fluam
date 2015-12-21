/*
  Raul P. Pelaez 2015, initial configuration of rods generator.
  compile using make g++ -std=c++11 -Ofast -ffast-math -march=native -funroll-loops rods.cp
  read.in contains:
  Nrods Np_per_rod Nrod
  Lx Ly
  Lbox
  Kspring

  ./mesh --random for random sizes between 3x3 and NpxNrod

readin is read with:
  ifstream readin("read.in");
  readin >> Nmesh >> Np >> Nrod;
  readin >> Lx >> Ly;
  readin >> Lbox;
  readin >> k;
  readin.close();
*/

#include<iostream>
#include<stdio.h>
#include<cmath>
#include<vector>
#include<fstream>
#include<cstring>
#include"matlib.h"

#define MIN_DIST 4.0
using namespace std;

class Rod{
public: 
  Rod(int Np, float Lrod, const float *pos, const float *dir){
    this->Lrod = Lrod;
    this->Np = Np;
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
    float dist2= 0.0;
    dist2 = vnorm2(&pos[0], &other_rod.pos[0], 3);
    
    if(dist2>pow(Lrod+other_rod.Lrod,2)) return false;
    
    fori(0,Np)forj(0,other_rod.Np){
      dist2 = vnorm2(&rod[3*i], &other_rod.rod[3*j], 3);
      if(dist2<MIN_DIST) return true;
    }
      return false;
  }
  void plist(ofstream &out){
    fori(0,Np){
      forj(0,3) out<<rod[3*i+j]<<" ";
      out<<"\n";
    }
  }
  void print_bond(ofstream &out, int ip, float k, float r0){
    fori(1,Np-1){
      forj(-1,2) out<<ip+i-j<<" ";
      out<<k<<" "<<r0;
      out<<"\n";
    }
  }
  
  vector<float> rod;
  float pos[3];
  float Lrod;
  int Np;
};


class Mesh{
public:
  Mesh(){}
  Mesh(int Np, float Nrod, float Lx, float Ly, const float pos[3], const float normal[3], const float dir[3]){
    fori(0,3)this->pos[i] = pos[i];
    fori(0,3)this->normal[i] = normal[i];
    fori(0,3)this->dir[i] = dir[i];
    
    this->Lx= Lx;
    this->Ly=Ly;
    this->Nrod = Nrod;
    this->Np = Np;
    this->Nt = Np*Nrod;
    this->Nbonds = (Np-2)*Nrod+Np*(Nrod-2);

    float pdir[3];
    cross(normal, dir, pdir); normalize(pdir);
    float ppos[3];
    float dr = Lx/(float)Nrod;
    forj(0,Nrod){
      fori(0,3) ppos[i] = pos[i] + pdir[i]*dr*(0.5*(float)Nrod-j);
      mesh.push_back(Rod(Np, Ly, ppos, dir));
    }
    

  }
  float radius()const{
    return  Lx>Ly?Lx:Ly; 
  }
  bool collides(const Mesh &a_mesh){
    float dist2= 0.0;
    
    dist2 = vnorm2(&pos[0], &a_mesh.pos[0], 3);
    float r1 = radius(); 
    float r2 = a_mesh.radius();

    if(dist2>pow(r1+r2,2)) return false;
    fori(0,mesh.size())forj(0, a_mesh.mesh.size())
      if(mesh[i].collides(a_mesh.mesh[j])) return true;
      
    return false;
  }
  bool collides_with_walls(float Lbox){
    if(mesh[0].collides_with_walls(Lbox)) return true;
    else if(mesh.back().collides_with_walls(Lbox))return true;
    return false;
  }

  void print(ofstream &out){
    fori(0,mesh.size()){
      mesh[i].plist(out);
    }
  }
  void print_3bond(ofstream &out, int i0, float k){
    float r0y = Ly/(float)Np;
    //int i0 = im*Np*Nrod;
    fori(0,Nrod)
      mesh[i].print_bond(out, i0+Np*i, k, r0y);
    float r0x = Lx/(float)mesh.size();
    fori(0, Nrod-2){
      forj(0,Np){  
	fork(0,3){
	  out<<i0+Np*k+i*Np+j<<" ";
	}
	out<<k<<" "<<r0x;
	out<<"\n";	
      }
    }
  }
  void print_2bond(ofstream &out, int i0, float k){
    float L = sqrt(Lx*Lx+Ly*Ly)-1;
    out<<i0<<" "<<i0+Np*Nrod-1<<" "<<k<<" "<<L<<endl;
    out<<i0+Np-1<<" "<<i0+Np*(Nrod-1)<<" "<<k<<" "<<L<<endl;
    out<<i0+Np*(Nrod-1)<<" "<<i0+Np-1<<" "<<k<<" "<<L<<endl;
    out<<i0+Np*Nrod-1<<" "<<i0<<" "<<k<<" "<<L<<endl;
  }

  vector<Rod> mesh;
  float pos[3];
  float normal[3];
  float dir[3];
  float Lx, Ly;
  int Np, Nrod, Nbonds, Nt;
};


Mesh random_mesh(int Np, int Nrod, float r0x, float r0y, float Lbox){
  unitv3 dir, normal(dir);  
  float pos[3]; fori(0,3) pos[i]=(RANDESP-0.5)*Lbox;

  Mesh a_mesh = Mesh(Np, Nrod, r0x*Nrod, r0y*Np, pos, dir.v, normal.v);
  if(a_mesh.collides_with_walls(Lbox-1.0)) a_mesh = random_mesh(Np, Nrod, r0x, r0y, Lbox);
  return a_mesh;
}

int main(int argc, char *argv[]){
  int Nmesh, Np, Nrod;
  float Lx, Ly, k, Lbox;
  
  ifstream readin("read.in");
  readin >> Nmesh >> Np >> Nrod;
  readin >> Lx >> Ly;
  readin >> Lbox;
  readin >> k;
  readin.close();
  srand(time(NULL));

  float r0y = Lx/(float)Np;
  float r0x = Ly/(float)Nrod;

  vector<Mesh> meshes;
   bool random = false;
  
#define INOP(x,y) if(strcmp(argv[i], x)==0) { y }
  fori(0,argc){
    INOP("--random", random=true;)
  }
  #define RAND3(x) rand()%(x-2)+3 
  
  if(random)  meshes.push_back(random_mesh(RAND3(Np), RAND3(Nrod), r0x, r0y, Lbox));
  else  meshes.push_back(random_mesh(Np, Nrod, r0x, r0y, Lbox));

  int trycount = 0, count = 0;
  bool accepted = true;
  int maxcount = 10000;
  while(meshes.size()<Nmesh && trycount<maxcount){
    count++;
    if(count%1000==0 || trycount == 0){
      printf("\r\t Introducing mesh %d...    ", (int)meshes.size());
      printf("Attempt %d of %d", trycount, maxcount);
      fflush(stdout);
    }
    accepted = true;
    Mesh a_mesh;

    if(random)  a_mesh = random_mesh(RAND3(Np), RAND3(Nrod), r0x, r0y, Lbox);
    else a_mesh = random_mesh(Np, Nrod, r0x, r0y, Lbox);
    fori(0, meshes.size()){
      if(meshes[i].collides(a_mesh)){
	trycount++;
	accepted = false;
	break;
      }
    }
    if(accepted){
      meshes.push_back(a_mesh);
      trycount = 0;
    }

  }
  printf("\r\t Introducing mesh %d...    ", (int)meshes.size());
  printf("Attempt %d of %d", trycount+1, maxcount);
  fflush(stdout);
  printf("\n");
  Nmesh = meshes.size();
  if(trycount > 0)
    printf("ERROR!: Could not achieve desired density, only %d were created\n", Nmesh);

//  printf("Approximate volume fraction: %.3f\n", Nrods*Np*M_PI*(Lx+3.0/4.0)/pow(L,3));
  
  ofstream pout("meshes.pos");
  #ifndef DEBUG
      int Nt=0;
      fori(0,Nmesh)Nt += meshes[i].Nt;
      pout<<Nt<<"\n";
  #else
      pout<<"#L="<<Lbox*0.5<<";\n";
  #endif
  fori(0,Nmesh) meshes[i].print(pout);
  
  ofstream bout("meshes.3bonds");
  ofstream b2out("meshes.2bonds");
  int Nbonds= 0;
  fori(0,Nmesh) Nbonds+=meshes[i].Nbonds; 
  //  bout<<Nmesh*(Np-2)*Nrod+Nmesh*Np*(Nrod-2)<<"\n";
  bout<<Nbonds<<"\n";
  b2out<<Nmesh*4<<"\n";
  int pcount = 0;
  fori(0,Nmesh){
    meshes[i].print_3bond(bout, pcount, k);
    meshes[i].print_2bond(b2out, pcount, k);
    pcount += meshes[i].Nt;
  }
  bout.close();
  
  return 0;
}

