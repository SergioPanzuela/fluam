#include<iostream>
#include<fstream>
#include<cmath>
#include<vector>
#define fori(x,y) for(int i=x; i<y;i++)
#define forj(x,y) for(int j=x; j<y;j++)
#define RANDESP  (rand()/(float)RAND_MAX)
using namespace std;

int N;
float L;
float r0;
float Rpat; //Radius size of a patchy
vector<float> ps;
vector<float> centers;

ofstream outpos;
ofstream outbond;

bool collides(float *a, float *b, float r){
  float c[2];
    c[0] = a[0]-b[0];
    c[1] = a[1]-b[1];
    float R = (c[0]*c[0]+c[1]*c[1]);
    //   cout<<R<<" "<<a[0]<<" "<<b[0]<<" "<<a[1]<<" "<<b[1]<<endl;
    if(R<=(r*r)) return true;
    return false;
}
void gen_centers(){
  srand(time(NULL));
  float pos[3];
  pos[2] = 0.0f;
  centers.resize(3*N);
  bool accept;
  int trycount =0;
  fori(0,N){
    accept = false;
    trycount = 0;
    printf("\r\t Introducing mesh %d...    ", i);
    fflush(stdout);
    do{
      trycount++;
      if(trycount%100==0){
	printf("\r\t Introducing mesh %d...    ", i);
	printf("Attempt %d", trycount);
	fflush(stdout);
      }
      accept = true;
      forj(0,2) pos[j] = (RANDESP-0.5)*(L-Rpat);
      forj(0,i){
	if(collides(&pos[0], &centers[3*j], 2*Rpat)){
          
	  accept = false;
	  break;
	} 
      }
    }while(!accept);

    forj(0,3)centers[3*i+j] = pos[j];
  }
}
typedef float pat[4*3];
void rotate_pat(const pat &x0, pat &x){
  float phi = 0;
  float theta = 0;
  float psi = 2*M_PI*RANDESP;
  //Rotate icosahedron, we use formula 4.46 from Goldstein
  for(int j=0;j<4;j++){
    x[3*j] = (cos(psi)*cos(phi)-cos(theta)*sin(phi)*sin(psi)) * x0[3*j] + 
      (cos(psi)*sin(phi)+cos(theta)*cos(phi)*sin(psi)) * x0[3*j+1] +
      (sin(psi)*sin(theta)) * x0[3*j+2] ;
    x[3*j+1] = (-sin(psi)*cos(phi)-cos(theta)*sin(phi)*cos(psi)) * x0[3*j] +
      (-sin(psi)*sin(phi)+cos(theta)*cos(phi)*cos(psi)) * x0[3*j+1] +
      (cos(psi)*sin(theta)) * x0[3*j+1];
    x[3*j+2] = (sin(theta)*sin(phi)) * x0[3*j] +
      (-sin(theta)*cos(phi)) * x0[3*j+1] +
      cos(theta) * x0[3*j+2];
  }
}

void spawn_patchys(){
  float r = 3.0f/sqrt(3);
  float h = r*sqrt(3)/2.0f;
  pat p;
  fori(0,3) p[i] = 0;
  
  p[3] = -r*0.5f;
  p[4] = -h/3.0f;
  p[5] = 0.0f;
  
  p[6] = r*0.5f;
  p[7] = -h/3.0f;
  p[8] = 0.0f;
  
  p[9] = 0.0f;
  p[10] = 2.0f*h/3.0f;
  p[11] = 0.0f;
  
  pat ptemp;
  fori(0,N){
    rotate_pat(p, ptemp);
    forj(0,12) ptemp[j] += centers[3*i+j%3];
    ps.insert(ps.end(), &ptemp[0], &ptemp[4*3]);
  }
}
void write_pat(int p){
  fori(0,4){
    outpos<<ps[4*3*p+3*i]<<" "<<ps[4*3*p+1+3*i]<<" "<<ps[4*3*p+2+3*i]<<" "<< ((i==0)?1:0 )<<endl;
  }

  outbond<<1+4*p<<" "<<0+4*p<<" "<<3+4*p<<" 1000 "<<r0<<endl;
  outbond<<3+4*p<<" "<<0+4*p<<" "<<2+4*p<<" 1000 "<<r0<<endl;
  outbond<<1+4*p<<" "<<0+4*p<<" "<<2+4*p<<" 1000 "<<r0<<endl;






}
int main(int argc, char *argv[]){
  
  N = atoi(argv[1]);
  L = stod(argv[2]);
  r0 = 1.0f;
  Rpat = 1.5f;
  //ps.resize(3*4*N);
  gen_centers();
  spawn_patchys();

  outpos.open("patchy.pos");
  outbond.open("patchy.3bond");
  outpos<<4*N<<endl;
  outbond<<3*N<<endl;
  fori(0,N){
    write_pat(i);
  }

  return 0;
}
