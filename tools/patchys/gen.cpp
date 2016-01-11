#include<iostream>
#include<fstream>
#include<cmath>
#include<vector>
#define fori(x,y) for(int i=x; i<y;i++)
#define RANDESP  (rand()/(float)RAND_MAX)
using namespace std;

int N;
float L;
vector<float> ps;
ofstream outpos;
ofstream outbond;

typedef float pat[4*3];

void spawn_patchy(){

  float r = 1.0f;

  float h = r*sqrt(3)/2.0f;
  
  pat p;
  
  fori(0,3) p[i] = 0;
  
  p[3] = -r*0.5f;
  p[4] = -h*0.5f;
  p[5] = 0.0f;
  
  p[6] = r*0.5f;
  p[7] = -h*0.5f;
  p[8] = 0.0f;
  
  p[9] = 0.0f;
  p[10] = 0.5f*h;
  p[11] = 0.0f;

  float pos[3];
  fori(0,3) pos[i] = (RANDESP-0.5)*(L-r);
  
  //fori(0,12){
    //p[i] += pos[i%3];
    //  }

  ps.insert(ps.end(), &p[0], &p[4*3]);

}
void write_pat(int p){
  fori(0,4){
    outpos<<ps[4*3*p+3*i]<<" "<<ps[4*3*p+1+3*i]<<" "<<ps[4*3*p+2+3*i]<<" "<< ((i==0)?1:0 )<<endl;
  }

  outbond<<1+4*p<<" "<<0+4*p<<" "<<3+4*p<<" 1000 1"<<endl;
  outbond<<3+4*p<<" "<<0+4*p<<" "<<2+4*p<<" 1000 1"<<endl;
  //outbond<<1+4*p<<" "<<0+4*p<<" "<<2+4*p<<" 1000 1"<<endl;






}
int main(int argc, char *argv[]){
  
  N = 1;
  L = 10.0f;
  //ps.resize(3*4*N);

  fori(0,N){
    spawn_patchy();
  }
  outpos.open("patchy.pos");
  outbond.open("patchy.3bond");
  outpos<<4*N<<endl;
  outbond<<2*N<<endl;
  fori(0,N){
    write_pat(i);



  }

  return 0;
}
