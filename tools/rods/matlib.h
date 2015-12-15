#ifndef MATLIB_H
#define MATLIB_H

#define RANDESP (rand()/(float)RAND_MAX)
#define fori(x,y) for(int i=x; i<y;i++)
#define forj(x,y) for(int j=x; j<y;j++)
#define fork(x,y) for(int k=x; k<y;k++)


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
void cross(const float *a, const float *b, float *res){
  res[0] = a[1]*b[2]-a[2]*b[1];
  res[1] = a[2]*b[0]-a[0]*b[2];
  res[2] = a[0]*b[1]-a[1]*b[0];
}

void normalize(float *v){
  float norm = sqrt(vnorm2(v,3));
  if(norm==0.0) return;
  float invnorm = 1.0/norm;
  fori(0,3) v[i] *= invnorm;
}

class unitv3{
 public:
  unitv3(float x, float y, float z){
    v[0]=x; v[1]=y; v[2]=z;
    float norm = sqrt(vnorm2(v,3));
    if(norm==0.0) norm = 1.0;
    float invnorm = 1.0/norm;
    fori(0,3) v[i] *= invnorm;
  };
  unitv3(){
    randomize();
  }
  unitv3(const unitv3 &a_v){ //Perpendicular to
    int ii = rand()%3;
    while(a_v[ii] == 0.0) ii = rand()%3;
    int j = (ii+1)%3;
    int k = (ii+2)%3;
    v[j] = RANDESP-0.5; v[k] = RANDESP-0.5;
    v[ii]=-(1.0/a_v[ii])*(a_v[j]*v[j]+a_v[k]*v[k]);
    float invnorm =1.0/ sqrt(vnorm2(v,3));
    fori(0,3) v[i] *= invnorm;
  }
  void randomize(){
    fori(0,3) v[i]=RANDESP-0.5;
    float norm = sqrt(vnorm2(v,3));
    if(norm==0.0) randomize();
    float invnorm = 1.0/norm;
    fori(0,3) v[i] *= invnorm;
  }
  const float operator[] (unsigned int i)const{return v[i];}
  float v[3];

};


#endif
