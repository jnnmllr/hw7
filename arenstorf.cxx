#include <fstream>
#include <iostream>
#include <cmath>

using namespace std;

void f(double* const y, const double mu, double* k);
void k(double* const y, double* k1, double* k2, double* k3, double* k4, double* k5, double* k6, double* k7, const double mu, const double dt);
void maxx(const double* const y, const double* const x, double& H, double& a, double& b);

int main(){
 double dt=1e-3;
 const double mu=0.012277471;
 const double tol = 1e-5;
 double a=0,b=0;
 double H,T,q=0.4,tau=17.065216560157;
 double x[4],y[4],k1[4],k2[4],k3[4],k4[4],k5[4],k6[4],k7[4];
 y[0]=0.994; y[1]=0; y[2]=0; y[3]=-2.00158510637908;
 x[0]=0.994; x[1]=0; x[2]=0; x[3]=-2.00158510637908;
 k1[0]=0; k1[1]=0; k1[2]=0; k1[3]=0;
 ofstream out("Data.txt");

 while(T<tau){
  out << T << " " << y[0] << " " << y[1] << endl;
  k(y,k1,k2,k3,k4,k5,k6,k7,mu,dt);
  for(int i=0;i<4;i++){y[i]+= dt*((35.0/384.0)*k1[i]+(500.0/1113.0)*k3[i]+(125.0/192.0)*k4[i]-(2187.0/6784.0)*k5[i]+(11.0/84.0)*k6[i]);}


  for(int i=0;i<4;i++){x[i]+=dt*((5179.0/57600.0)*k1[i]+(7571.0/16695.0)*k3[i]+(393.0/640.0)*k4[i]-(92097.0/339200.0)*k5[i]+(187.0/2100.0)*k6[i]+(1.0/40.0)*k7[i]);}

  maxx(y,x,H,a,b);
  dt = q*dt*pow((tol/H),(1.0/5.0));
  T += dt;

  cout << dt << " " << T << " " << a << " " << b << endl;

 }

 out.close();
 return 0;
}

void k(double* const y, double* k1, double* k2, double* k3, double* k4, double* k5, double* k6, double* k7, const double mu, const double dt){
  double kk[4];kk[0]=0;kk[1]=0;kk[2]=0;kk[3]=0;
  f(y,mu,k1);
  for(int i=0;i<4;i++){kk[i]=y[i]+dt*(1.0/5.0)*k1[i];}

  f(kk,mu,k2);
  for(int i=0;i<4;i++){kk[i]=y[i]+dt*((3.0/40.0)*k1[i]+(9.0/40.0)*k2[i]);}

  f(kk,mu,k3);
  for(int i=0;i<4;i++){kk[i]=y[i]+dt*((44.0/45.0)*k1[i]-(56.0/15.0)*k2[i]+(32.0/9.0)*k3[i]);}

  f(kk,mu,k4);
  for(int i=0;i<4;i++){kk[i]=y[i]+dt*((19372.0/6561.0)*k1[i]-(25360.0/2187.0)*k2[i]+(64448.0/6561.0)*k3[i]-(212.0/729.0)*k4[i]);}

  f(kk,mu,k5);
  for(int i=0;i<4;i++){kk[i]=y[i]+dt*((9017.0/3168.0)*k1[i]-(355.0/33.0)*k2[i]+(46732.0/5247.0)*k3[i]+(49.0/176.0)*k4[i]-(5103.0/18656.0)*k5[i]);}

  f(kk,mu,k6);
  for(int i=0;i<4;i++){kk[i]=y[i]+dt*((35.0/384.0)*k1[i]+(500.0/1113.0)*k3[i]+(125.0/192.0)*k4[i]-(2187.0/6784.0)*k5[i]+(11.0/84.0)*k6[i]);}

  f(kk,mu,k7);
}

void f(double* const y, const double mu, double* k){
double r,s;
  r = sqrt(pow((y[0]+mu),2)+pow(y[1],2));
  s = sqrt(pow((y[0]-1-mu),2)+pow(y[1],2));

  k[0]=y[2];
  k[1]=y[3];
  k[2]=y[0]+2*y[3]-(1-mu)*(y[0]+mu)/pow(r,3)-mu*(y[0]-1+mu)/pow(s,3);
  k[3]=y[1]-2*y[2]-(1-mu)*y[1]/pow(r,3)-mu*y[2]/pow(s,3);
}

void maxx(const double* const y, const double* const x, double& H, double& a, double& b){
  a = abs(y[0]-x[0]);
  b = abs(y[1]-x[1]);

  H = a;
  if(b>a){H = b;}
}
