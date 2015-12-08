#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <sstream>

using namespace std;

void function(double*, int, int, double, double, double, double);

int main(){

  double x,y,z;
  
  const int N = 3;
  double k1[N];
  double k2[N];
  double k3[N];
  double k4[N];
  
  const int tmin = 0;
  const int tmax = 100;
  double dtmin = 0.0001;
  double dt;
  double Nit;
    
  const int M = 5;
  string fname = "hw6_output_dt_";
  stringstream s;
      
  const int a    = 10;
  const int b    = 28;
  const double c = 8.0/3.0;

  for(int j=0; j<M; j++){
    
      //initial values
    x = 1;
    y = 1;
    z = 1;
    
    double t = tmin;
    double dt = dtmin * pow(10,j);
    const double Nit = (tmax-tmin)/dt;
    
    s.str("");
    s << fname << dt << ".dat"; 
    string filename = s.str();
    ofstream out (filename.c_str());
    
    out << t << "\t" << x << "\t" << y << "\t" << z << endl;
    
    for(int i=0; i<Nit; i++){
      
      t += dt;
      
      function(k1,a,b,c,x,y,z);
      
      
      double xk2 = x + 0.5*dt * k1[0];
      double yk2 = y + 0.5*dt * k1[1];
      double zk2 = z + 0.5*dt * k1[2];
      
      function(k2,a,b,c,xk2,yk2,zk2);
      
      
      double xk3 = x + 0.5*dt * k2[0];
      double yk3 = y + 0.5*dt * k2[1];
      double zk3 = z + 0.5*dt * k2[2];
      
      function(k3,a,b,c,xk3,yk3,zk3);
      
      
      double xk4 = x + dt * k3[0];
      double yk4 = y + dt * k3[1];
      double zk4 = z + dt * k3[2];
      
      function(k4,a,b,c,xk4,yk4,zk4);
            
    
      x += dt/6.0 * (k1[0] + 2*k2[0] + 2*k3[0] + k4[0]);
      y += dt/6.0 * (k1[1] + 2*k2[1] + 2*k3[1] + k4[1]);
      z += dt/6.0 * (k1[2] + 2*k2[2] + 2*k3[2] + k4[2]);
      
      out << t << "\t" << x << "\t" << y << "\t" << z << endl;
      }
    out.close();
  }
  return 0;
}

void function(double* k, int a, int b, double c, double x, double y, double z){

  k[0] = a * (y - x);
  k[1] = x * (b - z) - y;
  k[2] = x * y - c * z;
}