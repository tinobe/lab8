// Lorenz.cxx
// Runge Kutta 4th order
// GL, 4.12.2015
//--------------------
#include <cmath>
#include <iostream>
#include <fstream>
//--------------------
void f(double* const y0, const double x);
void RKstep(double* const yn, const double* const y0, const double x, const double dx, double*, double*, double*, double*);
double inter( const double* const k1, const double* const k2, const double* const k3, const double* const k4,const double dt, const double theta, const double yn);
//--------------------
using namespace std;
//--------------------

int main(void)
{
double pmax=5;
	ofstream out("solution");
  const int dim = 2;
	double dx = 0.1,x=0;
	const double L = 100;
  	  double y0[dim];
	double yn[dim];
	double k1[dim], k2[dim], k3[dim], k4[dim];
	double interpolation,theta=1/2.0;
  //out << x << "\t" << y0[0] << "\t" << y0[1] << endl;
	for(double p0=0.1;p0<pmax;p0+=0.1){
		  y0[0]=p0;y0[1]=0;theta=1/2.0;x=0;
	while(x<=L)
	{
	  //if(test >=0 && y0[1]<=0) break;
		//test=y0[1];
		x += dx;
		RKstep(yn, y0, x, dx,k1,k2,k3,k4);
	        if(y0[1] >0 && yn[1] <0) break;
	        for(int i=0; i<dim; i++) y0[i] = yn[i];
	        out << x << "\t" << y0[0] << "\t" << y0[1] << endl;
	}
	double a=0,b=1;
	interpolation=inter(k1,k2,k3,k4,dx,theta,y0[1]);
	while(abs(interpolation)>=1e-7)
	{
	  interpolation=inter(k1,k2,k3,k4,dx,theta,y0[1]);
	  if(interpolation >=0){
	    a=theta;
	  }else
	  {
	    b=theta;
	  }
	  theta=(a+b)/2.0;
// 	  cout << theta << "\t" << interpolation << endl;
	}
	cout << p0 << " " << (x+theta*dx-dx) << endl;
	}
	
	out.close();
	return(0);
}
//-------------------
void RKstep(double* const yn, const double* const y0,
            const double x, const double dx, double* k1, double* k2, double* k3, double* k4)
{
	const int dim = 2;


  for(int i=0;i<dim; i++) k1[i] = y0[i];
	f(k1, x);

	for(int i=0;i<dim; i++) k2[i] = y0[i] + 0.5 * dx * k1[i];
  f(k2, x+0.5*dx);

	for(int i=0;i<dim; i++) k3[i] = y0[i] + 0.5 * dx * k2[i];
	f(k3, x+0.5*dx);

  for(int i=0;i<dim; i++) k4[i] = y0[i] + dx * k3[i];
	f(k4,  x+dx);

	for(int i=0;i<dim; i++)
	 yn[i] = y0[i] + 1./6.*dx*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
}
//-------------------
// Lorenz model
void f(double* const y0, const double x)
{
  double temp=y0[0];
  y0[0]=y0[1]; y0[1]=-temp/sqrt(1+temp*temp); 
}

double inter( const double* const k1, const double* const k2, const double* const k3, const double* const k4,const double dt, const double theta, const double yn){
  const double b1=theta-3*theta*theta/2.0+2*theta*theta*theta/3.0;
  const double b2=theta*theta-2*theta*theta*theta/3.0;
  const double b4=-theta*theta/2.0+2*theta*theta*theta/3.0;
  return yn+dt*(b1*k1[1]+b2*k2[1]+b2*k3[1]+b4*k4[1]);
}
