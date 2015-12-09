// Lorenz.cxx
// Runge Kutta 4th order
// GL, 4.12.2015
//--------------------
#include <cmath>
#include <iostream>
#include <fstream>
//--------------------
void f(double* const y0, const double x);
void RKstep(double* const yn, const double* const y0, const double x, const double dx);
//--------------------
using namespace std;
//--------------------

int main(void)
{
	ofstream out("solution");
  const int dim = 3;
	double dx = 0.001,x=0;
	const double L = 100;
  double y0[dim] = {1.01, 1.0, 1.0};
	double yn[dim];

  out << x << "\t" << y0[0] << "\t" << y0[1] << "\t" << y0[2] << endl;
	while(x<=L)
	{
		x += dx;
		RKstep(yn, y0, x, dx);
    for(int i=0; i<dim; i++) y0[i] = yn[i];
		out << x << "\t" << y0[0] << "\t" << y0[1] << "\t" << y0[2] << endl;
	}
	out.close();
	return(0);
}
//-------------------
void RKstep(double* const yn, const double* const y0,
            const double x, const double dx)
{
	const int dim = 3;
	double k1[dim], k2[dim], k3[dim], k4[dim];

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
	const double a = 10;
	const double b = 28;
	const double c = 8.0/3.0;
	double y[3] = { y0[0], y0[1], y0[2] };

  y0[0] = a*(y[1] - y[0]);
	y0[1] = y[0]*(b - y[2]) - y[1];
	y0[2] = y[0]*y[1] - c*y[2];
}
