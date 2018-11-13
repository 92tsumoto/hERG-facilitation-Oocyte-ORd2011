// runge-kutta.c
#include "syspara.h"

void runge(int n,double h,double x[],double t)
{
    double k1[NN],k2[NN],k3[NN],k4[NN];
    double xtemp[NN];
    int i;  
	
    function(x,k1,t);
    for (i=0; i<NN; i++) xtemp[i]=x[i]+k1[i]/2.0*h;
    t+=h/2.0;

    function(xtemp,k2,t);
    for (i=0; i<NN; i++) xtemp[i]=x[i]+k2[i]/2.0*h;

    function(xtemp,k3,t);
    for (i=0; i<NN; i++) xtemp[i]=x[i]+k3[i]*h;
    t+=h/2.0;

    function(xtemp,k4,t);
    for (i=0; i<NN; i++) x[i]+=(k1[i]+2.0*(k2[i]+k3[i])+k4[i])/6.0*h;
}
