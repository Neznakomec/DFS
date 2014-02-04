#include "math.h"

#define Pi  3.14159265358979323846

double u(double x, double y, double z)
{
    double f=0.0;
    if ((x<5.0) && (y<5.0) && (z<5.0))
    {
       f=sin(x*Pi/5.0)+sin(y*Pi/5.0)+sin(z*Pi/5.0);
    }
    return f;
}

double f(double x, double y, double z)
{
    double fxyz=0.0;
    if (((x>=5.0) && (x<=5.2)) && ((y>=5.0) && (y<=5.2)) && ((z>=5.0) && (z<=5.2)))
        fxyz=4.0;
    return fxyz;
}