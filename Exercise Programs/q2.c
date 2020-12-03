#include <stdio.h>
#include <math.h>

int main()
{
    double rk4x(double F, double w,double y0, double z0, double h, int n) {
        double r, t0=0, c;
        double f(double t, double y, double z) {
            c = cos(w*t);
            r = F*c - (float)0.25*z + (float)0.5*y*(1-y*y);
            return r;
        }
        double k1y, k1z, k2y, k2z, k3y, k3z, k4y, k4z;
        double yi = y0;
        double zi = z0;
        double ti = t0;
        int i;
        for(i=0;i<n;i++){
            k1y = h*zi;
            k1z = h*f(ti,yi,zi);
            k2y = h*(zi+k1z/(float)2);
            k2z = h*f(ti+h/(float)2, yi+k1y/(float)2, zi+k1z/(float)2);
            k3y = h*(zi+k2z/(float)2);
            k3z = h*f(ti+h/(float)2, yi+k2y/(float)2, zi+k2z/(float)2);
            k4y = h*(zi+k3z);
            k4z = h*f(ti+h/(float)2, yi+k3y/(float)2, zi+k3z/(float)2);

            yi = yi + (k1y+2*k2y+2*k3y+k4y)/6;
            zi = zi + (k1z+2*k2z+2*k3z+k4z)/6;
            ti = ti + h;
        }
        return yi;
    }

    double rk4v(double F, double w,double y0, double z0, double h, int n) {
        double r, t0=0, c;
        double f(double t, double y, double z) {
            c = cos(w*t);
            r = F*c - (float)0.25*z + (float)0.5*y*(1-y*y);
            return r;
        }
        double k1y, k1z, k2y, k2z, k3y, k3z, k4y, k4z;
        double yi = y0;
        double zi = z0;
        double ti = t0;
        int i;
        for(i=0;i<n;i++){
            k1y = h*zi;
            k1z = h*f(ti,yi,zi);
            k2y = h*(zi+k1z/(float)2);
            k2z = h*f(ti+h/(float)2, yi+k1y/(float)2, zi+k1z/(float)2);
            k3y = h*(zi+k2z/(float)2);
            k3z = h*f(ti+h/(float)2, yi+k2y/(float)2, zi+k2z/(float)2);
            k4y = h*(zi+k3z);
            k4z = h*f(ti+h/(float)2, yi+k3y/(float)2, zi+k3z/(float)2);

            yi = yi + (k1y+2*k2y+2*k3y+k4y)/6;
            zi = zi + (k1z+2*k2z+2*k3z+k4z)/6;
            ti = ti + h;
        }
        return zi;
    }



    double F, xa, va, H, w=1;
    int trans=200000, i;
    float pi=3.14159265;

    FILE * fp;
    fp = fopen ("./q2.csv","w");
    fprintf (fp, "xa,Fa\n");
    for(F=0;F<=0.7;F+=0.0005){
        xa = 1;
        va = -1;
        H = 0.01*2*pi/w;
        xa = rk4x(F, w, xa, va, H, trans);
        va = rk4v(F, w, xa, va, H, trans);
        H = 0.001*2*pi/w;
        for(i=0;i<=100;i++){
            xa = rk4x(F, w, xa, va, H, 1000);
            va = rk4v(F, w, xa, va, H, 1000);
            fprintf (fp, "%.2f, %.2f\n", xa, F);
        }
    }

    fclose (fp);

    return 0;
}

