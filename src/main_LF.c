#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double xL = 0., xR = 1., hx, t, Tmax=0.3, CFL = 0.9, dt, a = 1.;
#define Nx 512

double 	Uarr[Nx];

double f(double u)
{
	return 0.5*u*u;
}

void lax_friedrich()
{
    //Nx = pow(2, 6);
    hx = (xR - xL)/(Nx+1);
    dt = CFL*hx/fabs(a);

    double UarrN[Nx], t = 0.;
    int i, nt=0;


    FILE *op;
    op = fopen("output.upwind", "w");

    for(t = 0.; t<Tmax; t+=dt)
    {
        for(i = 0; i < Nx; i++)
            UarrN[i] = 0.5*((Uarr[i+1]+Uarr[i-1])-dt/hx*(f(Uarr[i+1])-f(Uarr[i-1])));

        //UarrN[0] = UarrN[Nx-1];
        for(i = 0; i < Nx; i++)
        {
            Uarr[i] = UarrN[i];
            fprintf(op,"%f %f %f\n",(double)i*hx, t, Uarr[i]);

        }
        nt++;
    }
    //}

    //for(i = 0; i < Nx; i++)
    //fprintf(op,"%f %f\n",(double)i*hx,Uarr[i]);
    //fclose(op);
    printf("%d\n", nt);

}

void init_cond()
{
    double x;
    int i;
    hx = (xR-xL)/(Nx-1);
    for( i = 0; i < Nx; i++)
    {
        x = (double)i*hx + xL;
        if(x> 0.3 && x<0.5)
            Uarr[i] = 1.;
        else
            Uarr[i] = 0.;
    }
}

int main()
{
    init_cond();
    lax_friedrich();

    return 0;
}
