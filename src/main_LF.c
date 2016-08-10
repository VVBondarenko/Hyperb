#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double xL = 0., xR = 1., hx, t, Tmax=0.3, CFL = 0.9, dt, a = 1.;
#define Nx 256

double 	Uarr[Nx];
double 	Uexact[Nx];

double f(double u)
{
    return 0.5*u*u;
}

void exact_solution(double lL, double lR, double tval)
{
    double x, L;
    int i;

    L = lR - lL;


    for( i = 0; i < Nx; i++)
    {
        x = (double)i*hx + xL;
        if(x>=lL && x<=lL+tval && 0.<tval && tval<=2.*L)
            Uexact[i] 	= (x - lL) / tval;
        else if(x>=lL+tval && x<=lR+0.5*tval && 0.<tval && tval<=2.*L)
            Uexact[i] 	= 1.;
        else if(x>=lL && x<=lL+sqrt(2.*L*tval) && tval>2.*L)
            Uexact[i]		= (x - lL) / tval;
        else
            Uexact[i] 	= 0.;

    }
}

void lax_friedrichs()
{
    //Nx = pow(2, 6);
    hx = (xR - xL)/(Nx+1);
    dt = CFL*hx/fabs(a);


    double UarrN[Nx], t = 0.;
    int i, nt=0;

    for(i = 0; i < Nx; i++)
    {
        UarrN[i] = Uarr[i];
    }
    FILE *op;
    op = fopen("../dat/burgers/output_LF", "w");

    for(t = 0.; t<Tmax; t+=dt)
    {
        for(i = 1; i < Nx-1; i++)
            UarrN[i] = 0.5*((Uarr[i+1]+Uarr[i-1])-dt/hx*(f(Uarr[i+1])-f(Uarr[i-1])));

        // Ouput
        for(i = 0; i < Nx; i++)
        {
            Uarr[i] = UarrN[i];
            fprintf(op,"%f %f %f\n",(double)i*hx, t, Uarr[i]);

        }
        nt++;
    }

    fclose(op); // Close output in first dat file

    op = fopen("../dat/burgers/output_LF_last", "w");
    exact_solution(0.3,0.5,Tmax);
    for(i = 0; i < Nx; i++)
        fprintf(op,"%f %f %f\n", (double)i*hx, Uarr[i],
                Uexact[i]);
    fclose(op);

    printf("%d %f %f\n", nt, dt, Uarr[Nx-1]);

}
void lax_friedrichs_ch()
{
    //Nx = pow(2, 6);
    hx = (xR - xL)/(Nx+1);
    dt = CFL*hx/fabs(a);


    double UarrN[Nx], t = 0.;
    int i, nt=0;

    for(i = 0; i < Nx; i++)
    {
        UarrN[i] = Uarr[i];
    }
    FILE *op;
    op = fopen("../dat/burgers/output_LF", "w");

    for(t = 0.; t<Tmax; t+=dt)
    {
        for(i = 1; i < Nx-1; i+=2)
            UarrN[i+nt%2] = 0.5*((Uarr[i+nt%2+1]+Uarr[i+nt%2-1])-dt/hx*(f(Uarr[i+nt%2+1])-f(Uarr[i+nt%2-1])));

        // Ouput
        for(i = 0; i < Nx; i++)
        {
            Uarr[i] = UarrN[i];
            fprintf(op,"%f %f %f\n",(double)i*hx, t, Uarr[i]);

        }
        nt++;
    }

    fclose(op); // Close output in first dat file

    op = fopen("../dat/burgers/output_LF_last", "w");
    exact_solution(0.3,0.5,Tmax);
    for(i = 0; i < Nx; i++)
        fprintf(op,"%f %f %f\n", (double)i*hx, Uarr[i],
                Uexact[i]);
    fclose(op);

    printf("%d %f %f\n", nt, dt, Uarr[Nx-1]);

}

void init_cond()
{
    double x;
    int i;
    hx = (xR-xL)/(Nx-1);
    for( i = 0; i < Nx; i++)
    {
        x = (double)i*hx + xL;
        if(x>= 0.3 && x<=0.5)
            Uarr[i] = 1.;
        else
            Uarr[i] = 0.;
    }
}


int main()
{
    int ret_codes = 0;
    init_cond();
    lax_friedrichs_ch();
    ret_codes = system("../bin/plot_LF");
    ret_codes = system("../bin/plot_LF_last");


    return ret_codes;
}
