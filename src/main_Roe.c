#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double xL = 0., xR = 1., hx, t, Tmax=0.5, CFL = 0.9, dt, a = 1.;
#define Nx 127

double 	Uarr[Nx];
double 	Uexact[Nx];

#include "tasks.c"

double f(double u)
{
    return 0.5*u*u;
}

double df_du(double u)
{
    return u;
}

double fminus(double u)
{
    return fmin(df_du(u),0.);
}

double fplus(double u)
{
    return fmax(df_du(u),0.);
}

double integral(double (*f)(double), double x0, double x1)
{
    double	 step = (x1-x0)/128.;

    const int   q_of_layers = (int)((x1-x0)/step) + 1;
    int         k;

    double res = (*f)(x0) +
                 (*f)((double)(q_of_layers-1)*step + x0);

    for(k = 1; k < q_of_layers-1; k++)
    {
        res += 2.*(double)(1+k%2)*(*f)((double)k*step + x0);
    }
    return step/3.*res;
}

void roe()
{
    hx = (xR - xL)/(Nx+1);
    dt = CFL*hx/fabs(a);

    double UarrN[Nx], t = 0.;
    int i, nt=0;
    double Fplus, Fminus;

    /* Set init */
    for(i = 0; i < Nx; i++)
    {
        UarrN[i] = Uarr[i];
    }

    FILE *op;
    op = fopen("../dat/burgers/output_Roe", "w");

    for(t = 0.; t < Tmax; t += dt)
    {
        for(i = 2; i < Nx-2; i++)
        {
            /* Roe average speed*/
            a = df_du((Uarr[i]+Uarr[i-1])*0.5);
            /* Fluxes */
            if(a<0.)
            {
                Fplus = f(Uarr[i+1]);
                Fminus = f(Uarr[i]);
            }
            else
            {
                Fplus = f(Uarr[i]);
                Fminus = f(Uarr[i-1]);
            }
            UarrN[i] = Uarr[i] - dt/hx*(Fplus-Fminus);
        }

        /* Update */
        for(i = 0; i < Nx; i++)
        {
            Uarr[i] = UarrN[i];
            fprintf(op,"%f %f %f\n",(double)i*hx, t, Uarr[i]);
        }
        nt++;
    }
    fclose(op);

    op = fopen("../dat/burgers/output_Roe_last", "w");
    exact_solution(0.3,0.5,Tmax);
    for(i = 0; i < Nx; i++)
        fprintf(op,"%f %f %f\n", (double)i*hx, Uarr[i], Uexact[i]);
    fclose(op);
    printf("%d\n", nt);

}
void rusanov()
{
    hx = (xR - xL)/(Nx+1);
    dt = CFL*hx/fabs(a);

    double UarrN[Nx], t = 0.;
    int i, nt=0;
    double Fplus, Fminus;

    /* Set init */
    for(i = 0; i < Nx; i++)
    {
        UarrN[i] = Uarr[i];
    }

    FILE *op;
    op = fopen("../dat/burgers/output_Rusanov", "w");

    for(t = 0.; t < Tmax; t += dt)
    {
        for(i = 2; i < Nx-2; i++)
        {
            Fplus 		= (f(Uarr[i])+f(Uarr[i+1]))*0.5
                          - fmax(fabs(df_du(Uarr[i])),fabs(df_du(Uarr[i+1])))*(Uarr[i+1]-Uarr[i])*0.5;
            Fminus 		= (f(Uarr[i-1])+f(Uarr[i]))*0.5
                          - fmax(fabs(df_du(Uarr[i-1])),fabs(df_du(Uarr[i])))*(Uarr[i]-Uarr[i-1])*0.5;
            UarrN[i] = Uarr[i] - dt/hx*(Fplus-Fminus);
        }

        /* Update */
        for(i = 0; i < Nx; i++)
        {
            Uarr[i] = UarrN[i];
            fprintf(op,"%f %f %f\n",(double)i*hx, t, Uarr[i]);
        }
        nt++;
    }
    fclose(op);

    op = fopen("../dat/burgers/output_Rusanov_last", "w");
    exact_solution(0.3,0.5,Tmax);
    for(i = 0; i < Nx; i++)
        fprintf(op,"%f %f %f\n", (double)i*hx, Uarr[i], Uexact[i]);
    fclose(op);
    printf("%d\n", nt);

}
void llf()
{
    hx = (xR - xL)/(Nx+1);
    dt = CFL*hx/fabs(a);

    double UarrN[Nx], t = 0.;
    int i, nt=0;
    double Fplus, Fminus,alpha;

    /* Set init */
    for(i = 0; i < Nx; i++)
    {
        UarrN[i] = Uarr[i];
    }

    FILE *op;
    op = fopen("../dat/burgers/output_llf", "w");

    for(t = 0.; t < Tmax; t += dt)
    {
        for(i = 2; i < Nx-2; i++)
        {
            alpha		= fmax(fabs(df_du(Uarr[i])),fabs(df_du(Uarr[i+1])));
            Fplus 		= 0.5*(f(Uarr[i])+f(Uarr[i+1])-alpha*(Uarr[i+1]-Uarr[i]));

            alpha		= fmax(fabs(df_du(Uarr[i-1])),fabs(df_du(Uarr[i])));
            Fminus 		= 0.5*(f(Uarr[i-1])+f(Uarr[i])-alpha*(Uarr[i]-Uarr[i-1]));
            UarrN[i] = Uarr[i] - dt/hx*(Fplus-Fminus);
        }

        /* Update */
        for(i = 0; i < Nx; i++)
        {
            Uarr[i] = UarrN[i];

            fprintf(op,"%f %f %f\n",(double)i*hx, t, Uarr[i]);
        }
        nt++;
    }
    fclose(op);

    op = fopen("../dat/burgers/output_llf_last", "w");
    exact_solution(0.3,0.5,Tmax);
    for(i = 0; i < Nx; i++)
        fprintf(op,"%f %f %f\n", (double)i*hx, Uarr[i], Uexact[i]);
    fclose(op);
    printf("%d\n", nt);

}
void EO()
{
    hx = (xR - xL)/(Nx+1);
    dt = CFL*hx/fabs(a);

    double UarrN[Nx], t = 0.;
    int i, nt=0;
    double Fplus, Fminus;

    /* Set init */
    for(i = 0; i < Nx; i++)
    {
        UarrN[i] = Uarr[i];
    }

    FILE *op;
    op = fopen("../dat/burgers/output_EO", "w");

    for(t = 0.; t < Tmax; t += dt)
    {
        for(i = 2; i < Nx-2; i++)
        {
            Fplus 		= integral(fplus,0.,Uarr[i]) + integral(fminus,0.,Uarr[i+1]);
            Fminus 		= integral(fplus,0.,Uarr[i-1]) + integral(fminus,0.,Uarr[i]);
            UarrN[i] = Uarr[i] - dt/hx*(Fplus-Fminus);
        }

        /* Update */
        for(i = 0; i < Nx; i++)
        {
            Uarr[i] = UarrN[i];

            fprintf(op,"%f %f %f\n",(double)i*hx, t, Uarr[i]);
        }
        nt++;
    }
    fclose(op);

    op = fopen("../dat/burgers/output_EO_last", "w");
    exact_solution(0.3,0.5,Tmax);
    for(i = 0; i < Nx; i++)
        fprintf(op,"%f %f %f\n", (double)i*hx, Uarr[i], Uexact[i]);
    fclose(op);
    printf("%d\n", nt);

}


int main()
{
    int ret_codes = 0;
    init_task(4);
    init_cond(0.3,0.5);
    llf();
    ret_codes = system("../bin/plot_llf");
    ret_codes = system("../bin/plot_llf_last");

    return ret_codes;
}
