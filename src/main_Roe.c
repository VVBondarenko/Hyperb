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
        fprintf(op,"%f %f %f\n", (double)i*hx, Uarr[i],
                Uexact[i]);
    fclose(op);
    printf("%d\n", nt);

}


int main()
{
	int ret_codes = 0;
	init_task(1);
    init_cond(0.3,0.5);
    roe();
	ret_codes = system("../bin/plot_Roe");
	ret_codes = system("../bin/plot_Roe_last");
	
	
    return ret_codes;
}
