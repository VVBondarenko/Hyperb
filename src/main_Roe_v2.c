#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double xL = 0., xR = 1., hx, t, Tmax=0.1, CFL = 0.9, dt, a = 1.;
#define Nx 256

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

    double UarrN[Nx], t = 0., fp, fm;
    int i, nt=0;


    FILE *op;
    op = fopen("../dat/burgers/output_Roe", "w");

    for(t = 0.; t < Tmax; t += dt)
    {
        for(i = 2; i < Nx-2; i++)
        {
			a = df_du((Uarr[i]+Uarr[i-1])*0.5);
            if(a<0.)
            {
				//a = df_du((Uarr[i+1]+Uarr[i])*0.5);
				fp = f(Uarr[i+1]);
				fm = f(Uarr[i]);
                //UarrN[i] = Uarr[i] - dt/hx*(fp - fm);
			}
            else
            {
				fp = f(Uarr[i]);
				fm = f(Uarr[i-1]);                
			}
			UarrN[i] = Uarr[i] - dt/hx*(fp - fm);
        }
		//UarrN[Nx-1] 	= 0.;
		//UarrN[Nx-2] 	= 0.;
		//UarrN[Nx-3] 	= 0.;
		//UarrN[Nx-4] 	= 0.;
		
		/* BC */
		UarrN[0] = Uexact[0];
		
		/* Update */
        for(i = 0; i < Nx; i++)
        {
			/*if(UarrN[i]==UarrN[i])
				Uarr[i] = UarrN[i];
			else 
				Uarr[i] = 0.;*/
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
	init_task(3);
    init_cond(0.3,0.5);
    roe();
	ret_codes = system("../bin/plot_Roe");
	ret_codes = system("../bin/plot_Roe_last");
	
	
    return ret_codes;
}
