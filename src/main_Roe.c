#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double xL = 0., xR = 1., hx, t, Tmax=0.3, CFL = 0.9, dt, a = 1.;
#define Nx 128

double 	Uarr[Nx];
double 	Uexact[Nx];

double f(double u)
{
	return 0.5*u*u;
}

double df_du(double u)
{
	return u;
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


void roe()
{
    hx = (xR - xL)/(Nx+1);
    dt = CFL*hx/fabs(a);

    double UarrN[Nx], t = 0.;
    int i, nt=0;


    FILE *op;
    op = fopen("../dat/burgers/output_Roe", "w");

    for(t = 0.; t < Tmax; t += dt)
    {
        for(i = 1; i < Nx-1; i++)
        {
			a = df_du((Uarr[i]+Uarr[i-1])*0.5);
            if(a<0.)
            {
				a = df_du((Uarr[i+1]+Uarr[i])*0.5);
                UarrN[i] = Uarr[i] + a*dt/hx*(Uarr[i+1] - Uarr[i]);
			}
            if(a>0.)
                UarrN[i] = Uarr[i] - a*dt/hx*(Uarr[i] - Uarr[i-1]);
        }
		UarrN[Nx-1] 	= 0.;
        for(i = 0; i < Nx; i++)
        {
			if(UarrN[i]==UarrN[i])
				Uarr[i] = UarrN[i];
			else 
				Uarr[i] = 0.;
            
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
    init_cond();
    roe();
	ret_codes = system("../bin/plot_Roe");
	ret_codes = system("../bin/plot_Roe_last");
	
	
    return ret_codes;
}
