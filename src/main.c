#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double xL = 0., xR = 1., hx, t, Tmax=0.5, CFL = 0.9, dt, a = 1.;
int Nx = 64;

int main()
{
	hx = (xR-xL)/(Nx-1);
	dt = CFL*hx/a;
	
	
	int i;
	double Uarr[Nx], UarrN[Nx], x;
	
	for( i = 0; i < Nx; i++)
	{
		x = (double)i*hx + xL;
		if(x> 0.1 && x<0.3)
			Uarr[i] = 1.;
		else
			Uarr[i] = 0.;
	}
	printf("%f\n", Tmax);
	FILE *op;
	op = fopen("output","w");
	
	for( t = 0.; t < Tmax; t+=dt )
	{
		for (i = 1; i < Nx; i++)
		{
			UarrN[i] = Uarr[i] - CFL*(Uarr[i] - Uarr[i-1]); 
		}
		
		for (i = 0; i < Nx; i++)
		{
			Uarr[i] = UarrN[i];
			fprintf (op, "%f %f %f\n", (double)i*hx + xL, t, Uarr[i]);
		}
	}
	fclose(op);
	
	
	return 0;
}
