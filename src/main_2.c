#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double xL = 0., xR = 1., hx, t, Tmax=0.034, CFL = 0.9, dt, a = -1.;
#define Nx 64

double 	Uarr[Nx];

void CIP()
{

    
    dt = CFL*hx/fabs(a);


    int i;
    double	UarrN[Nx],
            dUarr[Nx], dUarrN[Nx];

    double A, B, C, D, ksi;



    for(i = 0; i < Nx; i++)
    {
        if(i == 0)
        {
            dUarr[0] = (Uarr[1]-Uarr[0])/hx;
        }
        else if(i==Nx-1)
        {
            dUarr[Nx-1] = (Uarr[Nx-1]-Uarr[Nx-2])/hx;
        }
        else
        {
            dUarr[i] = (Uarr[i-1]+Uarr[i+1]-2.*Uarr[i])/hx;
        }
    }


    printf("%f\n", dt);
    FILE *op;

    int ts;
    for( t = 0.; t < Tmax; t+=dt )
    {
        for (i = 0; i < Nx-1; i++)
        {
            A = -2.*(Uarr[i+1]-Uarr[i])/hx/hx/hx + (dUarr[i+1]+dUarr[i])/hx/hx;
            B = 3.*(Uarr[i+1]-Uarr[i])/hx/hx - (2.*dUarr[i]+dUarr[i+1])/hx;
            C = dUarr[i];
            D = Uarr[i];

            ksi =  -a*dt;

            UarrN[i]  = ((A*ksi + B)*ksi+C)*ksi+D;
            dUarrN[i] = ((3.*A*ksi + 2.*B)*ksi+C);

            //UarrN[i] = Uarr[i] - CFL*(Uarr[i] - Uarr[i-1]);
        }

        UarrN[Nx-1] = UarrN[0];
        dUarrN[Nx-1] = dUarrN[0];
        for (i = 0; i < Nx; i++)
        {
            Uarr[i]  = UarrN[i];
            dUarr[i] = dUarrN[i];
            //fprintf (op, "%f %f %f\n", (double)i*hx + xL, t, Uarr[i]);
        }
        ts++;
    }
    //fclose(op);
    printf("%d\n",ts);

    op = fopen("output.CIP","w");
    for (i = 0; i < Nx; i++)
    {
        fprintf (op, "%f %f\n", (double)i*hx + xL, Uarr[i]);
    }
    fclose(op);
}

void upwind()
{
    //Nx = pow(2, 6);
    //hx = (xR - xL)/(Nx+1);
    dt = CFL*hx/fabs(a);

    double UarrN[Nx], t = 0.;
    int i, nt=0;


    FILE *op;
    op = fopen("output.upwind", "w");
    for(t=0.; t<Tmax; t+=dt)
    {
        for(i = 0; i < Nx-1; i++)
        {
            UarrN[i] = Uarr[i] - CFL*(Uarr[i+1] - Uarr[i]);
        }

        UarrN[Nx-1] = UarrN[0];

        for(i = 0; i < Nx; i++)
        {
            Uarr[i] = UarrN[i];
        }

        nt++;
    }
    for(i = 0; i < Nx; i++)
    {
        fprintf(op,"%f %f\n",(double)i*hx,Uarr[i]);
    }
    fclose(op);
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
        if(x> 0.1 && x<0.3)
            Uarr[i] = 1.;
        else
            Uarr[i] = 0.;
    }	
}

int main()
{

    init_cond();
    CIP();
    init_cond();
    upwind();
    
	return 0;
}
