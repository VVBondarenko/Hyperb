#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double xL = 0., xR = 1., hx, t, Tmax=6, CFL = 0.5, dt, a = -2.;
#define Nx 128

double 	Uarr[Nx],
        dUarr[Nx];

double f_B_3(double X)
{
    double absV=fabs(X);
    if(absV<2)
    {
        if(absV>=1)
        {
            return 0.25*(2.0-absV)*(2.0-absV)*(2.0-absV);
        }
        else
            return 1.0 - 1.5*absV*absV*(1.0 - 0.5*absV);
    }
    return 0.0;
}

double f_d_B_3(double X)
{
    double absV=fabs(X);
    if(X>=0)
    {
        if(absV<2)
        {
            if(absV>=1)
                return -0.75*(2.0-absV)*(2.0-absV);
            else
                return -1.5*(2.0*absV - 1.5*absV*absV);
        }
    }
    else
    {
        if(absV<2)
        {
            if(absV>=1)
                return 0.75*(2.0-absV)*(2.0-absV);
            else
                return 1.5*(2.0*absV - 1.5*absV*absV);
        }
    }
    return 0.0;
}


void IDO()
{
    dt = CFL*hx/fabs(a);

    int i;
    double	UarrN[Nx],
//            dUarr[Nx],
            dUarrN[Nx];

    double A, B, C, D;

    for(i = 0; i < Nx; i++)
    {
        if(i == 0)
            dUarr[0] = (Uarr[1]-Uarr[0])/hx;
        else if(i==Nx-1)
            dUarr[Nx-1] = (Uarr[Nx-1]-Uarr[Nx-2])/hx;
        else
            dUarr[i] = (Uarr[i-1]+Uarr[i+1]-2.*Uarr[i])/hx;
    }


    FILE *op;
    int ts;

    for( t = 0.; t < Tmax; t+=dt )
    {
        i = 0;
        A = (-0.75/hx*(Uarr[i+1]-Uarr[Nx-1])+0.25*(dUarr[i+1]+4*dUarr[i]+dUarr[Nx-1]))/hx/hx/hx/hx;
        B = (-0.5/hx*(Uarr[i+1]-2*Uarr[i]+Uarr[Nx-1]) +0.25*(dUarr[i+1]-dUarr[Nx-1]))/hx/hx/hx;
        C = (1.25/hx*(Uarr[i+1]-Uarr[Nx-1])-0.25*(dUarr[i+1]+8*dUarr[i]+dUarr[Nx-1]))/hx/hx;
        D = (1./hx*(Uarr[i+1] -2.*Uarr[i]+Uarr[Nx-1]) -0.25*(dUarr[i+1]-dUarr[Nx-1]))/hx;
        dUarrN[i] = dUarr[i]-a*2.*dt*D + a*a*3.*dt*dt*C - a*a*a*4.*B*dt*dt*dt;
        UarrN[i]  = Uarr[i] -a*dt*dUarr[i] + a*a*dt*dt*D - a*a*a*dt*dt*dt*C
                        + a*a*a*a*dt*dt*dt*dt*B;

        for (i = 1; i < Nx-1; i++)
        {
            A = (-0.75/hx*(Uarr[i+1]-Uarr[i-1])+0.25*(dUarr[i+1]+4*dUarr[i]+dUarr[i-1]))/hx/hx/hx/hx;
            B = (-0.5/hx*(Uarr[i+1]-2*Uarr[i]+Uarr[i-1]) +0.25*(dUarr[i+1]-dUarr[i-1]))/hx/hx/hx;
            C = (1.25/hx*(Uarr[i+1]-Uarr[i-1])-0.25*(dUarr[i+1]+8*dUarr[i]+dUarr[i-1]))/hx/hx;
            D = (1./hx*(Uarr[i+1] -2.*Uarr[i]+Uarr[i-1]) -0.25*(dUarr[i+1]-dUarr[i-1]))/hx;

            //ToDo: add possibility to change a to values, that are different from 1
            dUarrN[i] = dUarr[i]-a*2.*dt*D + a*a*3.*dt*dt*C - a*a*a*4.*B*dt*dt*dt;
            UarrN[i]  = Uarr[i] -a*dt*dUarr[i] + a*a*dt*dt*D - a*a*a*dt*dt*dt*C
                            + a*a*a*a*dt*dt*dt*dt*B;
        }
        //так как рассматриваем случай a<0, то будем считать, что у нас всё движется налево,
        // а значит и нужна периодчиность только с лева. Просчитаем дополнительно окрестность
        // точки i=0

        i = Nx-1;
        A = (-0.75/hx*(Uarr[0]-Uarr[i-1])+0.25*(dUarr[0]+4*dUarr[i]+dUarr[i-1]))/hx/hx/hx/hx;
        B = (-0.5/hx*(Uarr[0]-2*Uarr[i]+Uarr[i-1]) +0.25*(dUarr[0]-dUarr[i-1]))/hx/hx/hx;
        C = (1.25/hx*(Uarr[0]-Uarr[i-1])-0.25*(dUarr[i+1]+8*dUarr[i]+dUarr[i-1]))/hx/hx;
        D = (1./hx*(Uarr[0] -2.*Uarr[i]+Uarr[i-1]) -0.25*(dUarr[0]-dUarr[i-1]))/hx;
        dUarrN[i] = dUarr[i]-a*2.*dt*D + a*a*3.*dt*dt*C - a*a*a*4.*B*dt*dt*dt;
        UarrN[i]  = Uarr[i] -a*dt*dUarr[i] + a*a*dt*dt*D - a*a*a*dt*dt*dt*C
                        + a*a*a*a*dt*dt*dt*dt*B;

        if(a<0)
        {
            UarrN[Nx-1] = UarrN[0];
            dUarrN[Nx-1] = dUarrN[0];
        }
        else
        {
            UarrN[0] = UarrN[Nx-1];
            dUarrN[0] = dUarrN[Nx-1];
        }
        for (i = 0; i < Nx; i++)
        {
            Uarr[i]  = UarrN[i];
            dUarr[i] = dUarrN[i];
        }
        ts++;
    }

    printf("%d\n",ts);

    op = fopen("output.IDO","w");
    for (i = 0; i < Nx; i++)
        fprintf (op, "%f %f\n", (double)i*hx + xL, Uarr[i]);
    fclose(op);
}

void CIP()
{
    dt = CFL*hx/fabs(a);

    int i;
    double	UarrN[Nx],
//            dUarr[Nx],
            dUarrN[Nx];

    double A, B, C, D, ksi;

    for(i = 0; i < Nx; i++)
    {
        if(i == 0)
            dUarr[0] = (Uarr[1]-Uarr[0])/hx;
        else if(i==Nx-1)
            dUarr[Nx-1] = (Uarr[Nx-1]-Uarr[Nx-2])/hx;
        else
            dUarr[i] = (Uarr[i-1]+Uarr[i+1]-2.*Uarr[i])/hx;
    }


    FILE *op;
    int ts;
    if(a > 0.)
    {
        for( t = 0.; t < Tmax; t+=dt )
        {
            for (i = 1; i < Nx; i++)
            {
                A = 2.*(Uarr[i-1]-Uarr[i])/hx/hx/hx + (dUarr[i-1]+dUarr[i])/hx/hx;
                B = 3.*(Uarr[i-1]-Uarr[i])/hx/hx + (2.*dUarr[i]+dUarr[i-1])/hx;
                C = dUarr[i];
                D = Uarr[i];

                ksi =  -a*dt;

                UarrN[i]  = ((A*ksi + B)*ksi+C)*ksi+D;
                dUarrN[i] = ((3.*A*ksi + 2.*B)*ksi+C);
            }
            UarrN[0] = UarrN[Nx-1];
            dUarrN[0] = dUarrN[Nx-1];
            for (i = 0; i < Nx; i++)
            {
                Uarr[i]  = UarrN[i];
                dUarr[i] = dUarrN[i];
            }
            ts++;
        }
    }


    if(a < 0.)
    {
        for(t = 0.; t < Tmax; t+=dt)
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
            }

            UarrN[Nx-1] = UarrN[0];
            dUarrN[Nx-1] = dUarrN[0];
            for (i = 0; i < Nx; i++)
            {
                Uarr[i]  = UarrN[i];
                dUarr[i] = dUarrN[i];
            }
            ts++;
        }
    }
    printf("%d\n",ts);

    op = fopen("output.CIP","w");
    for (i = 0; i < Nx; i++)
        fprintf (op, "%f %f\n", (double)i*hx + xL, Uarr[i]);
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
    if(a<0.)
    {
        for(t = 0.; t < Tmax; t += dt)
        {
            for(i = 0; i < Nx-1; i++)
                UarrN[i] = Uarr[i] + CFL*(Uarr[i+1] - Uarr[i]);

            UarrN[Nx-1] = UarrN[0];

            for(i = 0; i < Nx; i++)
                Uarr[i] = UarrN[i];
            nt++;
        }
    }
    if(a > 0.)
    {
        for(t = 0.; t<Tmax; t+=dt)
        {
            for(i = 1; i < Nx; i++)
                UarrN[i] = Uarr[i] - CFL*(Uarr[i] - Uarr[i-1]);

            UarrN[0] = UarrN[Nx-1];

            for(i = 0; i < Nx; i++)
                Uarr[i] = UarrN[i];
            nt++;
        }
    }

    for(i = 0; i < Nx; i++)
        fprintf(op,"%f %f\n",(double)i*hx,Uarr[i]);
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
        if(x> 0.7 && x<0.8)
        {
//            Uarr[i] = f_B_3((x-0.75)/0.1*2.);
            Uarr[i] = 1.;
//            dUarr[i]= f_d_B_3((x-0.75)/0.1*2.);
        }
        else
        {
            Uarr[i] = 0.;
//            dUarr[i] = 0.;
        }
    }
}

int main()
{
    init_cond();
    IDO();
    init_cond();
    CIP();
    init_cond();
//    CIP();
    upwind();

    return 0;
}
