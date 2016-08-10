#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "smplDynArray.c"


double xL = 0., xR = 1., hx, t, Tmax=0.001, CFL = 0.9, dt, a = 1.;
int Nx = 512;
void (*solve)();
//double 	Uarr[Nx];
//double 	dUarr[Nx];
//double 	Uexact[Nx];
double* 	Uarr;
double* 	dUarr;
double* 	Uexact;

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

double psi(double r)
{
    return 0.;
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

double errC()
{
    int i;
    double temp = 0.;
    for (i = 2; i < Nx-2; i++)
    {
        temp = fmax(temp, fabs(Uarr[i]-Uexact[i]));
    }
    return temp;
}
double errL1()
{
    int i;
    double temp = 0.;
    for (i = 2; i < Nx-2; i++)
    {
        temp += fabs(Uarr[i]-Uexact[i]);
    }
    return temp*hx;
}
void metr_uref_l1()
{
    int i;
    double temp = 0.;
    for (i = 0; i < Nx; i++)
    {
        temp += fabs(Uexact[i]);
    }
    FILE* urefs_out;
    urefs_out = fopen("../uref.dat", "w");
    fprintf(urefs_out,"%3.3f",temp);
    fclose(urefs_out);
}
double errL2()
{
    int i;
    double temp = 0.;
    for (i = 2; i < Nx-2; i++)
    {
        temp += fabs(Uarr[i]-Uexact[i])*fabs(Uarr[i]-Uexact[i]);
    }
    return sqrt(temp*hx);
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
            fprintf(op,"%3.3f %3.3f %3.3f\n",(double)i*hx, t, Uarr[i]);
        }
        nt++;
    }
    fclose(op);

    op = fopen("../dat/burgers/output_Roe_last", "w");
    exact_solution(0.3,0.5,Tmax);
    for(i = 0; i < Nx; i++)
        fprintf(op,"%3.3f %3.3f %3.3f %3.3f\n", (double)i*hx, Uarr[i], Uexact[i], fabs(Uarr[i]-Uexact[i]));
    fclose(op);
    printf("%d %3.3f %3.3f %3.3f\n", Nx, errC(), errL1(),errL2());
	//printf("%d", nt);
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
            fprintf(op,"%3.3f %3.3f %3.3f\n",(double)i*hx, t, Uarr[i]);
        }
        nt++;
    }
    fclose(op);

    op = fopen("../dat/burgers/output_Rusanov_last", "w");
    exact_solution(0.3,0.5,Tmax);
    for(i = 0; i < Nx; i++)
        fprintf(op,"%3.3f %3.3f %3.3f %3.3f\n", (double)i*hx, Uarr[i], Uexact[i], fabs(Uarr[i]-Uexact[i]));
    fclose(op);
    printf("%d %3.3f %3.3f %3.3f\n", Nx, errC(), errL1(),errL2());

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

            fprintf(op,"%3.3f %3.3f %3.3f\n",(double)i*hx, t, Uarr[i]);
        }
        nt++;
    }
    fclose(op);

    op = fopen("../dat/burgers/output_llf_last", "w");
    exact_solution(0.3,0.5,Tmax);
    for(i = 0; i < Nx; i++)
        fprintf(op,"%3.3f %3.3f %3.3f %3.3f\n", (double)i*hx, Uarr[i], Uexact[i], fabs(Uarr[i]-Uexact[i]));
    fclose(op);
    printf("%d %3.3f %3.3f %3.3f\n", Nx, errC(), errL1(),errL2());

}
void llf_ch()
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
        for(i = 2; i < Nx-2; i+=2)
        {
            alpha		= fmax(fabs(df_du(Uarr[i+nt%2])),fabs(df_du(Uarr[i+nt%2+1])));
            Fplus 		= 0.5*(f(Uarr[i+nt%2])+f(Uarr[i+nt%2+1])-alpha*(Uarr[i+nt%2+1]-Uarr[i+nt%2]));

            alpha		= fmax(fabs(df_du(Uarr[i+nt%2-1])),fabs(df_du(Uarr[i+nt%2])));
            Fminus 		= 0.5*(f(Uarr[i+nt%2-1])+f(Uarr[i+nt%2])-alpha*(Uarr[i+nt%2]-Uarr[i+nt%2-1]));
            UarrN[i+nt%2] = Uarr[i+nt%2] - dt/hx*(Fplus-Fminus);
        }

        /* Update */
        for(i = 0; i < Nx; i++)
        {
            Uarr[i] = UarrN[i];

            fprintf(op,"%3.3f %3.3f %3.3f\n",(double)i*hx, t, Uarr[i]);
        }
        nt++;
    }
    fclose(op);

    op = fopen("../dat/burgers/output_llf_last", "w");
    exact_solution(0.3,0.5,Tmax);
    for(i = 0; i < Nx; i++)
        fprintf(op,"%3.3f %3.3f %3.3f %3.3f\n", (double)i*hx, Uarr[i], Uexact[i], fabs(Uarr[i]-Uexact[i]));
    fclose(op);
    printf("%d %3.3f %3.3f %3.3f\n", Nx, errC(), errL1(),errL2());

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

            fprintf(op,"%3.3f %3.3f %3.3f\n",(double)i*hx, t, Uarr[i]);
        }
        nt++;
    }
    fclose(op);

    op = fopen("../dat/burgers/output_EO_last", "w");
    exact_solution(0.3,0.5,Tmax);
    for(i = 0; i < Nx; i++)
        fprintf(op,"%3.3f %3.3f %3.3f %3.3f\n", (double)i*hx, Uarr[i], Uexact[i], fabs(Uarr[i]-Uexact[i]));
    fclose(op);
    printf("%d %3.3f %3.3f %3.3f\n", Nx, errC(), errL1(),errL2());

}
void tvd()
{
    hx = (xR - xL)/(Nx+1);
    dt = CFL*hx/fabs(a);

    double UarrN[Nx], t = 0.;
    int i, nt=0;
    double Fplus, Fminus, r;

    /* Set init */
    for(i = 0; i < Nx; i++)
    {
        UarrN[i] = Uarr[i];
    }

    FILE *op;
    op = fopen("../dat/burgers/output_tvd", "w");

    for(t = 0.; t < Tmax; t += dt)
    {
        for(i = 2; i < Nx-2; i++)
        {
            r			= (Uarr[i+1]-Uarr[i])/(Uarr[i]-Uarr[i-1]);
            Fplus 		= Uarr[i] + 0.5*psi(r)*(Uarr[i]-Uarr[i-1]);
            Fminus 		= Uarr[i-1] + 0.5*psi(r)*(Uarr[i-1]-Uarr[i-2]);
            UarrN[i] = Uarr[i] - CFL*(Fplus-Fminus);
        }

        /* Update */
        for(i = 0; i < Nx; i++)
        {
            Uarr[i] = UarrN[i];

            fprintf(op,"%3.3f %3.3f %3.3f\n",(double)i*hx, t, Uarr[i]);
        }
        nt++;
    }
    fclose(op);

    op = fopen("../dat/burgers/output_tvd_last", "w");
    exact_solution(0.3,0.5,Tmax);
    for(i = 0; i < Nx; i++)
        fprintf(op,"%3.3f %3.3f %3.3f %3.3f\n", (double)i*hx, Uarr[i], Uexact[i], fabs(Uarr[i]-Uexact[i]));
    fclose(op);
    printf("%d %3.3f %3.3f %3.3f\n", Nx, errC(), errL1(),errL2());

}

double df_dt(int i)
{
    return -df_du(Uarr[i])*df_du(Uarr[i])*dUarr[i];
}
void bgp()
{
    hx = (xR - xL)/(Nx+1);
    dt = CFL*hx/fabs(a);

    double UarrN[Nx],dUarrN[Nx], t = 0.;
    int i, nt=0;
    //double Fplus, Fminus;
    dUarr[0] = 0.;
    dUarr[Nx-1] = 0.;
    for (i = 1; i < Nx-1; i++)
    {
        dUarr[i] = (Uarr[i+1]-Uarr[i-1])/hx*0.5;
    }


    /* Set init */
    for(i = 0; i < Nx; i++)
    {
        UarrN[i]	= Uarr[i];
        dUarrN[i]	= dUarr[i];
    }

    FILE *op;
    op = fopen("../dat/burgers/output_bgp", "w");

    for(t = 0.; t < Tmax; t += dt)
    {
        for(i = 2; i < Nx-2; i++)// nt%2
        {
            dUarrN[i] 	= 0.5*(dUarr[i-1]+dUarr[i+1]) - 0.5*dt/hx*(df_du(Uarr[i+1])*dUarr[i+1]-df_du(Uarr[i-1])*dUarr[i-1]);
            UarrN[i]	= 0.5*(Uarr[i-1]+Uarr[i+1]) - 0.5*dt/hx*(f(Uarr[i+1])-f(Uarr[i-1])) - 0.25*hx*(dUarr[i+1]-dUarr[i-1])
                          - dt*dt/hx*0.25*(df_dt(i+1)-df_dt(i-1));
        }

        /* Update */
        for(i = 0; i < Nx; i++)
        {
            Uarr[i] = UarrN[i];
            dUarr[i]= dUarrN[i];
            fprintf(op,"%3.3f %3.3f %3.3f\n",(double)i*hx, t, Uarr[i]);
        }
        nt++;
    }
    fclose(op);

    op = fopen("../dat/burgers/output_bgp_last", "w");
    exact_solution(0.3,0.5,Tmax);
    for(i = 0; i < Nx; i++)
        fprintf(op,"%3.3f %3.3f %3.3f %3.3f\n", (double)i*hx, Uarr[i], Uexact[i], fabs(Uarr[i]-Uexact[i]));
    fclose(op);
    printf("%d %3.3f %3.3f %3.3f\n", Nx, errC(), errL1(),errL2());

}
void bgp_ch()
{
    hx = (xR - xL)/(Nx+1);
    dt = CFL*hx/fabs(a);

    double UarrN[Nx],dUarrN[Nx], t = 0.;
    int i, nt=0;
    //double Fplus, Fminus;
    dUarr[0] = 0.;
    dUarr[Nx-1] = 0.;
    for (i = 1; i < Nx-1; i++)
    {
        dUarr[i] = (Uarr[i+1]-Uarr[i-1])/hx*0.5;
    }


    /* Set init */
    for(i = 0; i < Nx; i++)
    {
        UarrN[i]	= Uarr[i];
        dUarrN[i]	= dUarr[i];
    }

    FILE *op;
    op = fopen("../dat/burgers/output_bgp", "w");

    for(t = 0.; t < Tmax; t += dt)
    {
        for(i = 2; i < Nx-2; i+=2)// nt%2
        {
            dUarrN[i+nt%2] 	= 0.5*(dUarr[i+nt%2-1]+dUarr[i+nt%2+1]) - 0.5*dt/hx*(df_du(Uarr[i+nt%2+1])*dUarr[i+nt%2+1]-df_du(Uarr[i+nt%2-1])*dUarr[i+nt%2-1]);
            UarrN[i+nt%2]	= 0.5*(Uarr[i+nt%2-1]+Uarr[i+nt%2+1]) - 0.5*dt/hx*(f(Uarr[i+nt%2+1])-f(Uarr[i+nt%2-1])) - 0.25*hx*(dUarr[i+nt%2+1]-dUarr[i+nt%2-1])
                              - dt*dt/hx*0.25*(df_dt(i+nt%2+1)-df_dt(i+nt%2-1));
        }

        /* Update */
        for(i = 0; i < Nx; i++)
        {
            Uarr[i] = UarrN[i];
            dUarr[i]= dUarrN[i];
            fprintf(op,"%3.3f %3.3f %3.3f\n",(double)i*hx, t, Uarr[i]);
        }
        nt++;
    }
    fclose(op);

    op = fopen("../dat/burgers/output_bgp_last", "w");
    exact_solution(0.3,0.5,Tmax);
    for(i = 0; i < Nx; i++)
        fprintf(op,"%3.3f %3.3f %3.3f %3.3f\n", (double)i*hx, Uarr[i], Uexact[i], fabs(Uarr[i]-Uexact[i]));
    fclose(op);
    printf("%d %3.3f %3.3f %3.3f\n", Nx, errC(), errL1(),errL2());

}



int main(int argc, char** argv)
{
    int ret_codes = 0, methodN = 1;
    Nx = atoi(argv[2]);
    CFL = atof(argv[4]);
	methodN = atoi(argv[3]);
	if(methodN == 1)
		solve = &roe;
	if(methodN == 2)
		solve = &rusanov;
	if(methodN == 3)
		solve = &llf;
	if(methodN == 4)
		solve = &EO;
	if(methodN == 5)
		solve = &tvd;
	if(methodN == 6)
		solve = &bgp;
	if(methodN == 7)
		solve = &bgp_ch;
			
    /*	Initing arrays	*/
    init1DArr(&Uarr,Nx);
    init1DArr(&dUarr,Nx);
    init1DArr(&Uexact,Nx);
	
    init_task(atoi(argv[1]));
    init_cond(0.3,0.5);
	solve();
	metr_uref_l1();
    //ret_codes = system("../bin/plotter.py ../dat/burgers/output_llf&");
    //ret_codes = system("../bin/plot_bgp");
    //ret_codes = system("../bin/plot_Roe_last");

    /*	Cleaning arrays	*/
    free1DArr(&Uarr);
    free1DArr(&dUarr);
    free1DArr(&Uexact);

    return ret_codes;
}
