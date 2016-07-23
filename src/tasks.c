void (*exact_solution)(double, double, double);
void (*init_cond)(double, double);

void ic1(double lL, double lR)
{
    double x;
    int i;
    hx = (xR-xL)/(Nx-1);
    for( i = 0; i < Nx; i++)
    {
        x = (double)i*hx + xL;
        if(x>= lL && x<=lR)
            Uarr[i] = 1.;
        else
            Uarr[i] = 0.;
    }
}
void es1(double lL, double lR, double tval)
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

void ic2(double lL, double lR)
{
    double x;
    int i;
    hx = (xR-xL)/(Nx-1);
    for( i = 0; i < Nx; i++)
    {
        x = (double)i*hx + xL;
        if (x>= lL && x<(lR+lL)*0.5)
        {
            Uarr[i] = 2.*(x-lL)/(lR-lL);
        }
        else if (x>=(lR+lL)*0.5 && x<=lR)
        {
            Uarr[i] = 2.*(lR-x)/(lR-lL);
        }
        else
        {
            Uarr[i] = 0.;
        }
    }
}
void es2(double lL, double lR, double tval)
{
    double x, L=lR-lL;
    int i;
    for( i = 0; i < Nx; i++)
    {
        x = (double)i*hx + xL;
        if ((x>=lL && x<=lL+sqrt(L*(L+2.*tval)*0.5) && tval >L*0.5) ||
                (x>=lL && x<=(lL+lR)*0.5+tval && tval>=0 && tval<=L*0.5))
        {
            Uexact[i] = 2.*(x-lL)/(L+2.*tval);
        }
        else if (x>=(lL+lR)*0.5+tval && x<=lR && tval>0 && tval<=L*0.5)
        {
            Uexact[i] = 2.*(lR-x)/(lR-lL-2.*tval);
        }
        else if ((!(x>=lL && x<=lL+sqrt(L*(L+2.*tval)*0.5))&& tval>L*0.5) ||
                 (!(x>=lL && x<=lR) && tval>0 && tval<=L*0.5))
        {
            Uexact[i] = 0.;
        }


    }
}

void ic3(double lL, double lR)
{
    double x;
    int i;
    hx = (xR-xL)/(Nx-1);
    for( i = 0; i < Nx; i++)
    {
        x = (double)i*hx + xL;
        if(x<= lL)
            Uarr[i] = 1.;
        else
            Uarr[i] = 0.;
    }
}
void es3(double lL, double lR, double tval)
{
	double x;
    int i;

    for( i = 0; i < Nx; i++)
    {
        x = (double)i*hx + xL;
		if (x<=lL+tval*0.5)
		{
			Uexact[i] = 1.;
		}
		else
		{
			Uexact[i] = 0.;
		}
    }
}

void ic4(double lL, double lR)
{
    double x;
    int i;
    hx = (xR-xL)/(Nx-1);
    for( i = 0; i < Nx; i++)
    {
        x = (double)i*hx + xL;
        if(x<= lL)
            Uarr[i] = 0.;
        else
            Uarr[i] = 1.;
    }
}
void es4(double lL, double lR, double tval)
{
	double x;
    int i;

    for( i = 0; i < Nx; i++)
    {
        x = (double)i*hx + xL;
		if (x<=lL)
		{
			Uexact[i] = 0.;
		}
		else if (x>lL && x<=(lL+tval))
		{
			Uexact[i] = (x-lL)/tval;
		}
		else
		{
			Uexact[i] = 1.;
		}
		
    }
}


void init_task(int id)
{
    exact_solution	= &es1;			//rectangle
    init_cond		= &ic1;
    if (id == 2)
    {
        exact_solution	= &es2;		//triangle
        init_cond		= &ic2;
    }
    if (id == 3)
	{
		exact_solution	= &es3;
		init_cond		= &ic3;
	}
	if (id == 4)
	{
		exact_solution	= &es4;
		init_cond		= &ic4;
	}
	
	
	
}
