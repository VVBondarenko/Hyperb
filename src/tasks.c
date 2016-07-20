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
