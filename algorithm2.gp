/*
    Effective Bounds for Singular Units     
    Yuri Bilu, Philipp Habegger, Lars Kuehne 

    PARI code to verify computations in Proposition 8.1 Implements
    Algorithm 2
    
    Rigorously verifies correct rounding when using floating point
    arithmetic.

    If successful, the script prints a (potentially) proper supset
    of those discriminants for which the associated singular
    modulus is potentially a singular unit.

    It terminates and prints "FAIL" if the used floating point
    arithmetic does not round as expected. Although this is unlikely to
    happen, we need this as a safeguard against floating point
    errors.

    Tested using PARI 2.9.3.
*/

get_log_norm(Delta) = 
{
    local(a,b,c,clno,X,m,maxa,sqrtXdiva,jnorm,y);
    
    clno = 0;

    X = abs(Delta);

    m = 42700*min(2/(5*X),1/250)^3;
    
    jnorm = 1;
    
    maxa = floor(sqrt(X/3));

    /* Make sure that maxa = floor(sqrt(X/3)) */
    if(maxa^2 > X/3 && (maxa+1)^2 <= X/3,
        print("FAIL (Rounding error)");
	return(0);
    );

    for(a=1,maxa,
    
        sqrtXdiva = floor(sqrt(X)/a);

	/* Make sure there was no rounding error. */
	if(a^2*sqrtXdiva^2 > X,
	    print("FAIL (Rounding error)");
	    return(0);
	);

	/* e^{\pi} > 23 */
	/* y \le e^{\pi X^{1/2}/a} - 2079  \le |j(\tau)| */
	    
        y = 23^sqrtXdiva-2079;
	for(b=-a+1,a,
	    if((b^2-Delta)%(4*a) != 0, next;);
	    c = (b^2-Delta)/(4*a);
	    
	    /*
	      Note that b>-a, next two ifs check that
	      we are in the fundamental domain.
	    */	 

	    if(a>c,next;);
	    if( (a==c) && (b<0),next;);

	    /* (a,b,c) must be coprime. */
	    if(gcd(a,gcd(b,c))!=1,next;);
	    
	    
	    clno = clno + 1;
			
	    jnorm = jnorm * max(m,y);		 
	);
    );

/* Sanity check. Make sure our n matches with Pari's computation */
/* Used only for debugging purposes */

/*
    if(qfbclassno(Delta) != clno || qfbclassno(Delta,1) != clno,
        return("FAIL (internal error)");
    	return(0);
    ); */
    return(jnorm);
}


find_small_norm(B) = {
    local(X,Delta);
    gettime();
    for(X=4,B,
	if(X%10000==0,print("Reached ",X));
	Delta = -X;
	if((Delta%4!=0) && (Delta%4)!=1,next);
	if(abs(get_log_norm(Delta))<=1,print(Delta, " may come from a singular unit"););
    );
    print("Done after ",gettime(),"ms.");
    
}

{
    find_small_norm(300000);
}