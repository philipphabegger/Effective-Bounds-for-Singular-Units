/*
    Effective Bounds for Singular Units     
    Yuri Bilu, Philipp Habegger, Lars Kuehne 

    PARI code to verify claims of Lemma 6.4.

    Tested using PARI 2.9.3.

    Prints "FAIL" if the claims could not be verified.
    Prints a success message otherwise.     
*/

get_abundants(B) = {
    local(n,t,sigma0,sigmam1);

    gettime();

    /* Verify numerical claim in Lemma 6.4(i) */
    if(sigma(21621600,-1) != 3472/715,
        print("FAIL.");
    	return(0);
    );
    
    
    for(n=1,B,
	/* Compute \sigma_0(n) */
        sigma0 = sigma(n,0);

	/* Computer \sigma_{-1}(n) which equals \sigma_1(n)/n */
        sigmam1 = sigma(n,-1);

	/* Verify inequality in Lemma 6.4(i) */
	if(sigmam1 > 3472/715,
	    print("FAIL.");
	    return(0);	    
	);

	/* Verify inequality Lemma 6.4(ii) */
	/* We take the fourth power to avoid floating point */
	/* arithmetic. All computations are with rationals */
	/* and therefore rigorous. */
	
	if(sigma0^4 > (17/2)^4 * n,
	    print("FAIL.");
	    return(0);
	);        
    );

    t = gettime();
    print("Success up to ",B," after ",t,"ms.");
    return(1);
}

{
    get_abundants(32*10^6);
}