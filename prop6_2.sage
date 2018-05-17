### Effective Bounds for Singular Units     
### Yuri Bilu, Philipp Habegger, Lars Kuehne 

### Sage code to verify computations in Proposition 6.2

### Produces an AssertionError if it cannot verify the claimed
### bounds for c1, c2, c3, and c4.

### Tested using Sage 7.3

def Prop62():

    # The four constants that appear in Proposition 6.2 
    c1 = RIF(0.712)
    c2 = RIF(1.010)
    c3 = RIF(2.598)
    c4 = RIF(2.267)


    # A high precision computation with sage yields that
    # lambda0 = 1/zeta(2) = 0.607927101854026...
    # lambda1 = -2*zetaderiv(1,2) / zeta(2)^2 + (2*euler_gamma-1)/zeta(2) = 0.786872460166245...
    # Intervals containing the constants lambda0, lambda1
    lambda0 = RIF(0.607927101854026,0.607927101854027)  
    lambda1 = RIF(0.786872460166245,0.786872460166246)
    

    S = RIF(0)

    #   Loop up-to 2*10^7 
    for n in xrange(1,20000000+1) :

        # successively computer the sum over 2^{\omega(n)} 
        S = S + 2^(len(prime_divisors(n)))

	# Progress update
        if n % 10000 == 0 : 
            print('Reached '+str(n))

	# Approximate \lambda_0 n \log(n) + \lambda_1
        gn = lambda0*RIF(n)*log(RIF(n)) + lambda1*RIF(n)

	# Approximate \lambda_0 (n+1) \log(n+1) + \lambda_1
        gnp1 = lambda0*RIF(n+1)*log(RIF(n+1)) + lambda1*RIF(n+1)    
	if n >= 2 : 
            assert((RIF(S)-gn)/sqrt(RIF(n)) <= c1)
            assert((gnp1-RIF(S))/sqrt(RIF(n)) <= c2)
    
        if n >= 40000 :
            assert((RIF(S)-gn)*log(RIF(n))/sqrt(RIF(n)) <= c3)
            assert((gnp1-RIF(S))*log(RIF(n))/sqrt(RIF(n)) <= c4)


Prop62()