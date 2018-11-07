/*
    No singular modulus is a unit
    Yuri Bilu, Philipp Habegger, Lars Kuehne 

    C code to verify computations in Lemma 7.2.
    This is an implementation of Algorithm 1.
    
    If termination is successful, an upper bound on 
    C_epsilon(Delta) for Delta <= DISCRIMINANT_MAX is printed.
   
    The program aborts unsuccessfully if
        (1) not sufficient memory space can be reserved,
	(2) the macro DISC_BLOCK_SIZE is not devisible by 4,
	(3) the long long int type is less than 8 bytes long,
	(4) the floating point computation of minc and maxc in
	the function count_tau produce wrong results. 
    In each case the reason for abortion is printed. If (1) 
    occurs, DISC_BLOCK_SIZE should be lowered. Its default value 
    10^10 needs about 5 GB free memeory. (3) is a violation of the 
    C99 and upwards standard and overcautious. (4) should not occur 
    in the usual range, but safeguards the validity of the output.

    Using GCC, compile with "gcc algorithm1.c -lm"    

    Tested using GCC 5.4.0 on Ubuntu 16.04LTS
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

// Define macro LEMMA72i to check Lemma 7.2(i). If this macro is
// undefined, then the program checks Lemma 7.2(ii).

#define LEMMA72i

// MINA_FACTOR is floor of  1000/(1+\sqrt{3}\eps + \eps^2)
// MINB_FACTOR is floor of (1-2\eps)1000

#ifdef LEMMA72i

// Lemma 7.2(i)
// The standard block size, must be divisible by 4. 
// With this choice, the program requires about 5GB RAM. 
#define DISC_BLOCK_SIZE (long long int) 10000000000  // 10^10

#define DISCRIMINANT_MAX (long long int) 10000000000 // 10^10
#define MINA_FACTOR (long long int) 998              // \eps = 0.001
#define MINB_FACTOR (long long int) 998              // \eps = 0.001
#define LEMMA_MESSAGE "Lemma 7.2 (i)"

#else 

// Lemma 7.2(ii)
// The standard block size, must be divisible by 4. 
// With this choice, the program requires about 50MB RAM. 
#define DISC_BLOCK_SIZE (long long int) 10000000  // 10^7

#define DISCRIMINANT_MAX (long long int) 10000000 // 10^7

#define MINA_FACTOR (long long int) 993           // \eps = 0.004
#define MINB_FACTOR (long long int) 992           // \eps = 0.004
#define LEMMA_MESSAGE "Lemma 7.2 (ii)"

#endif


int count_tau(unsigned char *, long long int minX, long long int len);

int count_tau(unsigned char *counter, long long int minX, long long int len) {
  long long int a, b, c, X, minc, minb, mina, maxc, roundmax, maxX;

  maxX = minX + len - 1;
  
  minc = floor(sqrt((long double) minX)/2);
  maxc = floor(sqrt((long double) maxX));  

  // We check the above rounding to preclude floating point errors.
  // We make sure that minc \le minX^{1/2}/2.
  // We also make sure that maxc > maxX^{1/2}-1, which implies maxc \ge \floor{maxX^{1/2}}
  if(4*minc*minc > minX || (maxc+1)*(maxc+1) <= maxX){
    printf("%lld %lld %lld %lld\n",minc,minX,maxc,maxX);
    printf("Rounding error, terminating.\n");
    return 0;
  }
  
  // We check that minX \equiv 1 \mod 4
  // This is true when called from main below.
  if(minX%4!=1) {
    printf("minX is not 1 mod 4, terminating.\n");
    return 0;
  }

  if(minc == 0)
    minc = 1;
  
  roundmax = 0;

  for (c = minc; c <= maxc; c++) {
    // Integer division truncates towards zero according to C89 standard (and upwards)
    mina = MINA_FACTOR*c/1000; 

    for (a = mina; a <= c; a++) {
      // Integer division truncates towards zero according to C89 standard (and upwards).
      minb = MINB_FACTOR*a/1000;

      // We only count conjugates with positive real part.
      // We need to multiply the final result by 2.
      for (b = minb; b <= a; b++) {

	// X is the absolute value of the discriminant of (a,b,c).
	X = 4*a*c - b*b;

	// Skip cases where we are out of bounds.
	if(X > maxX || X < minX){
	  continue;
	}
	
	// We cut memory usage in half by using that X is equivalent to 0,3 mod 4.
        // Since minX is equivalent to 1 mod 4, we have tmplocation = 2 + 4k or 3 + 4k.
	// location will be 2k or 1+2k, respectively.

	long long int tmplocation = X - minX;
	long long int location;
	
	if(tmplocation%4==2){
	  location = (tmplocation - 2) >> 1;
	}
	if(tmplocation%4==3){
	  location = (tmplocation - 1) >> 1;
	}

	// Either case below should never occur. Here for debugging purposes. 
	if(tmplocation%4==0 || tmplocation%4==1 || location >= DISC_BLOCK_SIZE/2){
	  return -1;
	}
		
	unsigned char tmp = counter[location];

	// Make sure that there is no overflow.
	if(tmp == 254)
	  {
	    return -1;
	  }
	
	/* 
	   Count twice to accommodate negative real part. 
	   We may overcount, which is acceptable since we are
	   only interested in an upper bound for C_epsilon(Delta)
	*/
	counter[location] = tmp + 2;

	// Have we found a new record?
	if(counter[location] > roundmax){
	  roundmax = counter[location];
	}        
      }
    }
  }
  return roundmax;
}


int main(void) {
  long long int n,loops,count,globalcount;
  unsigned char *counter;
  clock_t t1,t2;

  t1 = clock();
  
  printf("Calculations needed in %s.\n",LEMMA_MESSAGE);
  if(sizeof(long long int) < 8){
    printf("We require long long int to be at least 8 bytes long, terminating.");
    return(0);
  }

  if(DISC_BLOCK_SIZE%4 != 0){
    printf("Block size must be divisible by 4, terminating.");
    return(0);
  }

  printf("Going up to absolute discriminant %lld\n", DISCRIMINANT_MAX);

  printf("Block size: %lld\n",DISC_BLOCK_SIZE);
  
  printf("Reserving %lld bytes of RAM.\n",sizeof(unsigned char)*DISC_BLOCK_SIZE);

  // Allocate memory
  counter = (unsigned char *) malloc(sizeof(unsigned char)*DISC_BLOCK_SIZE/2);

  // Exit gracefully if not enough memory available.
  if(counter == NULL){
    printf("Error allocating memory, terminating.\n");
    exit(0);
  }

  // Compute number of loops
  loops = DISCRIMINANT_MAX / DISC_BLOCK_SIZE;

  // If this quotient isn't an integer, round up.
  if(DISCRIMINANT_MAX % DISC_BLOCK_SIZE != 0)
    loops++;
   
  globalcount = 0;
 
  for(n=1;n<=loops;n++){
    // Set content of memory to 0.
    memset(counter,0,sizeof(unsigned char)*DISC_BLOCK_SIZE/2);

    long long int block_start = DISC_BLOCK_SIZE*(n-1)+1;
    long long int block_len = DISC_BLOCK_SIZE;

    // Truncate block at the end.
    if(block_start + block_len - 1 > DISCRIMINANT_MAX){
      block_len = DISCRIMINANT_MAX - block_start + 1;
    }

    printf("\nChecking absolute discriminant in range %lld to %lld\n", DISC_BLOCK_SIZE*(n-1)+1, block_start+block_len-1);

    count = count_tau(counter, block_start, block_len);

    // Exit on failure
    if(count == -1){
      printf("FAIL\n");
      exit(0);
    }
 
    printf("Found C_epsilon(Delta)<=%lld in current range.\n\n", count);

    if(count > globalcount){
      globalcount = count;
    }
  }
  t2 = clock();

  // Measure the time it took to do computation. 
  // For debugging purposes only.
  printf("Done after %fms.\n", ((double) (t2-t1))*1000/CLOCKS_PER_SEC);
  printf("Found C_epsilon(Delta)<=%lld in range 1 to %lld.\n", globalcount, DISCRIMINANT_MAX);  
  free(counter);
	
  return EXIT_SUCCESS;
}
