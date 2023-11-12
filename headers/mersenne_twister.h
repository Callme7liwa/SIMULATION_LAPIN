#ifndef __MT__
#define __MT__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Period parameters */  
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

void init_genrand(unsigned long s);
void init_by_array(unsigned long init_key[], int key_length);
unsigned long genrand_int32(void);
long genrand_int31(void);
double genrand_real1(void);
double genrand_real2(void);
double genrand_real3(void);
double genrand_res53(void);
double uniform(double a, double b); //generates a random number on [a,b]-real-interval
int uniform_int(int a, int b); //generates a random number on [a,b]-integer-interval
#endif
