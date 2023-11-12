#ifndef __LAPINS__
#define __LAPINS__

#include <stdio.h>
#include <math.h>
#include <stdlib.h> 
#include <string.h>
#include <time.h>
#include "mersenne_twister.h"


#define DURATION 10
#define NUMRABBITSTART 2
#define DEBUG 0


typedef struct Rabbit
{
    long long int males[16]; // Arrays in which to store the numbers of rabbits, each cell representing the age of the rabbits       
                            // [0] -> baby rabbits, [1] -> 1-year-old rabbits, ...
    long long int females[16];
}Rabbit;

int RabbitFibonacciSimulation(int n);
int generateRandomGender(void);
int calculateSurvivalRate(int ageOfRabbit);
double generateStandardGaussian(void);
int generateReproductionCount(void);
int generateNumberOfBabies(void);
long long int getRabbitsAdultFemales(Rabbit *rabbits);
int hasAdultMale(Rabbit *rabbits);
void simulateOneYear(Rabbit *rabbit);
void simulateRealisticRabbitPopulation(int numStartBaby, int numStartAdult, Rabbit *rabbit, int duration, unsigned long long int * meanF, unsigned long long int * meanM);
unsigned long long int simulateMultipleExperiments(int numStartBaby, int numStartAdult, Rabbit *rabbit, int duration, int numExp);

#endif
