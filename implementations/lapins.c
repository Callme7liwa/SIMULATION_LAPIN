/** ------------------------------------------------------------------- *
 * @todo       Rabbit colony simulation algorithm                      *
 * @author     Ayoub Seddiki                                           *                                    *
 *                                                                     *
 * Updates:                                                            *
 *     1. The simple rabbit swarm algorithm has been updated from      *
 *        a fictional version to a Fibonacci algorithm in order to     *
 *        simplify its complexity.                                     *    *
 *     2. The calculation has been updated to be performed on a        *
 *         yearly basis.                                               *
 * ------------------------------------------------------------------- */

#include "./../headers/lapins.h"


/** ------------------------------------------------------------------- *
 * @fn       RabbitFibonacciSimulation                                 *
 *                                                                     *
 * @brief    Calculates the number of rabbit pairs after n years       *
 *           using a Fibonacci-based series.                           *
 *                                                                     *
 * @param    n {int}    Duration of this algo.                         *
 *                                                                     *
 * @return   The total number of rabbit pairs after n years            *
 *           following a Fibonacci series.                             *
 * ------------------------------------------------------------------- */
int RabbitFibonacciSimulation(int n)
{
     // Initialize variables
    int currentPopulation, previousPopulation, temp;

    // Initialize the first two years
    currentPopulation = 1;
    previousPopulation = 1;

    // Special case for the first two years
    if (n == 1 || n == 2)
    {
        printf("Year: %d\tNumber of rabbit pairs: 1\n", n);
    }

    // Calculate the population for subsequent years using Fibonacci sequence
    for (int year = 3; year <= n; year++)
    {
        // Tumbling and adding method to update populations
        temp = currentPopulation;
        currentPopulation += previousPopulation;
        previousPopulation = temp;

        // Print the result for each year
        printf("Year: %d\tNumber of rabbit pairs: %d\n", year, currentPopulation);
    }

    return currentPopulation;
}

/** -------------------------------------------------------------------      *
 * @fn         generateRandomGender                                          *
 *                                                                           *
 * @brief       Generates a random gender with equal probability.            *
 *                                                                           *
 * @return     Gender code where 1 represents male and 0 represents female.  *
 * -------------------------------------------------------------------       */
int generateRandomGender(void)
{
    // Variable to store the generated random number
    double randomValue = uniform(0,1);
    // Determine gender based on the random value
    int gender = (randomValue < 0.5) ? 1 : 0;
    return gender;
}

/** ------------------------------------------------------------------- *
 * @fn         calculateSurvivalRate                                    *
 *                                                                      *
 * @brief      Calculates the survival rate of a rabbit                 *
 *             based on age and time to sexual maturity.                *
 *                                                                      *
 * @param       ageOfRabbit {int} The age of the rabbit.                *
 *                                                                      *
 * @return     Survival state where 1 represents survival               *
 *             and 0 represents death.                                  *
 * ------------------------------------------------------------------- */
int calculateSurvivalRate(int ageOfRabbit)
{
    // Variables for chance of survival, age adjustment, and survival state
    double survivalProbability = genrand_real2();
    double ageAdjustment = 0;
    int survivalStatus = 0;

    // If the rabbit is not mature
    if (ageOfRabbit == 0)
    {
        // If it survives
        if (survivalProbability <= 0.35)
            survivalStatus = 1;
    }
    else
    {
        // If the rabbit is adult and under ten years old
        if (ageOfRabbit < 10)
        {
            // If it survives
            if (survivalProbability <= 0.6)
                survivalStatus = 1;
        }
        else
        {
            ageAdjustment = ageOfRabbit - 9;

            // If it survives with age adjustment
            if (survivalProbability <= (0.6 - (ageAdjustment * 0.1)))
                survivalStatus = 1;
        }
    }

    return survivalStatus;
}

/** ------------------------------------------------------------------- *
 * @fn        generateStandardGaussian                                                *
 *                                                                     *
 * @brief     Generates a random value from a standard                 *
 *            Gaussian distribution.                                   *
 *                                                                     *
 * @return    A random value in a standard Gaussian distribution       *
 * ------------------------------------------------------------------- */
double generateStandardGaussian(void)
{
    static int hasSpare = 0;
    static double spare;

    // If there is a spare value, return it
    if (hasSpare)
    {
        hasSpare = 0;
        return spare;
    }

    // Generate two random values in the range [-1, 1]
    double u, v, s;
    do
    {
        u = 2.0 * genrand_real2() - 1.0;
        v = 2.0 * genrand_real2() - 1.0;
        s = u * u + v * v;
    } while (s >= 1.0 || s == 0.0);

    // Calculate a random value following a standard Gaussian distribution
    s = sqrt(-2.0 * log(s) / s);
    spare = v * s;
    hasSpare = 1;
    return u * s;
}

/** ------------------------------------------------------------------- *
 * @fn        generateReproductionCount                                *
 *                                                                     *
 * @brief     Randomly generates the number of times a female          *
 *            rabbit reproduces in a year.                             *
 *                                                                     *
 * @return    The number of reproduction in normal distribution            *
 * ------------------------------------------------------------------- */
int generateReproductionCount(void)
{
     int reproductionCount;
    double randomNumber = genrand_real2() + 6;

    // Determine the reproduction count based on the random number
    if (randomNumber < 4.5)
        reproductionCount = 4;
    else if (4.5 <= randomNumber && randomNumber < 5.5)
        reproductionCount = 5;
    else if (5.5 <= randomNumber && randomNumber < 6.5)
        reproductionCount = 6;
    else if (6.5 <= randomNumber && randomNumber < 7.5)
        reproductionCount = 7;
    else
        reproductionCount = 8;

    return reproductionCount;
}

/** ------------------------------------------------------------------- *
 * @fn        generateNumberOfBabies                                                *
 *                                                                     *
 * @brief     Generate an integer between 3 and 6 with equal           *
 *            probability                                              *
 *                                                                     *
 * @return    An integer variable between 3 and 6                      *
 * ------------------------------------------------------------------- */
int generateNumberOfBabies(void)
{
    // Generate a random value in the range [0, 1) and map it to [3, 6)
    return 3 + (genrand_real2() * 3);
}


/** -------------------------------------------------------------------     *
 * @fn       updateRabbitPopulation                                         *
 *                                                                          *
 * @brief    Updates the rabbit population after one year based on their    *
 *           age and survival rate.                                         *
 *                                                                          *
 * @param    start {Rabbit *} The array of rabbits.                         *
 * ------------------------------------------------------------------- */
void updateRabbitPopulation(Rabbit *start)
{
    long long int rabbitCount;

    // All 15-year-old rabbits die
    start->females[15] = 0;
    start->males[15] = 0;

    // For each age, starting from the oldest
    for (int age = 14; age >= 0; age--)
    {
        // Update survival and move to the next age for females
        for (rabbitCount = 0; rabbitCount < start->females[age]; rabbitCount++)
        {
            if (calculateSurvivalRate(age))
            {
                start->females[age + 1]++;
            }
        }

        // Reset the current age count for females
        start->females[age] = 0;

        // Update survival and move to the next age for males
        for (rabbitCount = 0; rabbitCount < start->males[age]; rabbitCount++)
        {
            if (calculateSurvivalRate(age))
            {
                start->males[age + 1]++;
            }
        }

        // Reset the current age count for males
        start->males[age] = 0;
    }
}


/** ------------------------------------------------------------------------*
 * @fn         hasAdultMale                                                 *
 *                                                                          *
 * @brief      Checks whether there is at least one male rabbit.            *
 *                                                                          *
 * @param      rabbits {Rabbit *} Struct representing the rabbits.          *
 * @return     State of the presence of adult male rabbit [1->Yes; 0->No].  *
 * ------------------------------------------------------------------------ */
int hasAdultMale(Rabbit *rabbits)
{
    int currentAge = 1;

    // Iterate through ages until an adult male is found or reach the maximum age
    while (currentAge < 16 && rabbits->males[currentAge] == 0)
    {
        currentAge++;
    }

    // Return 1 if an adult male is found, otherwise return 0
    return (currentAge == 16) ? 0 : 1;
}

/** ------------------------------------------------------------------- *
 * @fn         getAdultFemales                                         *
 *                                                                     *
 * @brief      Calculates the total number of adult female rabbits.    *
 *                                                                     *
 * @param      rabbit {Rabbit} The rabbit structure.                   *
 * @return     The total count of adult female rabbits.                *
 * ------------------------------------------------------------------- */
long long int getRabbitsAdultFemales(Rabbit *rabbits)
{
    long long int adultFemalesCount = 0;

    // Iterate through different age groups to accumulate adult females
    for (int age = 1; age < 16; age++)
    {
        adultFemalesCount += rabbits->females[age];
    }

    return adultFemalesCount;
}

/** ------------------------------------------------------------------------*
 * @fn       simulateOneYear                                                *
 *                                                                          *      
 * @brief    Simulates one year for the rabbit population.                  *
 *                                                                          *
 * @param    rabbit {Rabbit *} Struct representing the rabbits.             *
 * ------------------------------------------------------------------------ */
void simulateOneYear(Rabbit *rabbit)
{
    long int babyCountThisYear;
    int numberMale = 0, numberFemale=0;
    int hasAdultMaleRabbit;
    int reproductionTimeThisYear;
    int gender;
    long long int numAdultFemales;

    // Calculate the number of adult females in the rabbit population
    numAdultFemales = getRabbitsAdultFemales(rabbit);

    // Check if there is at least one adult male rabbit
    hasAdultMaleRabbit = hasAdultMale(rabbit);

    // Initialize variables for baby count
    babyCountThisYear = 0;

    /* If there is at least one adult male */
    if (hasAdultMaleRabbit)
    {
        /* For all the adult females */
        for (int i = 0; i < numAdultFemales; i++)
        {
            // Get the number of times a female reproduces in a year
            reproductionTimeThisYear = getRabbitsAdultFemales(rabbit);

            /* For all the litters */
            for (int j = 0; j < reproductionTimeThisYear; j++)
            {
                // Accumulate the number of babies born this year
                babyCountThisYear += generateNumberOfBabies(); 
            }
        }

        /* For all the babies born */
        for (int i = 0; i < babyCountThisYear; i++)
        {
             // Determine the gender of each baby
            gender = generateRandomGender();
            if(gender == 1)
                numberMale ++ ; 
            else
                numberFemale++;
        }
    }

    // Update the rabbit population based on deaths
     updateRabbitPopulation(rabbit);

    // Add the newborns to the population
    rabbit->males[0] = numberMale;
    rabbit->females[0] = numberFemale;
}

/** ------------------------------------------------------------------- *
 * @fn       simulateRealisticRabbitPopulation                          *
 *                                                                      *
 * @brief    An algorithm for realistic simulation of rabbit            *
 *           populations.                                               *
 *                                                                      *
 * @param    numStartBaby {int} Number of baby rabbit couples at the    *
 *           start.                                                     *
 * @param    numStartAdult {int} Number of adult rabbit couples at the  *
 *           start.                                                     *
 * @param    rabbit {Rabbit *} Struct representing the rabbit           *
 *           population.                                                *
 * @param    duration {int} Duration of the simulation.                 *
 * @param    meanF {unsigned long long int *} Array to count how many   *
 *           female rabbits there are for each year.                    *
 * @param    meanM {unsigned long long int *} Array to count how many   *
 *           male rabbits there are for each year.                      *
 * ------------------------------------------------------------------- **/
void simulateRealisticRabbitPopulation(int numStartBaby, int numStartAdult, Rabbit *rabbit, int duration, unsigned long long int *meanF, unsigned long long int *meanM)
{
    long long int numMales;
    long long int numFemales;
    int age;

    // Initialize the arrays for each age
    for (int i = 0; i < 16; i++)
    {
        rabbit->females[i] = 0;
        rabbit->males[i] = 0;
    }

    // Initialize the baby rabbits
    for (int i = 0; i < numStartBaby; i++)
    {
        rabbit->males[0]++;
        rabbit->females[0]++;
    }

    // Initialize the adult rabbits with age from 1 to 4
    for (int i = 0; i < numStartAdult; i++)
    {
        age = 1 + 3 * genrand_real1();
        rabbit->males[age]++;
        rabbit->females[age]++;
    }

    // Simulation for every year
    for (int year = 1; year <= duration; year++)
    {
        simulateOneYear(rabbit);

        // Display the rabbit counts for the year
        numMales = 0;
        numFemales = 0;

        // Count all the rabbits
        for (int i = 0; i < 16; i++)
        {
            numMales += rabbit->males[i];
            numFemales += rabbit->females[i];
        }

        if(DEBUG)printf("After year %d, there are %lld males and %lld females\n", year, numMales, numFemales);
        meanF[year - 1] += numFemales;
        meanM[year - 1] += numMales;
    }

    if (DEBUG) printf("\n");
}

/** ------------------------------------------------------------------- *
 * @fn       simulateMultipleExperiments                                *
 *                                                                      *
 * @brief    Simulates a certain number of rabbit population growth.    *
 *                                                                      *
 * @param    numStartBaby {int} Number of couples of baby rabbits.      *
 * @param    numStartAdult {int} Number of couples of adult rabbits.    *
 * @param    rabbit {Rabbit *} Struct for the rabbits.                  *
 * @param    duration {int} Duration of the simulation.                 *
 * @param    numExp {int} Number of simulations to be done.             *
 * @return   Mean for the population of rabbits.                        *
 * ------------------------------------------------------------------- **/
unsigned long long int simulateMultipleExperiments(int numStartBaby, int numStartAdult, Rabbit *rabbit, int duration, int numExp)
{
    unsigned long long int *meanF = (unsigned long long int *)malloc(sizeof(unsigned long long int) * duration);
    unsigned long long int *meanM = (unsigned long long int *)malloc(sizeof(unsigned long long int) * duration);
    unsigned long long int meanTotal = 0;

    // Initialize mean arrays
    for (int d = 0; d < duration; d++)
    {
        meanF[d] = 0;
        meanM[d] = 0;
    }

    // For each experiment
    for (int exp = 1; exp <= numExp; exp++)
    {
        if (DEBUG) printf("_________________________\n\tExperiment %d : \n", exp);

        simulateRealisticRabbitPopulation(numStartBaby, numStartAdult, rabbit, duration, meanF, meanM);

        if (DEBUG) printf("\n");
    }

    // Display means for each year
    for (int d = 0; d < duration; d++)
    {
        printf("_____________________\n *-* Means for year %d :\n\n", d + 1);
        printf("Females : %llu\n", meanF[d] / numExp);
        printf("Males : %llu\n", meanM[d] / numExp);
        meanTotal += meanM[d] + meanF[d];
    }

    // Free memory
    free(meanF);
    free(meanM);

    // Return mean population count over all experiments
    return meanTotal / numExp;
}
