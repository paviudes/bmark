#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "arrays.h"
#include "rgens.h"

double Random(const char *distribution, int nspecs, double *specifications){
  // Sample according to a given distribution
  double urand = ((double)rand()/(double)RAND_MAX);
  double randnum = urand;
  int i;
  if (strncmp(distribution, "uniform", 7) == 0)
    randnum = urand;

  else if (strncmp(distribution, "exponential", 11) == 0)
    randnum = (-1) * (log(1 - urand)) * specifications[0];
  
  else if (strncmp(distribution, "normal", 6) == 0){
    double sum = 0, normrand, *urands = malloc(sizeof(double) * 2);
    do{
      urands[0] = ((double)rand()/(double)RAND_MAX);
      urands[1] = ((double)rand()/(double)RAND_MAX);
      sum = pow(urands[0], 2) + pow(urands[1], 2);
    }while (sum >= 1);
    normrand = urands[0] * sqrt((-2) * log(sum)/sum);
    randnum = normrand * specifications[1] + specifications[0];
    free(urands);
  }
  
  else if (strncmp(distribution, "shuffle", 6) == 0){
    // Shuffle an array using Fisher-Yates shuffle algorithm
    // https://en.wikipedia.org/wiki/Fisher–Yates_shuffle#The_modern_algorithm
    double temp;
    int randint;
    for (i = nspecs - 1; i >= 0; i --){
      randint = (int) (((double)rand()/(double)RAND_MAX) * i);
      temp = specifications[i];
      specifications[i] = specifications[randint];
      specifications[randint] = temp;
    }
  }

  else if (strncmp(distribution, "custom", 6) == 0){
    double *cumulative = malloc(sizeof(double) * nspecs);
    FillDoubleArrayWithConstants(cumulative, nspecs, 0);
    cumulative[0] = specifications[0];
    for (i = 1; i < nspecs; i ++)
      cumulative[i] = cumulative[i - 1] + specifications[i];
    randnum = (double) (nspecs - 1);
    for (i = 0; i < nspecs; i ++){
      if (urand < cumulative[i]){
        randnum = i;
        break;
      }
    }
    free(cumulative);
  }
  else if (strncmp(distribution, "cumulative", 10) == 0){
    randnum = (double) (nspecs - 1);
    for (i = 0; i < nspecs; i ++){
      if (urand < specifications[i]){
        randnum = i;
        break;
      }
    }
  }
  else
    fprintf(stderr, "\033[2mDistribution type '%s' is unknown. Assuming a uniform distrubution instead.\033[0m\n", distribution);
  return randnum;
}

void RandBinWeight(int nbits, int weight, int* sequence){
  // Generate a random binary sequence of a given weight.
  // If the weight is not specified, i.e., it is -1, then generate a random weight.
  if (weight < 0)
    weight = (int) floor((double) rand()/(double) RAND_MAX * (double) nbits);
  // Choose k locations at random without collisions and place 1's
  int count = 0, loc = 0;
  while(count < weight){
    loc = (int) floor((double) rand()/(double) RAND_MAX * (double) nbits);
    if (sequence[loc] == 0){
      sequence[loc] = 1;
      count ++;
    }
  }
}

void BinomialBitString(int *result, int nbits, int *frozen, double prob){
  // Produce a n-bit string where the value of each bit is sampled according to a binomial distribution with given parameters.
  int i, pass;
  for (i = 0; i < nbits; i ++){
    pass = 0;
    if (frozen == NULL)
      pass = 1;
    else
      if (frozen[i] == 0)
        pass = 1;
    if (pass == 1)
      if (Random("uniform", 0, NULL) < prob)
        result[i] = 1;
  }
}

#if 0
int main()
{
  /*
    Testing gsl poisson random number generator
    to run this file, do
    gcc rgens.c arrays.c memory.c -o rgens.o
  */

  // srand(time(NULL));
  const gsl_rng_type *rng;
  gsl_rng *gen;

  gsl_rng_env_setup();

  rng = gsl_rng_default;
  gen = gsl_rng_alloc(rng);

  clock_t start = clock();
  int i, max = 100;
  double *specs = malloc(sizeof(double));
  specs[0] = 10;
  fprintf(stderr, "Poisson random numbers\n");
  for (i = 0; i < max; i ++)
    fprintf(stderr, "%d\t%u\n", i, gsl_ran_poisson(gen, specs[0]));
  gsl_rng_free(gen);

  clock_t end = clock();
  fprintf(stderr, "Total runtime is %g seconds.\n", ((double) (end - start))/(double) (CLOCKS_PER_SEC));
  free(specs);
  return 0;
}
#endif
