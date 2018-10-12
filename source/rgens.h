#ifndef RGENS_H
#define RGENS_H

/*
  Generate a radom number according to a given distribution
  Inputs:
          distribution: char* (name of a distribution)
          nspecs: integer (number of parameters to specify the distribution)
          specifications: double array (parameters to specify the distribution)

              distribution                nspecs          specification                                       algorithm
          1.    "uniform"                    0            "NULL"                                              Built-in "double rand(..)" function.
          2.    "poisson"                    1            double array of size 1: [mean]                      https://en.wikipedia.org/wiki/Poisson_distribution#Generating_Poisson-distributed_random_variables
          3.    "exponential"                1            double array of size 1: [mean]                      Inverse sampling transform. CDF of the exponential distributio is F(v) = 1 - exp( - 1/avg * v). So, output x = - log(1 - u) * avg, where u is a uniform random variate.
          4.    "normal"                     2            double array of size 2: [mean, variance]            https://en.wikipedia.org/wiki/Marsaglia_polar_method
          5.    "custom"                     m            double array of size m: [prob(i) : 1 <= i <= m]     Construct the cumulative distribution. generate a uniform random number u. Find out the index of the distribution that u falls into, i.e, next element is greater than u while the previous is less than u.
          (where m is the number of events for which the probability distribution is defined)
  Output:
          (double) real number
*/
extern double Random(const char *distribution, int nspecs, double *specifications);

/*
  Produce a n-bit string where the value of each bit is sampled according to a binomial distribution with given parameters.
  For each bit:
  1. Pick a real number in the range [0,1] according to the unform distribution
  2. If the number is less than prob, the bit value is 1
  3. Otherwise, the bit value is 0.
*/
extern void BinomialBitString(int *result, int nbits, int *frozen, double prob);

/*
  Generate a random binary sequence of a given weight.
  If the weight is not specified, i.e., it is -1, then generate a random weight.
*/

extern void RandBinWeight(int nbits, int weight, int* sequence);

#endif /* RGENS_H */
