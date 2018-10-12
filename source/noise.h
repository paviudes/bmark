#ifndef NOISE_H
#define NOISE_H

struct noise
{
	// All the parameters of a noise model
	/*
		model = 0 for uncorrelated (independent) noise model
			  > 1 for correlared noise
			  	1 -- Ball of erasure centered at some edge
			  	2 -- A cluster of erased edges formed by taking a random walk from some edge
			  	3 -- The erasure pattern corresponds to a configuration of spins on a random bond Ising model
	*/
	int model;
	
	/*
		params = Parameters that specify a noise
		The array is formatted as follows.
			params[0][0] = number of noise rates. A simulation must be performed for every noise rate.
			params[0][1] = number of parameters specifying a noise rate
			params[j] = list of noise rates that specify a noise process
						uncorrelated erasures: param[j] = qubit erasure rate
						ball erasure: params[j] = qubit erasure rate, ball radius
						walk erasure: params[j] = qubit erasure rate, steps in the walk
						Ising model erasure: params[j] = |interaction strength|, |field strength|, temperature, cooling steps
	*/
	double **params;

	// index of the current noise rate that is being simulated. This will be overwritten after every failure rate estimation.
	int current;

	// Maximum weight of an error that can be corrected, irrecpective of the error pattern. This is similar to the minimum distance of the code.
	int dist;

	// probability of errors whose weight is less than the distance. They are trivially correctable.
	double trivial;
	
	/* Probability distribution of weights of errors.
		Errors of weight less than the distance are trivially correctable. So, it is profitable to sample only those errors whose weight is larger than the distance.
		However, these errors need to be chosen in a specific way, rather than on a choose-and-reject basis. Hence we first sample a value for the weight of the error, from d to n, according to the distribution given in pn->biased.
		For i.i.d noise, the probability of selecting a weight w, given than it is greater than d, is \binom{n}{w}p^w (1-p)^(n-w)/Z where Z is a normalization that only depends on d, so it can be precomputed.
			pn->truedist = {Prob(a binary sequence of weight k) | k >= d}
			truedist[j] = probability of an error having weight j, given than j >= d.
	*/
	double *truedist;
	// Distribution for importance sampling -- Gaussian distribution over weights >= d, with mean at n/2.
	double *impdist;
	// Cumulative distribution of the importance distribution... for easy sampling
	double *impcumul;
	// Bias (correction) in importance sampling -- ratio between the true and the importance distributions.
	double *bias;
	// current value of the bias
	double currentbias;

	/*
		1. average empirical noise rate -- (number of erased qubits)/(total number of qubits), averaged over Montecarlo trials for a given noise rate.
		2. For Ising model type correlations, there are a few more properties.
			A. magnetization
			B. internal evergy
			C. number of clusters
			D. maximum cluster size
			E. Average cluster size
	*/
	double *properties;

	double *corrdist;

	// Random bond Ising model structure for designing a correlated noise model of type 4, where an erasure pattern corresponds to a congifuration of spins on the Ising model.
	struct RandomBondIsingModel *ising;
};

/*
	Initialize the parameters of a noise model
	1. A noise process can be described by one or more parameters, depending on the type of noise model chosen.
	2. The values of every parameter must be obtained by scanning a range associated to that parameter. 
	3. First we must expand all the ranges to obtain all the values of every parameter. 
	
	This is important for noise models where a noise process is described by more than one parameter.
	4. Then, we must compute all possible noise rates -- that is all possible conbinations of the expanded sets (Cartesian Product).
	5. Each such combination specifies one noise process.
*/
extern void InitializeNoise(struct noise *pn, struct tiling *pG, int dist, int model, double **prange);

/*
	Initialize the Random bond Ising model according to a lattice, with arbitrary coupling strengths, field strengths and temperature
	In this case, we make a choice of
		coupling strengths = field strengths = -0.1
		temperature = 1
*/
extern void InitializeRBIMNoise(struct noise *pn, struct tiling *pG);

/*
	Update the parameters of a noise model and assign new filenames for the results of the various montecarlo runs.
		1. Presently, the updating is effective only when noise model is of type 1 or 4.
		2. For any noise model, if the new empirical noise rate is less than the old one and for the old one we have no failures, then the failure rate for the new empirical noise rate will also be 0.
			A. if the noise model is uncorrelated, ball or erasures, the empirical noise rate is directly proportional to p, the qubit-erasure rate.
				Hence we always scan values of p from large to small.
				If the old value of p has no failures, then all subsequent values of p will result in no failures.
			B. if the noise model is Ising type
				We will scan J and B in the increasing order.
				If one of them is positive, then we know that the empirical noise rate decreases with an increase in the other quantity.
				In that case, we will assume that a failure rate 0 imples that the failure rate is 0 for all subsequent values of the other quantity.
*/
extern int UpdateNoise(struct noise *pn, struct mcresult *pmcr, char *tileName, int ne);

/*
	Construct an erasure pattern produced by an i.i.d channel, but having a certain number, w, of erased qubits.
	1. We can focus our sampling over errors whose weight is larger than the minimum distance. To do this, we need to do first pick a weight w, from the distribution {Prob(a weight k binary string) | k >= d} where Prob(…) is the Binomial distribution. So then the probability of picking w will be \binom{n}{w}p^w (1-p)^(n-w)/Z where Z is a normalization that only depends on d, so it can be precomputed.
	2. Now we are guaranteed that all errors we sample are not trivially correctable.
	3. Now, once we do MonteCarlo with this sampling and estimate a decoding success rate (or a failure rate), s, we know that
			S = Z + (1 -Z) s
	where S is the true decoding success rate had we included the trivial errors (that are of weight less than d). So once we have s from Monte-Carlo, we can report S given by the above formula.
*/
extern void ImportanceIID(int *erasure, struct tiling *pG, int distance, double *distribution);

// Free a noise structre
extern void FreeNoise(struct noise *pn);

// Initialize the distributions required for importance sampling. For now, it only applies to uncorrelated noise model (model 1).
void InitImportanceSampling(struct noise *pn, struct tiling *pG);

/* 
	Print the values of the noise parameters, that correspond to the current simulation when the function is called.
	"Current" noise parameters are in (pn->params)[pn->current]
*/
extern void PrintNoise(struct noise *pn);

#endif /* NOISE_H */
