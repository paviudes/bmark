#ifndef MCRESULT_H
#define MCRESULT_H

struct mcresult
{
	// This structure is used to store the results of the Montecarlo simulations for a given noise rate.
	/*
		fails (3 x 1) long -- number of trials that resulted in a decoding failure for every error type.
			  fails[0] -- number of failures with errors of type X or Z
			  fails[1] -- number of failures with errors of type X
			  fails[2] -- number of failures with errors of type Z
	*/
	long *fails;

	/*
		trials (3 x 1) long -- number of Montecarlo trials used to estimate the failure rate for every error type.
			  trials[0] -- number of trials for errors of type X or Z
			  trials[1] -- number of trials for errors of type X
			  trials[2] -- number of trials for errors of type Z
	*/
	long *trials;

	/*
		rates (3 x 1) double -- failure rates for errors of every type
			  rates[0] -- failure rate for errors of type X or Z
			  rates[1] -- failure rate for errors of type X
			  rates[2] -- failure rate for errors of type Z
	*/
	double *rates;

	/*
		bars (3 x 1) double -- error in the estimation of failure rates.
			 bars[0] -- error in the estimation of failure rate with X or Z errors.
			 bars[1] -- error in the estimation of failure rate with X errors.
			 bars[2] -- error in the estimation of failure rate with Z errors.
	*/
	double *bars;

	/*
		running (M x 1) double -- checkpoints at which the estimated running average by Monte Carlo needs to be stored.
			running[0] -- (2 x 1) array of double
				running[0][0] = number of check points.
				running[0][1] = index of the next check point at which the running average needs to be registered.
			running[i + 1] -- (4 x 1) array of double
				running[i + 1][0] -- the Monte Carlo trial at which the i-th running average needs to be stored.
				running[i + 1][1] -- the running average of X or Z errors, at the i-th checkpoint.
				running[i + 1][2] -- the running average of X errors, at the i-th checkpoint.
				running[i + 1][3] -- the running average of Z errors, at the i-th checkpoint.
	*/
	double **running;

	double runtime; // time taken to estimate the failure rate
	int *pass; // flag to indicate if the estimation of the failure rates are reliable.
};

// Initialize the parameters of the result structure
extern void InitializeResult(struct mcresult *pmcr, long *breakpoints);

/*
	Reset some parameters of the result structure which are to be recorded newly for the next Montecarlo simulation.
	Fill the rates, trials, fails and bars arrays with zeros.
*/
extern void ResetResults(struct mcresult *pmcr);

// Free memory allocated to the results structure
extern void FreeResult(struct mcresult *pmcr);

// Print all the values in the result of the current simulation
extern void PrintResult(struct mcresult *pmcr);


/*
	Compute the running estimated average at regular intervals in the Monte Carlo simulation.
	The intervals (called breakpoints) are given in the first column of pmcr->running.
	The next three columns contain the X & Z, X and Z failure rates.
	The index of the next breakpoint for registering the running averge is given in pmcr->nextbp
*/
extern void UpdateRunningAverage(long stat, struct mcresult *pmcr);

#endif /* MCRESULT_H */