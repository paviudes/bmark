#ifndef INFER_H
#define INFER_H INFER_H

/*
	Infer the physical noise rate that gives a particular logical failure rate, for a given code
	Inputs:
		G: structure tiling
		errortype: integer to denote the type of errors for which the physical noise rate must be found.
			errortype = 0 for both X and Z type errors
						1 for X type errors
						2 for Z type errors
		scaling: real number that represents the target logical failure rate for which the corresponding physical noise rate must be found.
				 The logical failure rate f is interpretted as: f = 10^scaling.
				 So, scaling is always a negative number.
	Output:
		double --> physical noise rate for which the target logical error rate is schieved.

	Algorithm:
		1. A few points are found on the performance curve by evaluating the logical error rate at a few values of the physical noise rate.
			These values are stored as points on a log-log scale in a set called the training set.
		2. Now, a best fit line is obtained for the points (on log-log scale) in the training set.
		4. Sometimes the points in the training set might be in a region where the performance curve has a small slope.
			In this case, the interpolation in Predictor(...) will yield is very small physical noise rate.
			To avoid this, after a training set has been built we will verify if the prediction is actually less than the target logical failure rate itself.
			If so, then we go back to step 1, redesign the training set by choosing points further down the performance curve, where holefully the curve is steeper.
		3. By interpolating this line to the target logical failure rate, we compute the corresponding physical noise rate.
			The computed physical nosie rate is called a prediction and it is computed by the function Predictor(...).
		4. After a prediction is obtained, we validate this prediction using Validator(...) to check if it indeed gives the target logical failure rate (up to some tolerance level).
		5. If so, we terminate.
		6. Else, we add the prediction to the training set along with the corresponding logical failure rate computed for this prediction and start over again from step 2.

	default values for constants:
		rtol:
			positive real number to denote the tolerance with which we accept a logical failure rate to be equal to the target failure rate.
			the tolerance is on a log-log scale and it is a relative tolerance.
			Presently, it is set to 0.1.
			This means that a target logical failure rate of 4 is assumed to be achieved whenever the measured failure rate is 10^-f where
			10^(-f * (1 + 0.1)) < target < 10^(-f * (1 - 0.1))
			10^ (-1.1 f) < target < 10^(-0.9 f)
		maxTrainingSetSize:
			positive integer to denote the maximum number of points to be on the training set during the entrie simulation.
			Presently, it is set to 0.1.
		maxLearningSteps:
			maximum number of initial points to sample on the performance curve before a prediction is made.

*/
extern double Infer(struct tiling G, int errortype, double target);

// Predictor function, also to be used by RBIMInfer in rbim.c
/*
	Given a set of (x,y) data and a new y-point that may not be in the data, predict the new x-value that would result in the given y-point, if the data were to follow the trend of its best fit line.
	Input: trainingset -- array of N+1 double pointers where N is the total number of points in the (x,y) data.
			trainingset[0][0] = N
			trainingset[j] = array of two real numbers, the x data followed by the y data for each point in the set.
		   test -- real number denoting the new y-point
	Output: real number denoting the predicted x-point.
*/
extern double Predictor(double **trainingSet, double test);

#endif /* INFER_H */
