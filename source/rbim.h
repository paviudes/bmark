#ifndef RBIM_H
#define RBIM_H

struct RandomBondIsingModel
{
	/* 
	Define a random bond Ising model.
	1. cToS -- Couplings to sites array which indicates which sites are coupled, for each coupling index.
	2. sToC -- Sites to complings array which indicates which couplings a site is involved in, for every site.
	3. interactions -- Coupling strength for each coupling
	4. fields -- Strength of magnetic field at each site
	5. config -- configuration on spins.
	 */
	int nSites;
	int range;
	int coordination;
	int **sToC;
	int **cToS;
	int *config;
	long nCouplings;
	double temperature;
	double *interactions;
	double *fields;
	double energy;
	double maxClusterSize;
	double nclusters;
	double avgClusterSize;
	double magnetization;
	double internalEnergy;
	double twopoint;
};

// Initialize a random bond Ising model
extern void InitRBIM(struct RandomBondIsingModel *prbim);

// Free memory allocated to the structure
extern void FreeRBIM(struct RandomBondIsingModel *prbim);

// Print the details of a random bond Ising model object
extern void PrintRBIM(struct RandomBondIsingModel *prbim, int isDetailed);

/*
Update a random bond Ising model by assigning new coupling strengths to interactions and few field strengths
The parameters are ordered as follows.
uniform interaction strength, uniform field strength, temperature
*/
extern void UpdateRBIM(struct RandomBondIsingModel *prbim, double interact, double field, double temp);

// Compute the energy of a configuration of a random bond Ising model.
extern void Energy(struct RandomBondIsingModel *prbim);

// Initialize a random bond Ising model using the specifications of a tiling that defined an error correcting code
extern void IntializeRBMWithCode(struct RandomBondIsingModel *prbim, struct tiling* pG, double coupStrength, double fieldStrength, double temperature);

/*
	Cool a random bond Ising model, using the following procedure.
	1. Pick a random site.
	2. Compute the energy cost associated to flipping a spin at this site.
	3. If the energy cost is negative (i.e., the spin flip lowers the energy), make the spin flip
	4. Else, make the spin flip, with probability exp(- energy cost).
*/
extern void Cool(struct RandomBondIsingModel *prbim, int nsweeps);

extern void PhysicalProperties(struct RandomBondIsingModel *prbim);

// decide if a simulation can be skipped based on the previous parameters
/*
	The order of parameters is: J, B, T and l (cooling steps)
	Skip rules:
	if (new J == old J){
		For J < 0, there is no advanced skipping rule.
		if (new B < old B){
			this means that B is negative
			if all fails for old B then all fails for the new B as well.
		}
		if (new B > old B){
			this means that B is positive
			if no fails for old B then no fails for the new B as well.
		}
	}
	if (new B == old B){
		For J < 0, there is no advanced skipping rule.
		if (newB is positive){
			if (new J > old J){
				if no fails for old J then no fails for new J as well.
			}
		}
		if (new B is negative){
			if (new J > old J){
				if all fails for old J then all fails for new J as well.
			}
		}
	}
*/	
extern int RBIMSkipTest(double *oldparams, long* fails, long*trials, double *newparams);

// Infer the parameters of the Ising model using mean-field theory
extern double MeanFieldInfer(double *available, int coordination, int toinfer);

// The mean field back-calculation does not work well if J > meanfieldValidity, which in this case is 0.2.
// In this case, we resort to making a few dry-runs and learning from the accumulated data.
// The learning set is simply drawing a best fit line.
extern void RBIMInfer(struct RandomBondIsingModel *prbim, double **params, int *sizes, int toinfer);

#endif /* RBIM_H */
