#ifndef SUBMISSION_H
#define SUBMISSION_H

struct submission
{
	// An index to denote the code family
	int family;
	// Number of code lengths in the simulation
	int ncodes;
	// The lengths of the code in the specified family.
	int *lengths;

	// Code distances
	int *distances;
	
	// Noise model type -- see noise.h for different noise models.
	int model;
	
	// Names of the noise model parameters. Refer to noise.h and notations.h for the names of different parameters of different noise models.
	char **pnames;
	
	/*
		Values of the different parameters specifying the noise model. The ranges are specified as a triple -- low, max, step.
		Refer to the variable params of noise.h.
		params[0][0] = Number of different parameters
		params[j] = Range for parameter j-1.
	*/
	double **params;
	
	// Number of source files. A source file is either a .c or a .h file that is required by makefile and perf.txt file.
	int nsource;

	/* 
		Names of the source files. 
		All the source files must be listed in a text file called "sources.txt".
		The first line of "sources.txt" must be the number of source files.
	*/
	char **source;
	
	/*
		Names of the files specifying the tilings for the respective code lengths of the specified family.
		lattices[j] -- consists of two strings. One is the file name containing the edge information and another is the filename containing the face information, of the underlying tiling.
	*/
	char ***lattices;

	// name of the job assigned to the submission -- A user input is propmpted to fill in this information.
	char *job;
	
	/* 
		Name of the host-cluster on which the simulation must be run.
		The host name is by default the name of the mammount host in which the c file is being run. This is obtained by "echo $BQMAMMOUTH"
	*/
	char *host;

	// name of the queue on which the simulation task must be held -- A user input is propmted to fill in this detail.
	char *queue;

	// Walltime required per node, in hours.
	int wall;

	// The current (system) time in milliseconds. This string is used to stamp the input, run, preprocessor and post processor files.
	char *timestamp;

	// Name of the file containing the inputs to the simulation
	char *input;

	// Name of the file containing the command to be executed
	char *run;

	// Name of the file containing the post-processing commands
	char *post;

	// Name of the file containing the pre-processing commands
	char *pre;

	// Email address to which job termination status and notifications are sent out.
	char *email;

	// Name of the log-file
	char *logfile;

	// Maximum number of nodes that can be run simultaneously
	int concurrent;

	// Total number of nodes required by the simulation
	int njobs;

	// Percentage of total resource allocation used by the simulation
	float usage;

	// file containing the submission instructions
	char *instructions;
};


/*
	Make a submission
	1. Prepare the values for the various variables in the submission data structure, in PrepareSubmission(...).
	2. Prepare the list of commands that must be executed before starting the job.
		These are executed right away after running bqsubmit and are in Preprocessor(...).
	3. Prepare the contents of the input file as well as a perf.txt file.
	4. Prepare a file called run.sh containing the main command to be launched for the simulation, at every node.
	5. Prepare the post processor file called post.sh which contains instructures to seggregate the output data from various BQ folders.
	6. Prepare a bqsubmit.dat file with the details of the submission
	7. Check if everything is OK with the bqsubmit file and all the included files
*/
extern void Submit(char *instructions, int family, int **lengths, int model, double **params, long stats, char *jobname, char *queue, int wall);

#endif /* SUBMISSION_H */
