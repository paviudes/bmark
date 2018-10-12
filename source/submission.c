#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <time.h>
#include <errno.h>
#include <unistd.h>
#include <string.h>
#include "memory.h"
#include "notations.h"
#include "submission.h"

#define MAX_JOBS 100

int PrepareSubmission(struct submission *psub, int fam, int **lengths, int model, double **ranges, char *jobname, char *queue, int wall, char *instructions){
	/*
		Prepare the values for the various variables in the submission data structure
		Set the queue, host, source, time stamp and the names for files containing preprossing, running, input and postprocessing information.
		1. The host name is by default the name of the mammount host in which the c file is being run.
			This is obtained by "echo $BQMAMMOUTH".
		2. Email address is by default pavithran.sridhar@gmail.com
		3. Logfile is called log.txt and it is always opened in append mode.
			It is once opened and a timestamp is added.
			The logging or copying of stdout to file is done using the tee command.
		4. All the source files must be listed in a text file called "sources.txt" and thier count must be provided in the first line of this file.
		5. The tiling information filenaes are read for the code lengths that are part of the simulation
		6. If the total number of jobs jobs is greater than 144, set the number of concurrent jobs to 144. Or else set it to the total number of jobs available.
		7. The job name, walltime and the queue are to be prompted from the user
		8. The usage is calculated in percentage.
			The total time allowance on MS is 100000 cores-hours. On MP2 it is 700000 cores-hours.
	*/
	int fi, pi;
	time_t tme = time(NULL);
	struct tm tm = *localtime(&tme);
	psub->timestamp = malloc(sizeof(char) * 100);
	sprintf(psub->timestamp, "%02d_%02d_%02d_%02d_%02d_%02d", tm.tm_mday, tm.tm_mon + 1, tm.tm_year + 1900, tm.tm_hour, tm.tm_min, tm.tm_sec);

	psub->instructions = malloc(sizeof(char) * 100);
	sprintf(psub->instructions, "%s", instructions);
	psub->pre = malloc(sizeof(char) * 100);
	sprintf(psub->pre, "pre_%s.sh", psub->timestamp);
	psub->input = malloc(sizeof(char) * 100);
	sprintf(psub->input, "launch.txt");
	psub->run = malloc(sizeof(char) * 100);
	sprintf(psub->run, "run_%s.sh", psub->timestamp);
	psub->post = malloc(sizeof(char) * 100);
	sprintf(psub->post, "post_%s.sh", psub->timestamp);

	// The host name is by default the name of the mammount host in which the c file is being run. This is obtained by "echo $BQMAMMOUTH"
	psub->host = malloc(sizeof(char) * 100);
	FILE *hfid = popen("echo $BQMAMMOUTH", "r");
	if (hfid == NULL)
		sprintf(psub->host, "host_unknown");
	else
		fgets(psub->host, 100, hfid);
	fclose(hfid);
	(psub->host)[strcspn((psub->host), "\n")] = 0;

	// Email address is by default pavithran.sridhar@gmail.com
	psub->email = malloc(sizeof(char) * 100);
	sprintf(psub->email, "pavithran.sridhar@gmail.com");

	// Logfile is called log.txt and it is always opened in append mode. It is once opened and a timestamp is added.
	// The logging or copying of stdout to file is done using the tee command.
	psub->logfile = malloc(sizeof(char) * 100);
	sprintf(psub->logfile, "log.txt");
	FILE *lfid = fopen(psub->logfile, "a");
	fprintf(lfid, "\n\n************************* %02d/%02d/%02d at %02d h %02d m %02d s *************************\n\n", tm.tm_mday, 1 + tm.tm_mon, tm.tm_year + 1900, tm.tm_hour, tm.tm_min, tm.tm_sec);
	fclose(lfid);

	psub->family = fam;
	psub->ncodes = lengths[0][0];
	(psub->lengths) = malloc(sizeof(int) * psub->ncodes);
	(psub->distances) = malloc(sizeof(int) * psub->ncodes);
	(psub->lattices) = malloc(sizeof(char **) * psub->ncodes);
	char *tile = malloc(sizeof(char) * 100);
	for (pi = 0; pi < psub->ncodes; pi ++){
		(psub->lengths)[pi] = lengths[1 + pi][0];
		(psub->distances)[pi] = lengths[1 + pi][2];
		TilingFile(tile, psub->family, (psub->lengths)[pi]);
		(psub->lattices)[pi] = malloc(sizeof(char *) * 2);
		(psub->lattices)[pi][0] = malloc(sizeof(char) * 100);
		sprintf((psub->lattices)[pi][0], "tilings/%s_edges.txt", tile);
		(psub->lattices)[pi][1] = malloc(sizeof(char) * 100);
		sprintf((psub->lattices)[pi][1], "tilings/%s_faces.txt", tile);
	}
	free(tile);

	// All the source files must be listed in a text file called "sources.txt" and thier count must be provided in the first line of this file.
	FILE *sfid = fopen("source/sources.txt", "r");
	fscanf(sfid, "%d", &(psub->nsource));
	(psub->source) = malloc(sizeof(char *) * (psub->nsource + 1));
	(psub->source)[0] = malloc(sizeof(char) * 100);
	sprintf((psub->source)[0], "Makefile");
	char *fname = malloc(sizeof(char) * 100);
	for (fi = 1; fi <= psub->nsource; fi ++){
		(psub->source)[fi] = malloc(sizeof(char) * 100);
		fscanf(sfid, "%s", fname);
		sprintf((psub->source)[fi], "source/%s", fname);
	}
	free(fname);
	fclose(sfid);

	// The tiling information filenaes are read for the code lengths that are part of the simulation
	psub->model = model;
	struct representation *repr = malloc(sizeof(struct representation));
	FixRepresentation(repr, psub->model, psub->family);
	psub->params = malloc(sizeof(double *) * (repr->nvars + 1));
	(psub->params)[0] = malloc(sizeof(double));
	(psub->params)[0][0] = ranges[0][0];
	(psub->params)[0][1] = ranges[0][1];
	psub->pnames = malloc(sizeof(char *) * repr->nvars);
	psub->njobs = psub->ncodes;
	for (pi = 1; pi <= repr->nvars; pi ++){
		(psub->pnames)[pi - 1] = malloc(sizeof(char) * 100);
		sprintf((psub->pnames)[pi - 1], "%s", (repr->varnames)[pi - 1]);
		(psub->params)[pi] = malloc(sizeof(double) * 3);
		(psub->params)[pi][0] = ranges[pi][0];
		(psub->params)[pi][1] = ranges[pi][1];
		(psub->params)[pi][2] = ranges[pi][2];
	}
	FreeRep(repr);

	int i, nrates, maxrange = 0, maxidx = 0;
	for (i = 1; i <= (psub->params)[0][0]; i ++){
		if ((i - 1) != ranges[0][1]){
			nrates = 1 + (int) (((psub->params)[i][1] - (psub->params)[i][0])/(psub->params)[i][2]);
			if (maxrange < nrates){
				maxrange = nrates;
				maxidx = i;
			}
		}
	}
	(psub->njobs) = (psub->njobs) * maxrange;

	// If the total number of jobs jobs is greater than 144, set the number of concurrent jobs to 144.
	// Or else set it to the total number of jobs available
	psub->concurrent = 144;
	if (psub->njobs < psub->concurrent)
		psub->concurrent = psub->njobs;

	// The job name, walltime and the queue are to be prompted from the user
	psub->job = malloc(sizeof(char) * 100);
	sprintf(psub->job, "%s", jobname);
	psub->queue = malloc(sizeof(char) * 100);
	sprintf(psub->queue, "%s", queue);
	psub->wall = wall;

	// The total time allowance on MS is 100000 cores-hours. On MP2 it is 700000 cores-hours.
	// The usage is calculated in percentage.
	if (strcmp(psub->host, "ms"))
		psub->usage = ((psub->njobs) * (psub->wall))/1000;
	else
		psub->usage = ((psub->njobs) * (psub->wall))/7000;
	return maxidx;
}


void FreeSubmission(struct submission *psub){
	// Free the submission structure
	int li;
	free(psub->instructions);
	free(psub->job);
	free(psub->host);
	free(psub->email);
	free(psub->queue);
	free(psub->timestamp);
	free(psub->pre);
	free(psub->input);
	free(psub->run);
	free(psub->post);
	free(psub->logfile);
	FreeCharArray(psub->pnames, (int) (psub->params)[0][0]);
	FreeDoubleArray(psub->params, 1 + (int) (psub->params)[0][0]);
	FreeCharArray(psub->source, psub->nsource);
	free(psub->lengths);
	for (li = 0; li < psub->ncodes; li ++)
		FreeCharArray((psub->lattices)[li], 2);
	free(psub->lattices);
}


void PrepareBQSubmitFile(struct submission *psub, int maxidx){
	// Prepare a bqsubmit.dat file with the details of the submission
	// Prepare an input file containing the parameter names, run file with the command to run and a post processing file called post.sh.
	int i;
	FILE *bqfid = fopen("bqsubmit.dat", "w");

	// batch name
	fprintf(bqfid, "batchName = %s\n\n", psub->job);

	// names of the files to be copies -- these are the source files alone
	fprintf(bqfid, "copyFiles = ");
	for (i = 0; i <= psub->nsource; i ++)
		fprintf(bqfid, "%s, ", (psub->source)[i]);
	fprintf(bqfid, "%s\n\n", psub->run);

	// names of the files to be linked -- these are the files that define the lattices used to define the codes which are part of the simulation
	// fprintf(bqfid, "linkFiles = %s, %s", (psub->lattices)[0][0], (psub->lattices)[0][1]);
	// for (fi = 1; fi < psub->ncodes; fi ++)
	// 	fprintf(bqfid, ", %s, %s", (psub->lattices)[fi][0], (psub->lattices)[fi][1]);
	// fprintf(bqfid, "\n\n");
	fprintf(bqfid, "linkFiles = tilings/\n\n");

	// name of the template file containing the input parameter names
	fprintf(bqfid, "templateFiles = %s\n\n", (psub->input));

	/*
	// These cooamnds are specific to mp2
	// 1. number of jobs to be run simulatenously on a node.
	// 2. Number of jobs to be accomulated for execution on a node.
	// Note that the second parameter must be a multiple of the first. For now, we will take both these as 1.
	// if (strncmp(psub->host, "mp2", 2))
	// 	fprintf(bqfid, "runJobsPerNode = 1\naccJobsPerNode = 1\n\n");
	*/

	// Preprossing
	fprintf(bqfid, "preBatch = ./%s\n\n", psub->pre);

	// Command to run on the compute node
	fprintf(bqfid, "command = ./%s\n\n", psub->run);

	// Post processing on the results file
	fprintf(bqfid, "postBatch = ./%s\n\n", psub->post);

	// Required resources for each node
	fprintf(bqfid, "submitOptions=-q %s@%s -l walltime=%d:00:00,nodes=1\n\n", psub->queue, psub->host, psub->wall);

	// Parameter ranges to explore in the simulation
	// The length and distance prameters of the code are explicitly written
	fprintf(bqfid, "param1 = length = [%d", (psub->lengths)[0]);
	for (i = 1; i < psub->ncodes; i ++)
		fprintf(bqfid, ", %d", (psub->lengths)[i]);
	fprintf(bqfid, "]\n");
	// For the other parameter that describe the noise model, it suffices to give a concise range
	fprintf(bqfid, "param2 = %s = %g:%g:%g\n", (psub->pnames)[maxidx - 1], (psub->params)[maxidx][0], (psub->params)[maxidx][2], (psub->params)[maxidx][1]);

	// Number of concurrent jobs that must be run at once. By default this is 144
	fprintf(bqfid, "concurrentJobs = %d\n\n", psub->concurrent);

	// Send notification to an email address
	fprintf(bqfid, "emailAddress = %s\n\n", psub->email);

	fclose(bqfid);
}

void PreProcessor(struct submission *psub){
	/*
		List of commands that must be executed before starting the job. These are executed right away after running bqsubmit.
		Copy the bqsubmit.dat file into a folder that is named using a code family and a timestamp
		These commands will be stored in a file called pre.sh.
		Set execute permissions to pre.sh
	*/
	FILE *prefp = fopen(psub->pre, "w");
	fprintf(prefp, "mkdir family%d/model%d/%s/\n", psub->family, psub->model, psub->timestamp);
	fprintf(prefp, "mkdir family%d/model%d/%s/results/\n", psub->family, psub->model, psub->timestamp);
	fprintf(prefp, "cp bqsubmit.dat %s %s family%d/model%d/%s/\n", psub->instructions, psub->input, psub->family, psub->model, psub->timestamp);
	fclose(prefp);
	chmod(psub->pre, S_IRWXU);
}

void Command(struct submission *psub){
	// Prepare a file containing the main command to be launched for the simulation, at every node.
	FILE *runfp = fopen(psub->run, "w");
	fprintf(runfp, "mkdir results/\n");
	fprintf(runfp, "mkdir source/\n");
	int i;
	for (i = 1; i <= psub->nsource; i ++)
		fprintf(runfp, "mv %s source/\n", (psub->source)[i] + 7);
	fprintf(runfp, "mkdir obj/\n");
	fprintf(runfp, "make clean | tee -a %s\n", psub->logfile);
	fprintf(runfp, "make | tee -a %s\n", psub->logfile);
	fprintf(runfp, "./benchmarking %s 2>&1 | tee -a %s\n", psub->input, psub->logfile);
	fprintf(runfp, "echo -e \"\\n ----------------------------------------------- \\n\" >> %s\n", psub->logfile);
	fclose(runfp);
	chmod(psub->run, S_IRWXU);
}

void PostProcessor(struct submission *psub){
	/*
		Prepare the contents of the post processing file
		In this file, we will collect all the results from the various .BQ folders and put them into a folder with the timestamp
		Finally, set permissions for this file to be executable.
	*/
	FILE *postfp = fopen(psub->post, "w");

	fprintf(postfp, "cp -r %s_*.BQ/results/* family%d/model%d/%s/results/\n", psub->job, psub->family, psub->model, psub->timestamp);

	fprintf(postfp, "cat %s_*.BQ/%s >> family%d/model%d/%s/SimulationDetails_%s.txt\n", psub->job, psub->logfile, psub->family, psub->model, psub->timestamp, psub->timestamp);

	struct representation *repr = malloc(sizeof(struct representation));
	FixRepresentation(repr, psub->model, psub->family);
	char *tileName = malloc(sizeof(char) * 100);
	int l;
	for (l = 0; l < psub->ncodes; l ++){
		TilingFile(tileName, psub->family, (psub->lengths)[l]);
		fprintf(postfp, "cat %s*length%d*.BQ/results/timestamp_%s_%s.txt >> family%d/model%d/%s/results/timestamp_%s_%s.txt\n", psub->job, (psub->lengths)[l], tileName, repr->mname, psub->family, psub->model, psub->timestamp, tileName, repr->mname);
	}
	free(tileName);
	FreeRep(repr);

	fprintf(postfp, "tar -zcvf family%d/model%d/%s.tar.gz family%d/model%d/%s\n", psub->family, psub->model, psub->timestamp, psub->family, psub->model, psub->timestamp);
	fclose(postfp);
	chmod(psub->post, S_IRWXU);
}


void InputFile(struct submission *psub, long stats, int maxidx){
	// Prepare the contents of the input file as well as a perf.txt file.
	// Return the index of the parameter whose range is the widest
	// maxidx has the index of the parameter whose values will be spread across nodes and other parameters will run locally in every node.
	FILE *inpfp = fopen(psub->input, "w");
	fprintf(inpfp, "mode simulate\n");
	fprintf(inpfp, "family %d\n", psub->family);
	fprintf(inpfp, "codes 1 ~~length~~\n");
	fprintf(inpfp, "noise %d\n", psub->model);
	int i;
	for (i = 1; i <= (psub->params)[0][0]; i ++){
		if (i != maxidx)
			fprintf(inpfp, "%s %g %g %g\n", (psub->pnames)[i - 1], (psub->params)[i][0], (psub->params)[i][1], (psub->params)[i][2]);
		else
			fprintf(inpfp, "%s ~~%s~~ ~~%s~~ 0.1\n", (psub->pnames)[i - 1], (psub->pnames)[i - 1], (psub->pnames)[i - 1]);
	}
	struct representation *repr = malloc(sizeof(struct representation));
	FixRepresentation(repr, psub->model, psub->family);
	if ((int) (psub->params)[0][1] != -1)
		fprintf(inpfp, "implicit %s\n", (repr->varnames)[(int) (psub->params)[0][1]]);
	FreeRep(repr);
	fprintf(inpfp, "trials %ld\n", stats);
	fclose(inpfp);
}

int BQCheck(struct submission *psub){
	/*
		Check if everything is OK with the bqsubmit file and all the included files
		1. Check if all the source and tiling files are present
		2. Check if run.sh, post.sh and input.txt are present.
		3. Check if job_*.BQ folders already exist, in which case bqsubmit will fail.
	*/
	int fi, pi, pass = 1, acc, overall = 1;
	for (fi = 0; fi <= psub->nsource; fi ++){
		acc = !(access((psub->source)[fi], R_OK) == -1);
		pass *= acc;
		if (acc == 0)
			fprintf(stderr, "\033[91m%s is missing.\033[0m\n", (psub->source)[fi]);
	}
	overall *= pass;
	if (pass == 0)
		fprintf(stderr, "\033[91mOne or more source files are missing.\033[0m\n");
	acc = 1;
	pass = 1;
	for (fi = 0; fi < psub->ncodes; fi ++){
		for (pi = 0; pi < 2; pi ++){
			acc = !(access((psub->lattices)[fi][pi], R_OK) == -1);
			pass *= acc;
			if (acc == 0)
				fprintf(stderr, "\033[91m%s is missing.\033[0m\n", (psub->lattices)[fi][pi]);
		}
	}
	overall *= pass;
	if (pass == 0)
		fprintf(stderr, "\033[91mOne or more lattice files are missing.\033[0m\n");
	if (overall == 1){
		printf("\033[92m_/ All looks good!\033[0m\n");
		printf("\033[92m_/ Created with the time-stamp %s.\033[0m\n", psub->timestamp);
		printf("\033[2m 	In all, %d nodes will be used for a maximum of %d hours.\n 	The entire wall-time will use up %g%% of the total quota on %s.\033[0m\n", psub->njobs, psub->wall, (double)(psub->usage), (psub->host));
	}

	if (psub->njobs > MAX_JOBS)
		printf("\033[2mThe total number of jobs is quiet high (= %d). Use 'microJobs = ...' to group jobs on a single node.\033[0m\n", psub->njobs);

	return overall;
}

void CopyLines(char *fromfname, FILE *tofp){
	// Copy contents from a file and append into another file
	const int LINE_WIDTH = 100;
	char line[LINE_WIDTH];
	FILE *infp = fopen(fromfname, "r");
	while (fgets(line, sizeof(line), infp))
		fprintf(tofp, "%s", line);
	fclose(infp);
}

void LogSubmission(struct submission *psub){
	// Append the parameters of this submission to the log file
	FILE *logfp = fopen(psub->logfile, "a");
	fprintf(logfp, "Submission was created with the time stamp: %s.\n", psub->timestamp);
	fprintf(logfp, "\n#### Input file: %s #### \n", psub->instructions);
	CopyLines(psub->instructions, logfp);
	fprintf(logfp, "\n#### bqsubmit file: bqsubmit.dat #### \n");
	CopyLines("bqsubmit.dat", logfp);
	fclose(logfp);
}

void Submit(char *instructions, int family, int **lengths, int model, double **params, long stats, char *jobname, char *queue, int wall){
	// Make a submission
	struct submission *pbs = malloc(sizeof(struct submission));
	int maxidx = PrepareSubmission(pbs, family, lengths, model, params, jobname, queue, wall, instructions);
	PreProcessor(pbs);
	InputFile(pbs, stats, maxidx);
	Command(pbs);
	PostProcessor(pbs);
	PrepareBQSubmitFile(pbs, maxidx);
	BQCheck(pbs);
	LogSubmission(pbs);
	FreeSubmission(pbs);
}
