#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include "memory.h"
#include "tiling.h"
#include "arrays.h"
#include "notations.h"
#include "save.h"
#include "plot.h"

void InitAxes(struct axes *ax){
	// Initialize the axes data structure
	ax->xlabel = malloc(sizeof(char) * 100);
	ax->ylabel = malloc(sizeof(char) * 100);
	ax->zlabel = malloc(sizeof(char) * 100);
	ax->scale = malloc(sizeof(char) * 100);
	ax->xlim = malloc(sizeof(float) * 2);
	(ax->xlim)[0] = -1000;
	(ax->xlim)[1] = 0;
	ax->ylim = malloc(sizeof(float) * 2);
	(ax->ylim)[0] = 0;
	(ax->ylim)[1] = 0;
	ax->err = -1;
}

void FreeAxes(struct axes *ax){
	// Free the axes data structure
	free(ax->xlabel);
	free(ax->ylabel);
	free(ax->zlabel);
	free(ax->scale);
	free(ax->xlim);
	free(ax->ylim);
}

void InitPlot(struct plot *plt, int nfigs){
	// Initialize the plot data structure
	plt->nfigs = nfigs;
	plt->ID = malloc(sizeof(char) * 100);
	plt->timestamp = malloc(sizeof(char) * 100);
	plt->title = malloc(sizeof(char) * 100);
	plt->titles = malloc(sizeof(char *) * nfigs);
	plt->data = malloc(sizeof(char *) * nfigs);
	plt->labels = malloc(sizeof(char *) * nfigs);
	int i;
	for (i = 0; i < nfigs; i ++){
		(plt->titles)[i] = malloc(sizeof(char) * 100);
		(plt->data)[i] = malloc(sizeof(char) * 100);
		(plt->labels)[i] = malloc(sizeof(char) * 100);
	}
}

void FreePlot(struct plot *plt){
	// Free the plot data structure
	free(plt->ID);
	free(plt->timestamp);
	free(plt->title);
	FreeCharArray(plt->titles, plt->nfigs);
	FreeCharArray(plt->data, plt->nfigs);
	FreeCharArray(plt->labels, plt->nfigs);
	free(plt);
}

void RunPython(char *fname, int verbose){
	// run a python file. Hide the outputs of the python program if verbose is 0.
	char* systemCommand = malloc(sizeof(char) * 300);
	if (verbose == 0)
		sprintf(systemCommand, "python %s>/dev/null 2>&1", fname);
	else
		sprintf(systemCommand, "python %s", fname);
	printf("\033[2m[Python] python %s\033[0m\n", fname);
	system(systemCommand);
}

void CompileTexFile(char* texFName, struct plot *plt){
	// Compile a tex file using pdflatex
	char* systemCommand = malloc(sizeof(char)*300);
	sprintf(systemCommand, "/Library/TeX/texbin/pdflatex -aux-directory=reports/ -output-directory=reports/ %s>/dev/null 2>&1", texFName);
	system(systemCommand);
	// Delete the auxiliary files, but keep the tex file.
	sprintf(systemCommand, "rm -rf reports/%s{.aux,.synctex.gz,.log}", plt->ID);
	system(systemCommand);
	free(systemCommand);
}

char* TeXFriendly(char *expression){
	// convert a string expression to TeX friendly format, by replacing the "_" with "\_".
	char *texExp = malloc(sizeof(char) * 100);
	long ei, ti = 0;
	for (ei = 0; ei < strlen(expression); ei ++){
		if (expression[ei] == '_')
			texExp[ti ++] = '\\';
		texExp[ti ++] = expression[ei];
	}
	texExp[ti] = 0;
	return texExp;
}

void Plot2D(struct plot *plt, struct axes *ax){
	// produce a 2D contour plot using numpy in Python.
	// For each data file in plt->data, do
	// read the plot file as numpy array
	// identify the x, y and the z axis
	// use http://stackoverflow.com/questions/9008370/python-2d-contour-plot-from-3-lists-x-y-and-rho
	// http://stackoverflow.com/questions/18892973/python-2d-contour-plot-from-3-lists-axes-not-generated-in-plot
	char *pythonFname = malloc(sizeof(char) * 300);
	sprintf(pythonFname, "reports/%s.py", plt->ID);
	char *plotfname = malloc(sizeof(char) * 300);
	sprintf(plotfname, "reports/%s.pdf", plt->ID);
	FILE *pyfp = fopen(pythonFname, "w");
	// preamble for python
	fprintf(pyfp, "############\n");
	fprintf(pyfp, "# This file was produced with the time stamp: %s.\n", plt->timestamp);
	fprintf(pyfp, "############\n\n");

	fprintf(pyfp, "try:\n");
	fprintf(pyfp, "\timport numpy as np\n");
	fprintf(pyfp, "\timport scipy.interpolate\n");
	fprintf(pyfp, "except Exception:\n");
	fprintf(pyfp, "\tsys.stderr.write(\"\\033[91mNumpy and/or Scipy do not exist, cannot do contour plots.\\n\\033[0m\")\n\n");

	fprintf(pyfp, "try:\n");
	fprintf(pyfp, "\timport matplotlib\n");
	fprintf(pyfp, "\tmatplotlib.use(\"Agg\")\n");
	fprintf(pyfp, "\tfrom matplotlib.backends.backend_pdf import PdfPages\n");
	fprintf(pyfp, "\timport matplotlib.pyplot as plt\n");
	fprintf(pyfp, "except Exception:\n");
	fprintf(pyfp, "\tsys.stderr.write(\"\\033[91mMatplotlib does not exist, cannot make plots.\\n\\033[0m\")\n\n");

	fprintf(pyfp, "try:\n");
	fprintf(pyfp, "\timport datetime as dt\n");
	fprintf(pyfp, "except Exception:\n");
	fprintf(pyfp, "\tsys.stderr.write(\"\\033[91mDatetime does not exist, cannot make timestamps.\\n\\033[0m\")\n\n");

	// fprintf(pyfp, "print(\"\\033[2m[Python] This file was produced with the time stamp: %s.\\033[0m\")\n", plt->timestamp);
	fprintf(pyfp, "npoints = 100\n");
	fprintf(pyfp, "(xcol, ycol, zcol) = (%d, %d, %d)\n\n", ax->xcol, ax->ycol, ax->zcol);
	fprintf(pyfp, "with PdfPages(\"%s\") as pdf:\n", plotfname);
	int i;
	for (i = 0; i < plt->nfigs; i ++){
		fprintf(pyfp, "\n\t######################\n");
		fprintf(pyfp, "\t# Figure %d\n", i + 1);
		fprintf(pyfp, "\t# Read the data in %s into a numpy array. Remove lines starting with \"#\"\n\n", plt->data[i]);

		fprintf(pyfp, "\tdataset = np.loadtxt(\"%s\", comments = \"#\")\n", plt->data[i]);
		fprintf(pyfp, "\trelevant = np.array([(dataset[ri, xcol], dataset[ri, ycol]) for ri in range(dataset.shape[0])])\n");
		fprintf(pyfp, "\t(__, uniqueIndices) = np.unique([\"%%g%%g\" %% (dataset[ri, xcol], dataset[ri, ycol]) for ri in range(dataset.shape[0])], return_index = True)\n");
		
		// set the x and y axis limits
		if ((ax->xlim)[0] != -1000){
			fprintf(pyfp, "\t# set the x and y axis limits\n");
			fprintf(pyfp, "\tuniqueIndices = uniqueIndices[np.where(dataset[uniqueIndices, xcol] >= %g)]\n", (ax->xlim)[0]);
			fprintf(pyfp, "\tuniqueIndices = uniqueIndices[np.where(dataset[uniqueIndices, xcol] <= %g)]\n", (ax->xlim)[1]);
			fprintf(pyfp, "\tuniqueIndices = uniqueIndices[np.where(dataset[uniqueIndices, ycol] >= %g)]\n", (ax->ylim)[0]);
			fprintf(pyfp, "\tuniqueIndices = uniqueIndices[np.where(dataset[uniqueIndices, ycol] <= %g)]\n", (ax->ylim)[1]);
		}

		if (strncmp(ax->scale, "log", 3) == 0){
			fprintf(pyfp, "\t# In order to implement log scale, we must remove zeros from the z-dataset.\n\t# Since the zeros are just empirical zeros, they can be set to the lowest finite number.\n");
			fprintf(pyfp, "\tnonzero = uniqueIndices[np.where(dataset[uniqueIndices, zcol] > 0)]\n");
			fprintf(pyfp, "\tsmallestNonzero = np.min(dataset[nonzero, zcol])\n");
			fprintf(pyfp, "\tzero = uniqueIndices[np.where(dataset[uniqueIndices, zcol] == 0)]\n");
			fprintf(pyfp, "\tdataset[zero, zcol] = smallestNonzero\n");
			fprintf(pyfp, "\tdataset[uniqueIndices, zcol] = -np.log10(dataset[uniqueIndices, zcol])\n");
		}
		fprintf(pyfp, "\tprint(\"\\033[2m[Python] The data set has %%d unique points.\\033[0m\" %% uniqueIndices.shape[0])\n");

		// Setup a regular grid of interpolated points in the input space
		fprintf(pyfp, "\tdatax, datay = np.linspace(np.min(dataset[uniqueIndices, xcol]), np.max(dataset[uniqueIndices, xcol]), npoints), np.linspace(np.min(dataset[uniqueIndices, ycol]), np.max(dataset[uniqueIndices, ycol]), npoints)\n");
		fprintf(pyfp, "\tmeshx, meshy = np.meshgrid(datax, datay)\n");

		// Interpolation
		fprintf(pyfp, "\tmeshz = scipy.interpolate.griddata((dataset[uniqueIndices, xcol], dataset[uniqueIndices, ycol]), dataset[uniqueIndices, zcol], (meshx, meshy), method='linear')\n");

		fprintf(pyfp, "\tfig = plt.figure()\n");
		fprintf(pyfp, "\tax = plt.gca()\n");
		fprintf(pyfp, "\tplt.title(\"%s\")\n", plt->titles[i]);
		fprintf(pyfp, "\tax.set_xlabel(\"%s\")\n", ax->xlabel);
		fprintf(pyfp, "\tax.set_ylabel(\"%s\")\n", ax->ylabel);
		fprintf(pyfp, "\tplt.imshow(meshz, vmin = np.min(dataset[uniqueIndices, zcol]), vmax = np.max(dataset[uniqueIndices, zcol]), extent = [np.min(dataset[uniqueIndices, xcol]), np.max(dataset[uniqueIndices, xcol]), np.min(dataset[uniqueIndices, ycol]), np.max(dataset[uniqueIndices, ycol])], origin = 'lower', aspect='auto')\n");
		fprintf(pyfp, "\tplt.scatter(dataset[uniqueIndices, xcol], dataset[uniqueIndices, ycol], c = dataset[uniqueIndices, zcol])\n");
		fprintf(pyfp, "\tplt.colorbar(ticks = np.linspace(np.min(dataset[uniqueIndices, zcol]), np.max(dataset[uniqueIndices, zcol]), 10))\n");

		// setting the legend and saving the plot
		fprintf(pyfp, "\tprint(\"\\033[2m[Python] Saving plot for figure %d\\033[0m\")\n", i + 1);
		fprintf(pyfp, "\t# save the current figure into a pdf page\n");
		fprintf(pyfp, "\tpdf.savefig(fig)\n");
		fprintf(pyfp, "\tplt.close()\n");
	}

	// Set the PDF attributes
	fprintf(pyfp, "\n\n\t######################\n# Set the PDF attributes\n");
	fprintf(pyfp, "\tpdfInfo = pdf.infodict()\n");
	fprintf(pyfp, "\tpdfInfo['Title'] = \"%s\"\n", plt->title);
	fprintf(pyfp, "\tpdfInfo['Author'] = \"Pavithran Iyer\"\n");
	fprintf(pyfp, "\tpdfInfo['ModDate'] = dt.datetime.today()\n");

	fprintf(pyfp, "\nprint(\"\\033[2m[Python] All plots saved to %s.\\033[0m\")\n", plotfname);

	fclose(pyfp);
	free(plotfname);

	// Run the python file
	RunPython(pythonFname, 1);
	free(pythonFname);
}


void Plot(struct plot *plt, struct axes *ax){
	// Plot using latex tikz package
	// Plot filename is the family index followed by noise model
	char* texFName = malloc(sizeof(char) * 300);
	sprintf(texFName, "reports/%s.tex", plt->ID);
	FILE* texFid = fopen(texFName, "w");
	// Tikz document preambles
	fprintf(texFid, "\\documentclass{article}\n");
	fprintf(texFid, "\\usepackage[left = 1cm, right = 2cm, top = 2cm, bottom = 2cm, landscape]{geometry}\n");
	fprintf(texFid, "\\pagestyle{empty}\n");
	fprintf(texFid, "\\usepackage{pgfplots}\n");
	fprintf(texFid, "\\pgfplotsset{scale=2.5, compat = newest, every tick/.append style = thin}\n");
	fprintf(texFid, "\\begin{document}\n");

	fprintf(texFid, "\\begin{center}\n");
	fprintf(texFid, "\\begin{tikzpicture}\n");
	fprintf(texFid, "\\begin{%s}[xlabel={%s}, ylabel={%s}, grid=major, x tick label style={rotate=90, /pgf/number format/fixed, /pgf/number format/precision=3}, title = {%s}]\n", ax->scale, ax->xlabel, ax->ylabel, TeXFriendly(plt->title));
	int i;
	for(i = 0; i < plt->nfigs; i ++){
		if(ax->err == -1)
			fprintf(texFid, "\\addplot table[x index=%d,y index=%d, col sep=space] {%s};\n", ax->xcol, ax->ycol, (plt->data)[i]);
		else
			fprintf(texFid, "\\addplot+[error bars/.cd, y dir=both,y explicit] table[x index=%d,y index=%d, y error index = %d, col sep=space] {%s};\n", ax->xcol, ax->ycol, ax->err, (plt->data)[i]);
		if (plt->nfigs > 1)
			fprintf(texFid, "\\addlegendentry{%s}\n", (plt->labels)[i]);
	}
	fprintf(texFid, "\\end{%s}\n", ax->scale);
	fprintf(texFid, "\\end{tikzpicture}\n");
	fprintf(texFid, "\\end{center}\n");

	fprintf(texFid, "\\end{document}\n");
	fclose(texFid);
	// Show the plot by compiling the tex file using pdflatex
	CompileTexFile(texFName, plt);
	free(texFName);
}

void CopyAndTime(char *fromfile, FILE* tofp, int nparams, double *runInfo){
	// Copy all the performance data from a local file to a global file.
	// The performance data can at most have 1000 characters.
	// Return the runtime data.
	// fprintf(stderr, "Reading %s\n", fromfile);

	int i;
	double *worstRun = malloc(sizeof(double) * 2), *physical = malloc(sizeof(double) * nparams), *failrates = malloc(sizeof(double) * 3), *ebrs = malloc(sizeof(double) * 3);
	long *totalTrials = malloc(sizeof(long) * 3), *nfails = malloc(sizeof(long) * 3);
	if (access(fromfile, R_OK) == -1){
		fprintf(stderr, "\033[2m!! Missing file %s.\033[0m\n", fromfile);
		fprintf(tofp, "0");
		for (i = 1; i < nparams + 12; i ++)
			fprintf(tofp, " 0");
		fprintf(tofp, "\n");
	}
	else{
		FILE *locfp = fopen(fromfile, "r");
		while(feof(locfp) == 0){
			for (i = 0; i < nparams; i ++)
				fscanf(locfp, "%lf", &(physical[i]));
			for(i = 0; i < 3; i ++)
				fscanf(locfp, "%lf %ld %ld %lf", &(failrates[i]), &(nfails[i]), &(totalTrials[i]), &(ebrs[i]));
			fscanf(locfp, "%lf %lf", &(worstRun[0]), &(worstRun[1]));
			// Write into another file
			fprintf(tofp, "%lf", physical[0]);
			for (i = 1; i < nparams; i ++)
				fprintf(tofp, " %lf", physical[i]);
			for (i = 0; i < 3; i ++)
				fprintf(tofp, " %.4e %ld %ld %.4e", failrates[i], nfails[i], totalTrials[i], ebrs[i]);
			fprintf(tofp, "\n");
			runInfo[0] = worstRun[0];
			runInfo[1] += worstRun[1];
		}
		fclose(locfp);
	}
	free(worstRun);
	free(physical);
	free(failrates);
	free(totalTrials);
	free(nfails);
	free(ebrs);
}


void Gather(int family, int **lengths, int model, char *simultime, char *current){
	// Read the timestamps.txt file. It contains the timestamp of the various simulations done for different ranges of noise rates.
	// Read every simulation result file and append the results to a global file, with another timestamp which will be used by a plot software.
	int l, p, ct;
	char *tileName = malloc(sizeof(char) * 100);
	char *allstamps = malloc(sizeof(char) * 100);
	char *stamp = malloc(sizeof(char) * 100);
	char *localresults = malloc(sizeof(char) * 200);
	char *globalresults = malloc(sizeof(char) * 200);
	double *runInfo = malloc(sizeof(double) * 3);
	char *runfile = malloc(sizeof(char) * 300);
	struct representation *repr = malloc(sizeof(struct representation));
	FixRepresentation(repr, model, family);

	FILE *allstampsfp, *globalfp, *runfp;

	sprintf(runfile, "results/runtimes_%s_%s_%s.txt", repr->fname, repr->mname, current);
	runfp = fopen(runfile, "w");
	for (l = 1; l <= lengths[0][0]; l ++){
		FillDoubleArrayWithConstants(runInfo, 3, 0);
		TilingFile(tileName, family, lengths[l][0]);
		sprintf(allstamps, "./../family%d/model%d/%s/results/timestamp_%s_%s.txt", family, model, simultime, tileName, repr->mname);
		sprintf(globalresults, "results/perf_%s_%s_%s.txt", tileName, repr->mname, current);
		globalfp = fopen(globalresults, "w");
		fprintf(globalfp, "# File to be used by some plot software.\n# Simulation done on %s.\n# %s", simultime, (repr->varnames)[0]);
		for (p = 1; p < repr->nvars; p ++)
			fprintf(globalfp, "\t%s", (repr->varnames)[p]);
		for (p = 0; p < repr->nprops; p ++)
			fprintf(globalfp, "\t%s", (repr->propnames)[p]);
		for (ct = 0; ct < 3; ct ++)
			fprintf(globalfp, "\tf_%d\tnf_%d\ttt_%d\te_%d", ct, ct, ct, ct);
		fprintf(globalfp, "\n");

		allstampsfp = fopen(allstamps, "r");
		while(feof(allstampsfp) == 0){
			fscanf(allstampsfp, "%s", stamp);
			sprintf(localresults, "./../family%d/model%d/%s/results/perf_%s_%s_%s.txt", family, model, simultime, tileName, repr->mname, stamp);
			CopyAndTime(localresults, globalfp, repr->nvars + repr->nprops, runInfo);
		}
		fclose(globalfp);
		fclose(allstampsfp);
		fprintf(runfp, "%d %g\n", (int) runInfo[0], (double) (runInfo[1]/(double) 3600));
	}
	free(tileName);
	free(allstamps);
	free(localresults);
	free(globalresults);
	FreeRep(repr);
	fclose(runfp);
	free(runfile);
	free(runInfo);
}

void PerformancePlot(int family, int model, int **lengths, struct axes *ax, char *timestamp){
	// Plot performance curves (surfaces)
	char *tileName = malloc(sizeof(char) * 100);
	struct representation *repr = malloc(sizeof(struct representation));
	FixRepresentation(repr, model, family);
	struct plot *plt = malloc(sizeof(struct plot));
	InitPlot(plt, lengths[0][0]);
	sprintf(plt->ID, "perf_%s_%s_%s", repr->fname, repr->mname, timestamp);
	sprintf(plt->timestamp, "%s", timestamp);
	sprintf(plt->title, "Performance of %s family with %s type noise.", repr->fname, repr->mname);
	int i;
	for (i = 0; i < plt->nfigs; i ++){
		sprintf((plt->titles)[i], "length %d %s code", lengths[i + 1][0], repr->fname);
		TilingFile(tileName, family, lengths[i + 1][0]);
		sprintf((plt->data)[i], "results/perf_%s_%s_%s.txt", tileName, repr->mname, timestamp);
	}

	Plot2D(plt, ax);

	free(tileName);
	FreeRep(repr);
	FreePlot(plt);
}

void RuntimePlot(int family, int model, char *timestamp){
	// Produce runtime plots for the code family and the error model.
	// The runtime for a given code length is just the worst case run time over all the possible noise parameter-values.
	struct representation *repr = malloc(sizeof(struct representation));
	FixRepresentation(repr, model, family);
	struct axes *ax = malloc(sizeof(struct axes));
	InitAxes(ax);
	struct plot *plt = malloc(sizeof(struct plot));
	InitPlot(plt, 1);
	sprintf(plt->timestamp, "%s", timestamp);
	sprintf(plt->ID, "runtimes_%s_%s_%s", repr->fname, repr->mname, timestamp);
	sprintf((plt->data)[0], "./../results/runtimes_%s_%s_%s.txt", repr->fname, repr->mname, timestamp);
	sprintf(plt->title, "Performance of %s family with %s type noise.", repr->fname, repr->mname);

	ax->xcol = 0;
	sprintf(ax->xlabel, "Number of qubits");
	ax->ycol = 1;
	sprintf(ax->ylabel, "Runtime in hours");
	sprintf(ax->scale, "axis");
	// there are no labels for the curves
	Plot(plt, ax);

	FreeRep(repr);
	FreePlot(plt);
	FreeAxes(ax);
}
