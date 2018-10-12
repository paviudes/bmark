#ifndef PLOT_H
#define PLOT_H

struct axes
{
	/*
		Specify the elements of all the axes of a plot
	*/

	// Column of the data that describes the values on the X-axis
	int xcol;

	// Column of the data that describes the values on the Y-axis
	int ycol;

	// Column of the data that describes the values on the Z-axis
	int zcol;

	// Column of the data that describes error bars
	int err;

	// Label for the x-axis
	char *xlabel;

	// Label for the y-axis
	char *ylabel;

	// Label for the z-axis
	char *zlabel;

	/*
		Scale used in the plot. The scale is specified by the corresponding tikz terminology. They are of three types.
		axis -- linear scale where both X and Y axis values are mentioned as such.
		semilogxaxis -- Only X axis is presented in the log-scale, while Y-axis is in linear scale.
		semilogyaxis -- Only Y axis is presented in the log-scale, while X-axis is in linear scale.
		loglogaxis -- Both X and Y are presented in the log-scale.
	*/
	char *scale;

	// Limits for the xaxis on the plot
	float *xlim;

	// Limits for the y axis on the plot
	float *ylim;
};


struct plot
{
	/*
		Specify the elements of the plot.
		A plot comprises of a given number of figures each of which have a given number of curves.
	*/

	// Number of figures in the plot
	int nfigs;

	// An ID assigned to the plot which is used to name the various files that are associated to the plot.
	char *ID;

	/*
		Data files for the plot.
		data[j][k] = name of the data file that is used to plot the k-th curve in the j-th figure.
	*/

	char **data;
	/*
		Titles for the figures in the plot
		titles[j] = title for the j-th figure in the plot.
	*/
	char **titles;

	// Global title for the entire set of plots
	char *title;

	/*
		Labels for the curves in each figure.
		labels[j] = label for the j-th curve in every figure.
	*/
	char **labels;

	// created date and time
	char *timestamp;
};

// Initialize an axes data structure
extern void InitAxes(struct axes *ax);

// Free an axes structure
extern void FreeAxes(struct axes *ax);

// Initialize a plot data structure
extern void InitPlot(struct plot *plt, int nfigs);

// Free a plot structure
extern void FreePlot(struct plot* plt);

/*
	Gather data from local files that contain performance data for a single noise rate, and copy them onto a global performance file.
		Input: family index, noise model, simulation time stamp, plotting time stamp.
		Output: All combinations of fixed parameters.

		Performs:
		The global performance files are produced for fixed values of a set of parameters. These parameter indices are to be specified in "fixed".
		1. Compute all combinations of fixed parameters and all combinations of variable parameters.
			A. Expand the fixed parameter and variable parameter ranges
			B. Compute the cartesian products of the expanded parameter ranges.
		2. For every combination of the fixed parameters, open a global file.
			A. Combine this with every combination of the variable parameters, open the resulting local data file and copy the contents onto the global file.
*/
extern void Gather(int family, int **lengths, int model, char *simultime, char *current);

/*
	Arrange the datafiles for performance plots. There will one plot for each combination of the fixed parameters.
	Each plot will contain a curve for every size of the tiling.
*/
extern void PerformancePlot(int family, int model, int **lengths, struct axes *ax, char *timestamp);

/*
	Produce runtime plots for the code family and the error model.
	The runtime for a given code length is just the worst case run time over all the possible noise parameter-values.
*/
extern void RuntimePlot(int family, int model, char *timestamp);

// Plot using latex tikz package
extern void Plot(struct plot *plt, struct axes *ax);

#endif	/* PLOT_H */
