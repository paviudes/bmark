#ifndef DRAW_H
#define DRAW_H

//create a pdf figure of the tiling
//c is 'b' for blue or 'r' for red
extern void DrawMin(struct tiling G, char *pdfName, char c);

//create a pdf figure of the tiling with indices
//pdfName is the name of the output pdf file
//c = 'b' for blue or 'r' for red
extern void Draw(struct tiling G, char *pdfName, char c);

//create a pdf figure of the tiling with indices
//c is 'b' for blue or 'r' for red
//color in red a subset of edges
//the second parameter is the indicator vector of this subset of egdes
extern void HighlightEdges(struct tiling, int *G, char c);

//Write the latex file and plot the performance curves
//the number of elements of the code family is familySize
//errorType = 1 for X and Z error, 2 for Z error, 3 for X error
extern void PlotPerformance(struct storePerformance *perf, int familySize, int errorType);

//generate a pdf with
//(i) a figure of the tiling and its dual
//(ii) the code parameters and overhead
//(iii) the weights of the measurements
//(iv) the performance plotted for X and Z errors, X-errors and Z-errors
extern void Report(struct tiling G, struct tiling H, struct storePerformance *perf);

//
extern void ComparisonReport(struct storePerformance *perf, int familySize);

#endif	/* DRAW_H */
