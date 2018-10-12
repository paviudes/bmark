#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "tiling.h"

//====================
// 	Latex header and footer
//	for a minimal tikz figure
//====================


//Tikz header to generated only the figure
void WriteTikzMinHeader(FILE *f, char *jobName, float scale)
{
	char *header = malloc(500*sizeof(char));
	sprintf(header, "\\documentclass{article}\n\
\\usepackage{graphics}\n\
\\usepackage{tikz}\n\
\\usepackage{pgfplots}\n\n\
\\pgfrealjobname{minimal}\n\n\
\\begin{document}\n\n\
\\beginpgfgraphicnamed{%s}\n\
\\begin{tikzpicture}[scale=%f, every node/.style={scale=%f}]\n\n", jobName, scale, scale*.25);

	fprintf(f, "%s", header);
//	free(header);
}


//Tikz footer to generated only the figure
void WriteTikzMinFooter(FILE *f)
{
	fprintf(f, "\\end{tikzpicture}\n\
\\endpgfgraphicnamed\n\n\
\\end{document}\n");
}


//====================
// 	Latex header and footer
//====================


//full article latex header
void WriteLatexHeader(FILE *f)
{
	fprintf(f, "\\documentclass[12pt]{article}\n\
\\usepackage{fullpage,amsmath,amssymb,amsfonts,amsthm}\n\
\\usepackage[utf8]{inputenc}\n\
\\usepackage[english]{babel}\n\
\\usepackage{graphics}\n\
\\usepackage{tikz}\n\
\\usepackage{pgfplots}\n\n\
\\pgfplotsset{compat=1.8}\n\
\\title{Performance Report}\n\n\
\\begin{document}\n\n\
\\maketitle\n\n");
}


void WriteLatexFooter(FILE *f)
{
	fprintf(f, "\\end{document}");
}


void BeginFigure(FILE *f)
{
	fprintf(f, "\\begin{figure}[h!]\n");
}


void EndFigure(FILE *f, char *caption)
{
	fprintf(f, "\\caption{%s}\n\\end{figure}\n", caption);
	fprintf(f, "\\vspace{.2cm}\n\n");
}


//Tikz header inside a latex document
void WriteTikzHeader(FILE *f, float figureScale, float nodeScale)
{
	char *header = malloc(100*sizeof(char));
	sprintf(header, "\\begin{tikzpicture}[scale=%.2f, every node/.style={scale=%.2f}]\n\n", figureScale, nodeScale);
	fprintf(f, "%s", header);
	free(header);
}


//Tikz footer inside a latex document
void WriteTikzFooter(FILE *f)
{
	fprintf(f, "\\end{tikzpicture}\n");
}


//====================================
//			Tikz vertices, edges and faces
//====================================


//define the tikz positions of vertices in the plane
void DefineLabels(struct tiling G, FILE *f)
{
	fprintf(f, "%%vertex labels\n");
	fprintf(f, "\\draw");
	int i;
	for (i=0; i<G.v; i++)
		fprintf(f, "\n	(%f,%f) node (v%d) {}", G.Coord[i][0], G.Coord[i][1], i);
	fprintf(f, ";\n\n");
}


//plot indexed vertices
void PlotIndexedVertices(struct tiling G, FILE *f)
{
	fprintf(f, "%%vertex indices\n");
	fprintf(f, "\\draw");
	int i;
	for (i=0; i<G.v; i++)
	{
		if (G.doV[i] == 0)
			fprintf(f, "\n	(%f,%f) node[draw, circle, inner sep=2pt, fill=black, label = above right:v%d] () {}", G.Coord[i][0], G.Coord[i][1], i);
		else
		fprintf(f, "\n	(%f,%f) node[draw, circle, inner sep=2pt, fill=white, label = above right:v%d] () {}", G.Coord[i][0], G.Coord[i][1], i);
	} 
	fprintf(f, ";\n\n");
}


//plot non-open vertices using their labels
void PlotNonOpenVertices(struct tiling G, FILE *f)
{	
	fprintf(f, "%%plot non-open vertices\n");
	fprintf(f, "\\draw");
	int i;
	for (i=0; i<G.v; i++)
		if (G.doV[i] == 0)
			fprintf(f, "\n	(v%d) node[draw, circle, inner sep=2pt, fill=black] {}", i);
	fprintf(f, ";\n\n");
}


//plot non-open vertices using their labels
void PlotOpenVertices(struct tiling G, FILE *f)
{	
	fprintf(f, "%%plot open vertices\n");
	fprintf(f, "\\draw");
	int i;
	for (i=0; i<G.v; i++)
		if (G.doV[i] == 1)
			fprintf(f, "\n	(v%d) node[draw, circle, inner sep=2pt, fill=white] {}", i);
	fprintf(f, ";\n\n");
}


//plot edges from the nodes labels of the vertices
void PlotEdges(struct tiling G, FILE *f)
{
	int i;
	//closed egdes
	fprintf(f, "%%plot non-open edges\n");	
	fprintf(f, "\\draw");
	for (i=0; i<G.e; i++)
		if (G.doE[i] == 0)
			fprintf(f, "\n	(v%d)--(v%d)", G.E[i][0], G.E[i][1]);
	fprintf(f, ";\n\n");

	//open edges
	fprintf(f, "%%plot open edges\n");
	fprintf(f, "\\draw[dashed]");
	for (i=0; i<G.e; i++)
		if (G.doE[i] == 1)
			fprintf(f, "\n	(v%d)--(v%d)", G.E[i][0], G.E[i][1]);
	fprintf(f, ";\n\n");	
}


//plot in red the edges given by their indicator vector subsetE
void PlotRedEdges(struct tiling G, FILE *f, int *subsetE)
{
	int i;
	//closed egdes
	fprintf(f, "%%plot non-open edges\n");
	fprintf(f, "\\draw[color=red, line width=1.5pt]");
	for (i=0; i<G.e; i++)
		if ((subsetE[i] == 1) && (G.doE[i] == 0))
			fprintf(f, "\n	(v%d)--(v%d)", G.E[i][0], G.E[i][1]);
	fprintf(f, ";\n\n");

	//open edges
	fprintf(f, "%%plot open edges\n");
	fprintf(f, "\\draw[dashed, color=red, line width=1.5pt]");
	for (i=0; i<G.e; i++)
		if ((subsetE[i] == 1) && (G.doE[i] == 1))
			fprintf(f, "\n	(v%d)--(v%d)", G.E[i][0], G.E[i][1]);
	fprintf(f, ";\n\n");
}



//plot faces from the nodes labels of the vertices
void PlotFaces(struct tiling G, FILE *f,char c)
{
	int *fiVertexList;
	int i, j;
	fprintf(f, "%%fill faces\n");
	if (c == 'r')
		fprintf(f, "\\fill[color=red!20, opacity=.5] ");
	else 
		fprintf(f, "\\fill[color=blue!20, opacity=.5] ");		
	for (i=0; i<G.f; i++)
	{
		fiVertexList = FaceToVertexList(G, G.F[i]);
		fprintf(f, "\n	");
		for (j=1; j<=fiVertexList[0]; j++)
			fprintf(f, "(v%d.center)--", fiVertexList[j]);
		fprintf(f, "cycle");
	}
	fprintf(f, ";\n\n");
}


//plot edges from the nodes labels of the vertices
void PlotIndexedEdges(struct tiling G, FILE *f)
{
	int i;
	//closed egdes
	fprintf(f, "%%plot non-open edges\n");	
	fprintf(f, "\\draw");
	for (i=0; i<G.e; i++)
		if (G.doE[i] == 0)
			fprintf(f, "\n	(v%d)--(v%d) node [midway] {e%d}", G.E[i][0], G.E[i][1], i);
	fprintf(f, ";\n\n");

	//open edges
	fprintf(f, "%%plot open edges\n");
	fprintf(f, "\\draw[dashed]");
	for (i=0; i<G.e; i++)
		if (G.doE[i] == 1)
			fprintf(f, "\n	(v%d)--(v%d) node [midway] {e%d}", G.E[i][0], G.E[i][1], i);
	fprintf(f, ";\n\n");	
}


//plot face indices
void PlotFaceIndices(struct tiling G, FILE *f)
{
	int i, j, ej;
	float centerx, centery, ux, uy, vx, vy;
	fprintf(f, "%%plot faces indices\n");
	fprintf(f, "\\draw");
	for (i=0; i<G.f; i++)
	{
		//compute the coordinate of the center fo the face
		centerx = 0;
		centery = 0;
		for (j=1; j<=G.F[i][0]; j++)
		{
			ej = G.F[i][j];
			ux = G.Coord[G.E[ej][0]][0];
			vx = G.Coord[G.E[ej][1]][0];
			uy = G.Coord[G.E[ej][0]][1];
			vy = G.Coord[G.E[ej][1]][1];
			centerx += ux + vx;
			centery += uy + vy;
		}
		centerx = centerx/(2*G.F[i][0]);
		centery = centery/(2*G.F[i][0]);
		fprintf(f, "\n   (%f,%f) node[color=blue!80] {f%d}", centerx, centery, i);
	}
	fprintf(f, ";\n\n");
}



//the scale is chosen to be 1/10 time the size of the figure
float FindScale(struct tiling G)
{
	float xmax = G.Coord[0][0];
	float xmin = G.Coord[0][0];
	float ymax = G.Coord[0][1];
	float ymin = G.Coord[0][1];
	int i;
	for (i=0; i<G.v; i++)
	{
		if (G.Coord[i][0] > xmax) xmax = G.Coord[i][0];
		if (G.Coord[i][0] < xmin) xmin = G.Coord[i][0];
		if (G.Coord[i][1] > ymax) ymax = G.Coord[i][1];
		if (G.Coord[i][1] < ymin) ymin = G.Coord[i][1];
	}
	
	if (xmax-xmin > ymax-ymin) return 20.0/(xmax-xmin);
	else return 20.0/(ymax-ymin);
}


//====================
// 		tikz data plot
//====================


//
void PlotAxis(FILE *f, struct storePerformance *perf, int familySize, char *dataFolder, int errorType)
{
	int i;
	
	//compute xmin and xmax
	double xmin = 1;
	double xmax = 0;
	for (i=0; i<familySize; i++)
	{
		if (perf[i].noiseRange[0] < xmin)
			xmin = perf[i].noiseRange[0];
		if (perf[i].noiseRange[1] > xmax)
			xmax = perf[i].noiseRange[1];
	}

	//compute the nb of non-empty data files
	int numberOfFiles = 0;
	for (i=0; i<familySize; i++)
		if (perf[i].dataFile != NULL)
			numberOfFiles++;

	fprintf(f, "\\begin{semilogyaxis}\n\
[\n\
    grid= major,\n\
    xlabel = {erasure probability},\n\
    ylabel = {logical error rate after decoding},\n");
	fprintf(f, "    xmin = %.2f, xmax = %.2f,\n", xmin, xmax);
	fprintf(f, "    ymax = 1,\n");
	fprintf(f, "    legend entries={");

	//plot index
	int fileIndex = 0;
	for (i=0; i<familySize; i++)
		if (perf[i].dataFile != NULL)
		{
			fprintf(f, "$%s$", perf[i].tilingName);
			fileIndex++;
			if (fileIndex == numberOfFiles)
				fprintf(f, "},\n");
			else
				fprintf(f, ", ");
		}

	fprintf(f, "    legend pos={south east},\n]\n\n");

	//plot data
	for (i=0; i<familySize; i++)
		if (perf[i].dataFile != NULL)
		{
			fprintf(f, "\\addplot+[\n");
			fprintf(f, "    error bars/.cd,\n");
			fprintf(f, "    y dir=both,\n");
			fprintf(f, "    y explicit,\n");
			fprintf(f, "    ]table[x index=0, y index=%d, y error index=%d]{%s%s};\n", errorType, errorType+7, dataFolder, perf[i].dataFile);
		}

	fprintf(f, "\n\\end{semilogyaxis}\n\n\n");
}


//====================================
//			table of code parameters
//====================================


void WriteCodeParameters(FILE *f, int n, int k)
{
	fprintf(f, "\\begin{table}[h]\n\
\\centering\n\
\\begin{tabular}{c c}\n");
	fprintf(f, "\\hline\n");
	fprintf(f, "Number of physical qubits & $n = %d$ \\\\\n", n);
	fprintf(f, "\\hline\n");
	fprintf(f, "Number of logical qubits & $k = %d$\\\\\n", k);
	fprintf(f, "\\hline\n");
	fprintf(f, "Overhead & $n/k = %.2f$\\\\\n", ((float)n)/((float) k));
	fprintf(f, "\\hline\n");
	fprintf(f, "\\end{tabular}\n\
\\caption{Error-correction overhead}\n\
\\end{table}\n\
\\vspace{.3cm}\n\n\n");
}


void WriteComparisonCodeParameters(FILE *f, struct storePerformance *perf, int familySize)
{
	int i;
	fprintf(f, "\\begin{table}[h]\n\
\\centering\n\
\\begin{tabular}{c c c c}\n");
	fprintf(f, "\\hline\n");
	fprintf(f, "name & Number of physical qubits & Number of logical qubits & Overhead\\\\\n");
	fprintf(f, "\\hline\n");
	for (i=0; i<familySize; i++)
	{
		fprintf(f, "$%s$ & $n = %d$ & $k = %d$ & $n/k = %.2f$\\\\\n", perf[i].tilingName, perf[i].n, perf[i].k, ((float) perf[i].n)/((float) perf[i].k));
		fprintf(f, "\\hline\n");
	}
	fprintf(f, "\\end{tabular}\n\
\\caption{Error-correction overhead}\n\
\\end{table}\n\
\\vspace{.3cm}\n\n\n");
}


//====================================
//			table of measurement weights
//====================================


void WriteWeights(FILE *f, int *ZWeight, int *XWeight, int maxWeight)
{
	fprintf(f, "\\begin{table}[h]\n\
\\centering\n\
\\begin{tabular}{c c}\n");
	fprintf(f, "\\hline\n");
	fprintf(f, "weight & number of $Z$-measurements\\\\\n");
	fprintf(f, "\\hline\n");
	int i;
	for (i=0; i<=maxWeight; i++)
		if (ZWeight[i] != 0)
			fprintf(f, "%d & %d\\\\\n", i, ZWeight[i]);
	fprintf(f, "\\hline\n");
	fprintf(f, "\\hline\n");
	fprintf(f, "weight & number of $X$-measurements\\\\\n");
	fprintf(f, "\\hline\n");
	for (i=0; i<=maxWeight; i++)
		if (XWeight[i] != 0)
			fprintf(f, "%d & %d\\\\\n", i, XWeight[i]);
	fprintf(f, "\\hline\n");
	fprintf(f, "\\end{tabular}\n\
\\caption{Measurements weight distribution}\n\
\\end{table}\n\
\\vspace{.3cm}\n\n\n");
}


//
void WriteComparisonWeights(FILE *f, struct storePerformance *perf, int familySize)
{
	int i, j;

	//compute the max weight over all codes
	int maxWeight = 0;
	for (i=0; i<familySize; i++)
		if (perf[i].maxWeight > maxWeight)
			maxWeight = perf[i].maxWeight;

	//indicator vector of weight giving non-trivial column in the table
	int *TrivialColumnZ = malloc((maxWeight+1)*sizeof(int));
	for (i=0; i<=maxWeight; i++)
	{
		TrivialColumnZ[i] = 1;
		for (j=0; j<familySize; j++)
			if (i <= perf[j].maxWeight)
				if (perf[j].ZWeight[i] != 0)
					TrivialColumnZ[i] = 0;
	}

	//indicator vector of weight giving non-trivial column in the table
	int *TrivialColumnX = malloc((maxWeight+1)*sizeof(int));
	for (i=0; i<=maxWeight; i++)
	{
		TrivialColumnX[i] = 1;
		for (j=0; j<familySize; j++)
			if (i <= perf[j].maxWeight)
				if (perf[j].XWeight[i] != 0)
					TrivialColumnX[i] = 0;
	}


	//begin the table
	fprintf(f, "\\begin{table}[h]\n\
\\centering\n\
\\begin{tabular}{ c ");
	for (i=0; i<=maxWeight; i++)
		if ((!TrivialColumnZ[i]) || (!TrivialColumnX[i]))
			fprintf(f, "c ");
	fprintf(f, "}\n");
	fprintf(f, "\\hline\n");

	//write first line with the list of Z-weights
	fprintf(f, "names $\\backslash$ $Z$-weights ");
	for (i=0; i<=maxWeight; i++)
		if ((!TrivialColumnZ[i]) || (!TrivialColumnX[i]))
			fprintf(f, "& %d ", i);
	fprintf(f, "\\\\\n");

	//add a Z-line for each code of the family
	for (i=0; i<familySize; i++)
	{
		fprintf(f, "\\hline\n");
		fprintf(f, "$%s$ ",perf[i].tilingName);
		for (j=0; j<=maxWeight; j++)
			if ((!TrivialColumnZ[j]) || (!TrivialColumnX[j]))
			{
				if (i <= perf[j].maxWeight)
					fprintf(f, "& %d ", perf[i].ZWeight[j]);
				else
					fprintf(f, "& 0 ");
			}

		fprintf(f, "\\\\\n");
	}
	fprintf(f, "\\hline\n");
	fprintf(f, "\\hline\n");

	//write first line with the list of X-weights
	fprintf(f, "names $\\backslash$ $X$-weights ");
	for (i=0; i<=maxWeight; i++)
		if ((!TrivialColumnZ[i]) || (!TrivialColumnX[i]))
			fprintf(f, "& %d ", i);
	fprintf(f, "\\\\\n");

	//add a X-line for each code of the family
	for (i=0; i<familySize; i++)
	{
		fprintf(f, "\\hline\n");
		fprintf(f, "$%s$ ",perf[i].tilingName);
		for (j=0; j<=maxWeight; j++)
			if ((!TrivialColumnZ[j]) || (!TrivialColumnX[j]))
			{
				if (i <= perf[j].maxWeight)
					fprintf(f, "& %d ", perf[i].XWeight[j]);
				else
					fprintf(f, "& 0 ");
			}
		fprintf(f, "\\\\\n");
	}
	fprintf(f, "\\hline\n");

	//finish the table
	fprintf(f, "\\end{tabular}\n\
\\caption{Measurements weight distribution}\n\
\\end{table}\n\
\\vspace{.3cm}\n\n\n");
}


//====================================
//			Remove and Compile
//====================================


//
void CleanLatexFile(char *fileName)
{
	char *command = malloc(100*sizeof(char));
	sprintf(command, "rm -f %s", fileName); //rm latex file
	system(command);

	free(command);
}


//compile latex minimal tikz figure
void GeneratePDF(char *jobName, char *fileName)
{
	char *command = malloc(200*sizeof(char));
	sprintf(command, "pdflatex --jobname=%s %s > texOutput.txt", jobName, fileName);
	system(command);

//	sprintf(command, "rm %s", fileName); //rm latex file
//	system(command);
	system("rm texOutput.txt");
	sprintf(command, "rm %s.log", jobName);	//rm latex log file
	system(command);
	
	free(command);
}


//compile latex file in folder report
void GenerateReportPDF(char *fileName)
{
	char *command = malloc(200*sizeof(char));
	sprintf(command, "cd reports && pdflatex %s > texOutput.txt", fileName);
	system(command);
	system("cd reports && rm texOutput.txt");
	sprintf(command, "cd reports && rm *.aux");
	system(command);
	sprintf(command, "cd reports && rm *.log");
	system(command);

	free(command);
}


//====================================
//			BASIC FIGURE
//====================================


//draw the graph structure without indices 
//with color c="r" for red or "b" for blue
void DrawMin(struct tiling G, char *pdfName, char c)
{
	char *fileName = malloc(100*sizeof(char));
	sprintf(fileName, "%s.tex", pdfName);
	
	char *jobName = malloc(100*sizeof(char));
	sprintf(jobName, "%s", pdfName);

	CleanLatexFile(fileName);
	
	FILE *f;
	f = fopen(fileName, "w+");

	//find scale and plot header
	WriteTikzMinHeader(f, jobName, FindScale(G));

	//define labels for vertices
	DefineLabels(G, f);

	//plot edges using vertices labels
	PlotEdges(G, f);

	//color face
	PlotFaces(G, f, c);

	//plot non-open vertices
	PlotNonOpenVertices(G, f);
	
	//plot open vertices
	PlotOpenVertices(G, f);

	//close file
	WriteTikzMinFooter(f);
	fclose(f);
	
	//generate pdf file and clean
	GeneratePDF(jobName, fileName);
}


//====================================
//			INDEXED FIGURE
//====================================


//draw the graph structure with indices 
//with color c="r" for red or "b" for blue
void Draw(struct tiling G, char *pdfName, char c)
{
	char *fileName = malloc(100*sizeof(char));
	sprintf(fileName, "%s.tex", pdfName);
	
	char *jobName = malloc(100*sizeof(char));
	sprintf(jobName, "%s", pdfName);

	CleanLatexFile(fileName);
	
	FILE *f;
	f = fopen(fileName, "w+");

	//find scale and plot header
	WriteTikzMinHeader(f, jobName, FindScale(G));

	//define labels
	DefineLabels(G, f);

	//plot indexed edges
	PlotIndexedEdges(G, f);

	//color face
	PlotFaces(G, f, c);
		
	//plot face indices	
	PlotFaceIndices(G, f);

	//plot vertex indices	
	PlotIndexedVertices(G, f);
	
	//close file
	WriteTikzMinFooter(f);
	fclose(f);
	
	//generate pdf file and clean
	GeneratePDF(jobName, fileName);
}



//====================================
//			HIGHLIGHT EDGES
//====================================


//draw the graph structure without indices 
//with color c="r" for red or "b" for blue
//draw in red a subset of edges
void HighlightEdges(struct tiling G, int *subsetE, char c)
{
	char *fileName = malloc(100*sizeof(char));
	sprintf(fileName, "figure.tex");

	char *jobName = malloc(100*sizeof(char));
	sprintf(jobName, "figure");

	CleanLatexFile(fileName);
	
	FILE *f;
	f = fopen(fileName, "w+");

	//find scale and plot header
	WriteTikzMinHeader(f, jobName, FindScale(G));

	//define labels for vertices
	DefineLabels(G, f);

	//plot edges using vertices labels
	PlotEdges(G, f);

	//color face
	PlotFaces(G, f, c);

	//plot non-open vertices
	PlotNonOpenVertices(G, f);
	
	//plot open vertices
	PlotOpenVertices(G, f);
	
	//plot in red the subset of edges
	PlotRedEdges(G, f, subsetE);
	
	//close file
	WriteTikzMinFooter(f);
	fclose(f);
	
	//generate pdf file and clean
	GeneratePDF(jobName, fileName);
}



//====================================
//		 Plot
//====================================

//Write the latex file and plot the performance curves
//the performance is tested in the range p in [xmin, xmax]
//the number of elements of the code family is mmax
//errorType = 1 for X or Z error, 2 for Z-error, 3 for X-error
void PlotPerformance(struct storePerformance *perf, int familySize, int errorType)
{
	int i;

	//check that there exist at least one non-empty file to plot
	int numberOfFiles = 0;
	for (i=0; i<familySize; i++)
		if (perf[i].dataFile != NULL)
			numberOfFiles++;

	if (numberOfFiles != 0)
	{
		char *fileName = malloc(50*sizeof(char));
		sprintf(fileName, "figure.tex");
		char *jobName = malloc(50*sizeof(char));
		sprintf(jobName, "figure");

		CleanLatexFile(fileName);
		FILE *f;
		f = fopen(fileName, "w+");

		//plot header
		float scale = 2;
		WriteTikzMinHeader(f, jobName, scale);

		//body
		char *dataFolder = "results/";
		PlotAxis(f, perf, familySize, dataFolder, errorType);

		//close file
		WriteTikzMinFooter(f);
		fclose(f);

		//generate pdf file and clean
		GeneratePDF(jobName, fileName);
	}
}


//====================================
//			Report
//====================================


//generate a pdf with
//(i) a figure of the tiling and its dual
//(ii) the code parameters and overhead
//(iii) the weights of the measurements
//(iv) the performance ploted for X and Z errors, X-errors and Z-errors
void Report(struct tiling G, struct tiling H, struct storePerformance *perf)
{
	char *reportFile = malloc(200*sizeof(char));
	sprintf(reportFile, "bench_%s_%ld.tex", perf->tilingName, time(0));
	char *addressReportFile = malloc(200*sizeof(char));
	sprintf(addressReportFile, "reports/%s", reportFile);
	CleanLatexFile(addressReportFile);
	
	FILE *f;
	f = fopen(addressReportFile, "w+");

	//write latex document header
	WriteLatexHeader(f);

	//draw the tiling and its dual when it contains at most 1000 edges
	if (G.e <= 1000)
	{
		//tiling DrawMin:
		BeginFigure(f);
		WriteTikzHeader(f, FindScale(G)/2.5, 0.25*FindScale(H)/2.5);
		DefineLabels(G, f);
		PlotEdges(G, f);
		PlotFaces(G, f, 'b');
		PlotNonOpenVertices(G, f);
		PlotOpenVertices(G, f);
		WriteTikzFooter(f);
		fprintf(f, "\\hspace{.2cm}\n");
	
		//dual tiling DrawMin:
		WriteTikzHeader(f, FindScale(H)/2.5, 0.25*FindScale(H)/2.5);
		DefineLabels(H, f);
		PlotEdges(H, f);
		PlotFaces(H, f, 'r');
		PlotNonOpenVertices(H, f);
		PlotOpenVertices(H, f);
		WriteTikzFooter(f);
		EndFigure(f, "The tiling and its dual tiling");	
	}
	
	//write code parameters
	fprintf(f, "\\newpage\n");
	WriteCodeParameters(f, perf->n, perf->k);
	
	//write measurement weights
	WriteWeights(f, perf->ZWeight, perf->XWeight, perf->maxWeight);

	//write performance plots if k>0
	if (perf->k != 0)
	{
		fprintf(f, "\\newpage\n");
		float scale = 1;
		BeginFigure(f);
		fprintf(f, "\\centering\n\n");
		WriteTikzHeader(f, scale, scale);
		int errorType = 1;
		PlotAxis(f, perf, 1, "../results/", errorType);
		WriteTikzFooter(f);
		EndFigure(f, "Logical error rate after decoding");
	
		BeginFigure(f);
		WriteTikzHeader(f, scale, scale);
		errorType = 2;
		PlotAxis(f, perf, 1, "../results/", errorType);
		WriteTikzFooter(f);
	
		WriteTikzHeader(f, scale, scale);
		errorType = 3;
		PlotAxis(f, perf, 1, "../results/", errorType);
		WriteTikzFooter(f);
		EndFigure(f, "Left: Z-logical error rate after decoding. Right: X-logical error rate after decoding");
	}
	
	//close file
	WriteLatexFooter(f);
	fclose(f);
	
	//generate pdf file and clean
	GenerateReportPDF(reportFile);
}

//====================================
//			Compare Codes
//====================================

void ComparisonReport(struct storePerformance *perf, int familySize)
{
	int i;

	char *comparisonFile = malloc(200*sizeof(char));
	sprintf(comparisonFile, "comparison_%ld.tex", time(0));
	char *addressComparisonFile = malloc(200*sizeof(char));
	sprintf(addressComparisonFile, "reports/%s", comparisonFile);
	CleanLatexFile(addressComparisonFile);

	FILE *f;
	f = fopen(addressComparisonFile, "w+");

	//write latex document header
	WriteLatexHeader(f);

	//write code parameters
	WriteComparisonCodeParameters(f, perf, familySize);

	//write measurement weights
	WriteComparisonWeights(f, perf, familySize);

	//plot performance if non-empty
	int numberOfFiles = 0;
	for (i=0; i<familySize; i++)
		if (perf[i].dataFile != NULL)
			numberOfFiles++;
	if (numberOfFiles != 0)
	{
		fprintf(f, "\\newpage\n");
		float scale = 1;
		int errorType = 1;
		BeginFigure(f);
		fprintf(f, "\\centering\n\n");
		WriteTikzHeader(f, scale, scale);
		PlotAxis(f, perf, familySize, "../results/", errorType);
		WriteTikzFooter(f);
		EndFigure(f, "Logical error rate after decoding");

		errorType = 2;
		BeginFigure(f);
		WriteTikzHeader(f, scale, scale);
		PlotAxis(f, perf, familySize, "../results/", errorType);
		WriteTikzFooter(f);

		errorType = 3;
		WriteTikzHeader(f, scale, scale);
		PlotAxis(f, perf, familySize, "../results/", errorType);
		WriteTikzFooter(f);
		EndFigure(f, "Left: Z-logical error rate after decoding. Right: X-logical error rate after decoding");
	}

	//close file
	WriteLatexFooter(f);
	fclose(f);

	//generate pdf file and clean
	GenerateReportPDF(comparisonFile);
}
