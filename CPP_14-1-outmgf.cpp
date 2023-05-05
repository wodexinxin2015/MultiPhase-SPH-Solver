/************************************************************************************
Multi-Phase and Parallelized SPH program
--parallelized and optimized for SMP (Symmetric multiprocessing ) architecture.
--Copyright (c) 2018-2023, Weijie ZHANG, GeoHohai, Hohai University.
************************************************************************************/

#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include "Class_Functions.h"
#include "Header_Parameters.h"

void clOutput_Fun::mgf(Particle *pVelPar, const Para_Pro &pPPro, const Cell_Con &pCellc, FILE *fp)
{
	int i;
	int ndim = pPPro.ndim;
	double bmax[3], bmin[3];
	int ntotal = pPPro.ntotal;
	double dr = 1.1 * pPPro.dr;

	bmax[0] = pCellc.xmax;
	bmax[1] = pCellc.ymax;
	bmax[2] = pCellc.zmax;
	bmin[0] = pCellc.xmin;
	bmin[1] = pCellc.ymin;
	bmin[2] = pCellc.zmin;

	//write step information
	fprintf(fp, "step%d\n", (int)(pPPro.l / pPPro.fop));

	//write boundary
	if (ndim == 2)
	{
		fprintf(fp, " disjoint line\n boundary\n color\n 8\n");

		fprintf(fp, "%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n", bmin[0], bmin[1], 0.0, 0.5, 0.5, 0.5);
		fprintf(fp, "%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n", bmax[0], bmin[1], 0.0, 0.5, 0.5, 0.5);

		fprintf(fp, "%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n", bmax[0], bmin[1], 0.0, 0.5, 0.5, 0.5);
		fprintf(fp, "%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n", bmax[0], bmax[1], 0.0, 0.5, 0.5, 0.5);

		fprintf(fp, "%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n", bmax[0], bmax[1], 0.0, 0.5, 0.5, 0.5);
		fprintf(fp, "%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n", bmin[0], bmax[1], 0.0, 0.5, 0.5, 0.5);

		fprintf(fp, "%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n", bmin[0], bmax[1], 0.0, 0.5, 0.5, 0.5);
		fprintf(fp, "%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n", bmin[0], bmin[1], 0.0, 0.5, 0.5, 0.5);
	}
	else if (ndim == 3)
	{
		fprintf(fp, " disjoint line\n boundary\n color\n 24\n");
		//bottom
		fprintf(fp, "%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n", bmin[0], bmin[1], bmin[2], 0.5, 0.5, 0.5);
		fprintf(fp, "%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n", bmax[0], bmin[1], bmin[2], 0.5, 0.5, 0.5);

		fprintf(fp, "%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n", bmax[0], bmin[1], bmin[2], 0.5, 0.5, 0.5);
		fprintf(fp, "%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n", bmax[0], bmax[1], bmin[2], 0.5, 0.5, 0.5);

		fprintf(fp, "%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n", bmax[0], bmax[1], bmin[2], 0.5, 0.5, 0.5);
		fprintf(fp, "%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n", bmin[0], bmax[1], bmin[2], 0.5, 0.5, 0.5);

		fprintf(fp, "%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n", bmin[0], bmax[1], bmin[2], 0.5, 0.5, 0.5);
		fprintf(fp, "%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n", bmin[0], bmin[1], bmin[2], 0.5, 0.5, 0.5);
		//top
		fprintf(fp, "%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n", bmin[0], bmin[1], bmax[2], 0.5, 0.5, 0.5);
		fprintf(fp, "%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n", bmax[0], bmin[1], bmax[2], 0.5, 0.5, 0.5);

		fprintf(fp, "%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n", bmax[0], bmin[1], bmax[2], 0.5, 0.5, 0.5);
		fprintf(fp, "%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n", bmax[0], bmax[1], bmax[2], 0.5, 0.5, 0.5);

		fprintf(fp, "%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n", bmax[0], bmax[1], bmax[2], 0.5, 0.5, 0.5);
		fprintf(fp, "%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n", bmin[0], bmax[1], bmax[2], 0.5, 0.5, 0.5);

		fprintf(fp, "%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n", bmin[0], bmax[1], bmax[2], 0.5, 0.5, 0.5);
		fprintf(fp, "%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n", bmin[0], bmin[1], bmax[2], 0.5, 0.5, 0.5);
		//vertical
		fprintf(fp, "%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n", bmin[0], bmin[1], bmin[2], 0.5, 0.5, 0.5);
		fprintf(fp, "%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n", bmin[0], bmin[1], bmax[2], 0.5, 0.5, 0.5);

		fprintf(fp, "%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n", bmax[0], bmin[1], bmin[2], 0.5, 0.5, 0.5);
		fprintf(fp, "%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n", bmax[0], bmin[1], bmax[2], 0.5, 0.5, 0.5);

		fprintf(fp, "%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n", bmax[0], bmax[1], bmin[2], 0.5, 0.5, 0.5);
		fprintf(fp, "%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n", bmax[0], bmax[1], bmax[2], 0.5, 0.5, 0.5);

		fprintf(fp, "%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n", bmin[0], bmax[1], bmin[2], 0.5, 0.5, 0.5);
		fprintf(fp, "%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n", bmin[0], bmax[1], bmax[2], 0.5, 0.5, 0.5);
	}

	fprintf(fp, " sphere\n particles\n color\n %d\n", ntotal);

	for (i = 0; i < ntotal; i++)
	{
		if (pVelPar[i].type == 0) //water color: blue 0.0, 0.0, 1.0
			fprintf(fp, "%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n", pVelPar[i].xp[0], pVelPar[i].xp[1], pVelPar[i].xp[2], dr, 0.80, 0.80, 0.80);
		if (pVelPar[i].type == 1) //water color: blue 0.0, 0.0, 1.0
			fprintf(fp, "%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n", pVelPar[i].xp[0], pVelPar[i].xp[1], pVelPar[i].xp[2], dr, 0.2, 0.47, 0.87);
		else if (pVelPar[i].type == 2) //soil color: yellow 1.0, 1.0, 0.0
			fprintf(fp, "%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n", pVelPar[i].xp[0], pVelPar[i].xp[1], pVelPar[i].xp[2], dr, 0.87, 0.87, 0.2);
		else if (pVelPar[i].type == 3) //air color: cyan 0.0, 1.0, 1.0
			fprintf(fp, "%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n", pVelPar[i].xp[0], pVelPar[i].xp[1], pVelPar[i].xp[2], dr, 0.0, 1.0, 0.67);
		else if (pVelPar[i].type == 4) //air color: cyan 0.0, 1.0, 1.0
			fprintf(fp, "%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n", pVelPar[i].xp[0], pVelPar[i].xp[1], pVelPar[i].xp[2], dr, 0.73, 0.47, 0.07);
		else if (pVelPar[i].type == 7) //rain color: blue 0.0, 0.0, 1.0
			fprintf(fp, "%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n", pVelPar[i].xp[0], pVelPar[i].xp[1], pVelPar[i].xp[2], dr, 0.5, 0.5, 0.5);
	}
}