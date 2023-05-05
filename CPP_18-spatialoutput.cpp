/************************************************************************************
Multi-Phase and Parallelized SPH program
--parallelized and optimized for SMP (Symmetric multiprocessing ) architecture.
--Copyright (c) 2018-2023, Weijie ZHANG, GeoHohai, Hohai University.
************************************************************************************/
#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <time.h>
#include <string.h>
#include "Header_Parameters.h"
#include "Class_particle.h"
#include "Class_ParticleProperty.h"
#include "Class_Problem.h"
#include "Class_spatial_variability.h"
#include "Header_Option.h"

void clFParaOutput_Fun::paramters_output(Particle *pPar, Para_Fluid *pVFluid, Para_Soil *pParti_ConsPara,
										 const Para_Pro &pPPro, char *argv)
{

	char file_name[100];
	int i, err;
	FILE *fpt;

	//open file
	memset(file_name, 0, 100);
	strcpy(file_name, argv);
	if (win32)
		strcat(file_name, "\\Parameters");
	else
		strcat(file_name, "/Parameters");
	strcat(file_name, ".txt");
	fpt = fopen(file_name, "w");

	//write file
	err = fprintf(fpt, "No.---- X----Y----Z----fai----c----cop----ds----type----matype----\n");
	for (i = 0; i < pPPro.ntotal; i++)
	{
		if (pPar[i].type == 2 || pPar[i].type == 4)
		{
			err = fprintf(fpt, "   %6d   %10.6e   %10.6e   %10.6e   %10.6e   %10.6e   %10.6e   %10.6e   %2d   %2d\n",
				i + 1, pPar[i].xp[0], pPar[i].xp[1], pPar[i].xp[2], pParti_ConsPara[i].fai,
				pParti_ConsPara[i].c, pParti_ConsPara[i].cop, pParti_ConsPara[i].ds, pPar[i].type, pPar[i].matype);
		}
		else if (pPar[i].type == 1 || pPar[i].type == 3)
		{
			err = fprintf(fpt, "   %6d   %10.6e   %10.6e   %10.6e   %10.6e   %10.6e   %10.6e   %10.6e   %2d   %2d\n",
				i + 1, pPar[i].xp[0], pPar[i].xp[1], pPar[i].xp[2], pVFluid[i].fai,
				pVFluid[i].c, pParti_ConsPara[i].cop, pParti_ConsPara[i].ds, pPar[i].type, pPar[i].matype);
		}
	}
	err = fprintf(fpt, "No.---- X----Y----Z----fai----c----cop----ds----type----matype----\n");

	/*Close files*/
	if (fpt != NULL)
		fclose(fpt);
}
