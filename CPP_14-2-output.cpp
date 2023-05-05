/************************************************************************************
Multi-Phase and Parallelized SPH program
--parallelized and optimized for SMP (Symmetric multiprocessing ) architecture.
--Copyright (c) 2018-2023, Weijie ZHANG, GeoHohai, Hohai University.
************************************************************************************/

#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include <cmath>
#include "string.h"
#include "Class_Functions.h"
#include "Header_Option.h"

/*Save the information to files*/
void clOutput_Fun::output_total(Particle *pVelPar, const Para_Pro &pPPro, char *argv)
{

	char file_name[100], temp[100];
	int i, temp1;
	double ux;

	/*File name*/
	FILE *fpt;

	memset(file_name, 0, 100);
	memset(temp, 0, 100);

	if (pPPro.l < 10)
	{
		sprintf(temp, "000000000%d", pPPro.l);
	}
	else if (pPPro.l < 100)
	{
		sprintf(temp, "00000000%d", pPPro.l);
	}
	else if (pPPro.l < 1000)
	{
		sprintf(temp, "0000000%d", pPPro.l);
	}
	else if (pPPro.l < 10000)
	{
		sprintf(temp, "000000%d", pPPro.l);
	}
	else if (pPPro.l < 100000)
	{
		sprintf(temp, "00000%d", pPPro.l);
	}
	else if (pPPro.l < 1000000)
	{
		sprintf(temp, "0000%d", pPPro.l);
	}
	else if (pPPro.l < 10000000)
	{
		sprintf(temp, "000%d", pPPro.l);
	}
	else if (pPPro.l < 100000000)
	{
		sprintf(temp, "00%d", pPPro.l);
	}

	/*Open files*/

	fpt = NULL;

	memset(file_name, 0, 100);
	strcpy(file_name, argv);
	if (win32)
		strcat(file_name, "\\data-");
	else
		strcat(file_name, "/data-");
	strcat(file_name, temp);
	strcat(file_name, ".odb");
	fpt = fopen(file_name, "w");
	fprintf(fpt, "*************************************************************************************************\n");
	fprintf(fpt, "%10d %14.6e\n", pPPro.l, pPPro.t);

	///< data//
	temp1 = 0;
	for (i = 0; i < pPPro.ntotal; i++)
	{
		if (pVelPar[i].type != 0)
		{
			ux = pVelPar[i].ux[0] * pVelPar[i].ux[0] + pVelPar[i].ux[1] * pVelPar[i].ux[1] + pVelPar[i].ux[2] * pVelPar[i].ux[2];
			ux = sqrt(ux);

			temp1 = temp1 + 1;
			fprintf(fpt, "%5d %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e \
%14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e \
%14.6e %14.6e %14.6e %14.6e %14.6e %14.6e \
%14.6e %14.6e %14.6e %14.6e %14.6e %14.6e \
%14.6e %14.6e %14.6e %14.6e %14.6e %14.6e \
%14.6e %14.6e %14.6e %14.6e %14.6e %14.6e \
%14.6e %14.6e %14.6e %14.6e  %14.6e %14.6e \
%14.6e %14.6e %14.6e \
%14.6e %4d %4d\n",
					temp1, pVelPar[i].xp[0], pVelPar[i].xp[1], pVelPar[i].xp[2], pVelPar[i].ux[0], pVelPar[i].ux[1], pVelPar[i].ux[2],
					ux, pVelPar[i].ax[0], pVelPar[i].ax[1], pVelPar[i].ax[2], pVelPar[i].vx[0], pVelPar[i].vx[1], pVelPar[i].vx[2],
					pVelPar[i].pre, pVelPar[i].prea, pVelPar[i].prew,
					pVelPar[i].sig[0], pVelPar[i].sig[1], pVelPar[i].sig[2], pVelPar[i].sig[3], pVelPar[i].sig[4], pVelPar[i].sig[5],
					pVelPar[i].strep[0], pVelPar[i].strep[1], pVelPar[i].strep[2],
					pVelPar[i].eps[0], pVelPar[i].eps[1], pVelPar[i].eps[2], pVelPar[i].eps[3], pVelPar[i].eps[4], pVelPar[i].eps[5],
					pVelPar[i].adps[0], pVelPar[i].adps[1], pVelPar[i].adps[2], pVelPar[i].adps[3], pVelPar[i].adps[4], pVelPar[i].adps[5],
					pVelPar[i].vsig, pVelPar[i].veps, pVelPar[i].meps, pVelPar[i].divq, pVelPar[i].divr, pVelPar[i].satu,
					pVelPar[i].suct, pVelPar[i].interfss[0], pVelPar[i].interfss[1], pVelPar[i].interfss[2], pVelPar[i].type, pVelPar[i].matype);
		}
		else
		{
			ux = pVelPar[i].ux[0] * pVelPar[i].ux[0] + pVelPar[i].ux[1] * pVelPar[i].ux[1] + pVelPar[i].ux[2] * pVelPar[i].ux[2];
			ux = sqrt(ux);

			temp1 = temp1 + 1;
			fprintf(fpt, "%5d %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e \
%14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e \
%14.6e %14.6e %14.6e %14.6e %14.6e %14.6e \
%14.6e %14.6e %14.6e %14.6e %14.6e %14.6e \
%14.6e %14.6e %14.6e %14.6e %14.6e %14.6e \
%14.6e %14.6e %14.6e %14.6e %14.6e %14.6e \
%14.6e %14.6e %14.6e %14.6e  %14.6e %14.6e \
%14.6e %14.6e %14.6e \
%14.6e %4d %4d\n",
					temp1, pVelPar[i].xp[0], pVelPar[i].xp[1], pVelPar[i].xp[2], pVelPar[i].ux[0], pVelPar[i].ux[1], pVelPar[i].ux[2],
					ux, pVelPar[i].interfss[0], pVelPar[i].interfss[1], pVelPar[i].interfss[2], pVelPar[i].interf[0], pVelPar[i].interf[1], pVelPar[i].interf[2],
					pVelPar[i].pre, pVelPar[i].prea, pVelPar[i].prew,
					pVelPar[i].sig[0], pVelPar[i].sig[1], pVelPar[i].sig[2], pVelPar[i].sig[3], pVelPar[i].sig[4], pVelPar[i].sig[5],
					pVelPar[i].strep[0], pVelPar[i].strep[1], pVelPar[i].strep[2],
					pVelPar[i].eps[0], pVelPar[i].eps[1], pVelPar[i].eps[2], pVelPar[i].eps[3], pVelPar[i].eps[4], pVelPar[i].eps[5],
					pVelPar[i].adps[0], pVelPar[i].adps[1], pVelPar[i].adps[2], pVelPar[i].adps[3], pVelPar[i].adps[4], pVelPar[i].adps[5],
					pVelPar[i].vsig, pVelPar[i].veps, pVelPar[i].meps, pVelPar[i].divq, pVelPar[i].divr, pVelPar[i].satu,
					pVelPar[i].suct, pVelPar[i].interfss[0], pVelPar[i].interfss[1], pVelPar[i].interfss[2], pVelPar[i].type, pVelPar[i].matype);
		}
	}
	fprintf(fpt, "*************************************************************************************************\n");
	/*Close files*/
	if (fpt != NULL)
		fclose(fpt);
}

/*V 2.9 Save the information of impact to files*/
void clOutput_Fun::output_impact(FILE *fpo, Particle *pVelPar, const Para_Pro &pPPro, char *argv)
{

	char file_name[100], temp[100];
	int i, temp1;
	int ndim = pPPro.ndim;
	double total_f;
	double dr = pPPro.dr;

	memset(file_name, 0, 100);
	memset(temp, 0, 100);

	/*File name*/
	FILE *fpi;

	if (pPPro.l < 10)
	{
		sprintf(temp, "000000000%d", pPPro.l);
	}
	else if (pPPro.l < 100)
	{
		sprintf(temp, "00000000%d", pPPro.l);
	}
	else if (pPPro.l < 1000)
	{
		sprintf(temp, "0000000%d", pPPro.l);
	}
	else if (pPPro.l < 10000)
	{
		sprintf(temp, "000000%d", pPPro.l);
	}
	else if (pPPro.l < 100000)
	{
		sprintf(temp, "00000%d", pPPro.l);
	}
	else if (pPPro.l < 1000000)
	{
		sprintf(temp, "0000%d", pPPro.l);
	}
	else if (pPPro.l < 10000000)
	{
		sprintf(temp, "000%d", pPPro.l);
	}
	else if (pPPro.l < 100000000)
	{
		sprintf(temp, "00%d", pPPro.l);
	}

	/*Open files*/
	fpi = NULL;

	strcpy(file_name, argv);
	if (win32)
		strcat(file_name, "\\Impact-force");
	else
		strcat(file_name, "/Impact-force");
	strcat(file_name, temp);
	strcat(file_name, ".inp");
	fpi = fopen(file_name, "w");

	temp1 = 1;
	total_f = 0.0;

	fprintf(fpi, "No.  X   Y  Z  Pre_Dens  Pre_water  Pre_Air  SIGXX  SIGYY  SIGZZ  SIGXY  SIGYZ  SIGZX\n");

	for (i = 0; i < pPPro.ntotal; i++)
	{
		if (pVelPar[i].type == 0 && pVelPar[i].matype >= 1)
		{
			fprintf(fpi, "%5d %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e \n", temp1, pVelPar[i].xp[0],
					pVelPar[i].xp[1], pVelPar[i].xp[2], pVelPar[i].pre, pVelPar[i].prew, pVelPar[i].prea,
					pVelPar[i].sig[0], pVelPar[i].sig[1], pVelPar[i].sig[2], pVelPar[i].sig[3],
					pVelPar[i].sig[4], pVelPar[i].sig[5]);
			temp1 += 1;
			total_f += pVelPar[i].prew * pow(dr, ndim - 1);
		}
	}
	fprintf(fpo, "Step%9d, Total force is %14.6e N\n", pPPro.l, total_f);
	/*Close files*/
	fclose(fpi);
}

/*V 2.9 Save the information of flowing distance of each phase to files*/
void clOutput_Fun::output_flowdist(FILE *fpi, Particle *pVelPar, const Para_Pro &pPPro)
{
	int i;
	double water_min[3], water_max[3];
	double soil_min[3], soil_max[3];
	double air_min[3], air_max[3];
	double str_min[3], str_max[3];

	//initializing
	for (i = 0; i < 3; i++)
	{
		water_min[i] = 10000.0;
		water_max[i] =-10000.0;

		soil_min[i] = 10000.0;
		soil_max[i] = -10000.0;

		air_min[i] = 10000.0;
		air_max[i] = -10000.0;

		str_min[i] = 10000.0;
		str_max[i] = -10000.0;
	}

	//calculate the maximum and minimum value of water, soil and air
	for (i = 0; i < pPPro.ntotal; i++)
	{
		if (pVelPar[i].type == 1)
		{
			water_max[0] = fmax(pVelPar[i].xp[0], water_max[0]);
			water_min[0] = fmin(pVelPar[i].xp[0], water_min[0]);

			water_max[1] = fmax(pVelPar[i].xp[1], water_max[1]);
			water_min[1] = fmin(pVelPar[i].xp[1], water_min[1]);

			water_max[2] = fmax(pVelPar[i].xp[2], water_max[2]);
			water_min[2] = fmin(pVelPar[i].xp[2], water_min[2]);
		}
		else if (pVelPar[i].type == 2)
		{
			soil_max[0] = fmax(pVelPar[i].xp[0], soil_max[0]);
			soil_min[0] = fmin(pVelPar[i].xp[0], soil_min[0]);

			soil_max[1] = fmax(pVelPar[i].xp[1], soil_max[1]);
			soil_min[1] = fmin(pVelPar[i].xp[1], soil_min[1]);

			soil_max[2] = fmax(pVelPar[i].xp[2], soil_max[2]);
			soil_min[2] = fmin(pVelPar[i].xp[2], soil_min[2]);
		}
		else if (pVelPar[i].type == 3)
		{
			air_max[0] = fmax(pVelPar[i].xp[0], air_max[0]);
			air_min[0] = fmin(pVelPar[i].xp[0], air_min[0]);

			air_max[1] = fmax(pVelPar[i].xp[1], air_max[1]);
			air_min[1] = fmin(pVelPar[i].xp[1], air_min[1]);

			air_max[2] = fmax(pVelPar[i].xp[2], air_max[2]);
			air_min[2] = fmin(pVelPar[i].xp[2], air_min[2]);
		}
		else if (pVelPar[i].type == 4)
		{
			str_max[0] = fmax(pVelPar[i].xp[0], str_max[0]);
			str_min[0] = fmin(pVelPar[i].xp[0], str_min[0]);

			str_max[1] = fmax(pVelPar[i].xp[1], str_max[1]);
			str_min[1] = fmin(pVelPar[i].xp[1], str_min[1]);

			str_max[2] = fmax(pVelPar[i].xp[2], str_max[2]);
			str_min[2] = fmin(pVelPar[i].xp[2], str_min[2]);
		}
	}

	for (i = 0; i < 3; i++)
	{
		if (fabs(water_min[i] - 10000.0) < 0.0001) water_min[i] = 0.0;
		if (fabs(water_max[i] + 10000.0) < 0.0001) water_max[i] = 0.0;
		if (fabs(soil_min[i] - 10000.0) < 0.0001) soil_min[i] = 0.0;
		if (fabs(soil_max[i] + 10000.0) < 0.0001) soil_max[i] = 0.0;
		if (fabs(air_min[i] - 10000.0) < 0.0001) air_min[i] = 0.0;
		if (fabs(air_max[i] + 10000.0) < 0.0001) air_max[i] = 0.0;
		if (fabs(str_min[i] - 10000.0) < 0.0001) str_min[i] = 0.0;
		if (fabs(str_max[i] + 10000.0) < 0.0001) str_max[i] = 0.0;
	}

	//output
	fprintf(fpi, "\nStep%9d            X                       Y                      Z\n", pPPro.l);
	fprintf(fpi, "-----------------------------------------------------------------------------------\n");
	fprintf(fpi, "                  max           min        max         min        max         min   \n");

	fprintf(fpi, "Water Phase:   %5.3e    %5.3e   %5.3e   %5.3e   %5.3e   %5.3e \n", water_max[0],
			water_min[0], water_max[1], water_min[1], water_max[2], water_min[2]);
	fprintf(fpi, "Soil Phase:    %5.3e    %5.3e   %5.3e   %5.3e   %5.3e   %5.3e \n", soil_max[0],
			soil_min[0], soil_max[1], soil_min[1], soil_max[2], soil_min[2]);
	fprintf(fpi, "Air Phase:     %5.3e    %5.3e   %5.3e   %5.3e   %5.3e   %5.3e \n", air_max[0],
			air_min[0], air_max[1], air_min[1], air_max[2], air_min[2]);
	fprintf(fpi, "Structure:      %5.3e    %5.3e   %5.3e   %5.3e   %5.3e   %5.3e \n", str_max[0],
			str_min[0], str_max[1], str_min[1], str_max[2], str_min[2]);
}

void clOutput_Fun::output_status(int step, int start, int end)
{
	//printf("It is %8d step of total %8d steps at present for output.\n", pPPro->l, pPPro->loop);//
	int i, loop;

	loop = (int)(21 * (step - start) / (end - start));

	printf("Step: %8d", step);
	printf("[%8d]", start);

	for (i = 0; i < 21; i++)
	{
		if (i <= loop)
			printf("=");
		else
			printf(" ");
	}
	printf("[%8d]\r", end);
}
