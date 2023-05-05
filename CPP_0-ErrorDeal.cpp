/************************************************************************************
Multi-Phase and Parallelized SPH program
--parallelized and optimized for SMP (Symmetric multiprocessing ) architecture.
--Copyright (c) 2018-2023, Weijie ZHANG, GeoHohai, Hohai University.
************************************************************************************/

#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include "Class_Functions.h"

void clErr_Fun::Err_Debug(int step, int pNum, FILE *flog)
{
	time_t rawtime;
	struct tm *timeinfo;

	time(&rawtime);
	timeinfo = localtime(&rawtime);

	printf("\n");

	printf("Current time: %s", asctime(timeinfo));
	printf("-------------------------------------------------------------------\n");
	printf("----||Error in Step NO. %9d and Particle No. %6d.\n", step, pNum);
	printf("-------------------------------------------------------------------\n");

	fprintf(flog, "Current time: %s", asctime(timeinfo));
	fprintf(flog, "-------------------------------------------------------------------\n");
	fprintf(flog, "----||Error in Step NO. %9d and Particle No. %6d.\n", step, pNum);
	fprintf(flog, "-------------------------------------------------------------------\n");
}

void clErr_Fun::Err_Deal(int err_type, FILE *flog)
{
	time_t rawtime;
	struct tm *timeinfo;

	time(&rawtime);
	timeinfo = localtime(&rawtime);

	printf("\n");

	if (err_type == 1)
	{
		printf("Current time: %s", asctime(timeinfo));
		printf("-------------------------------------------------------------------\n");
		printf("----||E-F01 There is no input.dat file or vibra.dat file.\n");
		printf("-------------------------------------------------------------------\n");

		fprintf(flog, "Current time: %s", asctime(timeinfo));
		fprintf(flog, "-------------------------------------------------------------------\n");
		fprintf(flog, "----||E-F01 There is no input.dat file or vibra.dat file.\n");
		fprintf(flog, "-------------------------------------------------------------------\n");

		fclose(flog);
	}
	else if (err_type == 2)
	{
		printf("Current time: %s", asctime(timeinfo));
		printf("-------------------------------------------------------------------\n");
		printf("----||E-M02 The number of particles is less than one.\n");
		printf("-------------------------------------------------------------------\n");

		fprintf(flog, "Current time: %s", asctime(timeinfo));
		fprintf(flog, "-------------------------------------------------------------------\n");
		fprintf(flog, "----||E-M02 The number of particles is less than one.\n");
		fprintf(flog, "-------------------------------------------------------------------\n");

		fclose(flog);
	}
	else if (err_type == 3)
	{
		printf("Current time: %s", asctime(timeinfo));
		printf("-------------------------------------------------------------------\n");
		printf("----||E-C03 Cell ID of particle is wrong.\n");
		printf("-------------------------------------------------------------------\n");

		fprintf(flog, "Current time: %s", asctime(timeinfo));
		fprintf(flog, "-------------------------------------------------------------------\n");
		fprintf(flog, "----||E-C03 Cell ID of particle is wrong.\n");
		fprintf(flog, "-------------------------------------------------------------------\n");

		fclose(flog);
	}
	else if (err_type == 4)
	{
		printf("Current time: %s", asctime(timeinfo));
		printf("-------------------------------------------------------------------\n");
		printf("----||E-C04 Number of particles in cell is too large.\n");
		printf("-------------------------------------------------------------------\n");

		fprintf(flog, "Current time: %s", asctime(timeinfo));
		fprintf(flog, "-------------------------------------------------------------------\n");
		fprintf(flog, "----||E-C04 Number of particles in cell is too large.\n");
		fprintf(flog, "-------------------------------------------------------------------\n");

		fclose(flog);
	}
	else if (err_type == 5)
	{
		printf("Current time: %s", asctime(timeinfo));
		printf("-------------------------------------------------------------------\n");
		printf("----||E-S05 Number of supporting particles is larger than 99.\n");
		printf("-------------------------------------------------------------------\n");

		fprintf(flog, "Current time: %s", asctime(timeinfo));
		fprintf(flog, "-------------------------------------------------------------------\n");
		fprintf(flog, "----||E-S05 Number of supporting particles is larger than 99.\n");
		fprintf(flog, "-------------------------------------------------------------------\n");

		fclose(flog);
	}
	else if (err_type == 6)
	{
		printf("Current time: %s", asctime(timeinfo));
		printf("-------------------------------------------------------------------\n");
		printf("----||E-R06 The type of rainfall is larger than 7 or less than 1.\n");
		printf("-------------------------------------------------------------------\n");

		fprintf(flog, "Current time: %s", asctime(timeinfo));
		fprintf(flog, "-------------------------------------------------------------------\n");
		fprintf(flog, "----||E-R06 The type of rainfall is larger than 7 or less than 1.\n");
		fprintf(flog, "-------------------------------------------------------------------\n");

		fclose(flog);
	}
	else if (err_type == 7)
	{
		printf("Current time: %s", asctime(timeinfo));
		printf("-------------------------------------------------------------------\n");
		printf("----||E-V07 The total step of vibration is less than 1.\n");
		printf("-------------------------------------------------------------------\n");

		fprintf(flog, "Current time: %s", asctime(timeinfo));
		fprintf(flog, "-------------------------------------------------------------------\n");
		fprintf(flog, "----||E-V07 The total step of vibration is less than 1.\n");
		fprintf(flog, "-------------------------------------------------------------------\n");

		fclose(flog);
	}
	else if (err_type == 8)
	{
		printf("Current time: %s", asctime(timeinfo));
		printf("-------------------------------------------------------------------\n");
		printf("----||E-EB08 The number of excavation or backfill stage is less than 1.\n");
		printf("-------------------------------------------------------------------\n");

		fprintf(flog, "Current time: %s", asctime(timeinfo));
		fprintf(flog, "-------------------------------------------------------------------\n");
		fprintf(flog, "----||E-EB08 The number of excavation or backfill stage is less than 1.\n");
		fprintf(flog, "-------------------------------------------------------------------\n");

		fclose(flog);
	}
	else if (err_type == 9)
	{
		printf("Current time: %s", asctime(timeinfo));
		printf("-------------------------------------------------------------------\n");
		printf("----||E-WW09: The type of wind or flowing water is not correct.\n");
		printf("-------------------------------------------------------------------\n");

		fprintf(flog, "Current time: %s", asctime(timeinfo));
		fprintf(flog, "-------------------------------------------------------------------\n");
		fprintf(flog, "----||E-WW09: The type of gas flow or liquid flow is not correct.\n");
		fprintf(flog, "-------------------------------------------------------------------\n");

		fclose(flog);
	}
	else if (err_type == 10)
	{
		printf("Current time: %s", asctime(timeinfo));
		printf("-------------------------------------------------------------------\n");
		printf("----||E-WD10: The working directory is not correct.\n");
		printf("-------------------------------------------------------------------\n");
	}
	else if (err_type == 11)
	{
		printf("Current time: %s", asctime(timeinfo));
		printf("-------------------------------------------------------------------\n");
		printf("----||E-WD11: The working directory is not correct.\n");
		printf("-------------------------------------------------------------------\n");

		fprintf(flog, "Current time: %s", asctime(timeinfo));
		fprintf(flog, "-------------------------------------------------------------------\n");
		fprintf(flog, "----||E-WD11: The working directory is not correct.\n");
		fprintf(flog, "-------------------------------------------------------------------\n");

		fclose(flog);
	}
	else if (err_type == 12)
	{
		printf("Current time: %s", asctime(timeinfo));
		printf("-------------------------------------------------------------------\n");
		printf("----||E-time12: The loops for designated time is too large.\n");
		printf("-------------------------------------------------------------------\n");

		fprintf(flog, "Current time: %s", asctime(timeinfo));
		fprintf(flog, "-------------------------------------------------------------------\n");
		fprintf(flog, "----||E-time12: The loops for designated time is too large.\n");
		fprintf(flog, "-------------------------------------------------------------------\n");

		fclose(flog);
	}
	else if (err_type == 13)
	{
		printf("Current time: %s", asctime(timeinfo));
		printf("-------------------------------------------------------------------\n");
		printf("----||E-SpatialV13: The direction of distribution is not correct.\n");
		printf("-------------------------------------------------------------------\n");

		fprintf(flog, "Current time: %s", asctime(timeinfo));
		fprintf(flog, "-------------------------------------------------------------------\n");
		fprintf(flog, "----||E-SpatialV13: The direction of distribution is not correct.\n");
		fprintf(flog, "-------------------------------------------------------------------\n");

		fclose(flog);
	}
	else if (err_type == 14)
	{
		printf("Current time: %s", asctime(timeinfo));
		printf("-------------------------------------------------------------------\n");
		printf("----||E-SpatialV14: The distribution (SV->type) must be either Standard\
Normal (1) or Log normal (2).\n");
		printf("-------------------------------------------------------------------\n");

		fprintf(flog, "Current time: %s", asctime(timeinfo));
		fprintf(flog, "-------------------------------------------------------------------\n");
		fprintf(flog, "----||E-SpatialV14: The distribution (SV->type) must be either Standard\
Normal (1) or Log normal (2).\n");
		fprintf(flog, "-------------------------------------------------------------------\n");

		fclose(flog);
	}
	else if (err_type == 15)
	{
		printf("Current time: %s", asctime(timeinfo));
		printf("-------------------------------------------------------------------\n");
		printf("----||E-SpatialV15: The distribution (SV->ddtype) is not correct.\n");
		printf("-------------------------------------------------------------------\n");

		fprintf(flog, "Current time: %s", asctime(timeinfo));
		fprintf(flog, "-------------------------------------------------------------------\n");
		fprintf(flog, "----||E-SpatialV15: The distribution (SV->ddtype) is not correct.\n");
		fprintf(flog, "-------------------------------------------------------------------\n");

		fclose(flog);
	}
	else if (err_type == 16)
	{
		printf("Current time: %s", asctime(timeinfo));
		printf("-------------------------------------------------------------------\n");
		printf("----||E-Kernel16: The flknl is not correct (must be 1 or 2).\n");
		printf("-------------------------------------------------------------------\n");

		fprintf(flog, "Current time: %s", asctime(timeinfo));
		fprintf(flog, "-------------------------------------------------------------------\n");
		fprintf(flog, "----||E-Kernel16: The flknl is not correct (must be 1 or 2)..\n");
		fprintf(flog, "-------------------------------------------------------------------\n");

		fclose(flog);
	}
	else if (err_type == 42)
	{
	printf("Current time: %s", asctime(timeinfo));
	printf("-------------------------------------------------------------------\n");
	printf("----||E-SPVRF42: Too many particles in the cells for the random field generation.\n");
	printf("-------------------------------------------------------------------\n");

	fprintf(flog, "Current time: %s", asctime(timeinfo));
	fprintf(flog, "-------------------------------------------------------------------\n");
	fprintf(flog, "----||E-SPVRF42: Too many particles in the cells for the random field generation.\n");
	fprintf(flog, "-------------------------------------------------------------------\n");

	fclose(flog);
	}
}
