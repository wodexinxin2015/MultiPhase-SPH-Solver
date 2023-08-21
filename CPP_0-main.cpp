/************************************************************************************
Multi-Phase and Parallelized SPH program
--parallelized and optimized for SMP (Symmetric multiprocessing) architecture.
--optimized for AVX1 or AVX2 instructions.
--Copyright (c) 2018-2023, Weijie ZHANG, GeoHohai, Hohai University.
************************************************************************************/
/************************************************************************************
revision:
2023-03-20: V4.13.1.0
	reduce the memory consumption of NNPS variables.
2023-03-24: V4.13.2.0
	add the correction algorithm of CSPM into the density, strain and XSPH calculation.
2023-03-28: V4.13.3.0
	revise the constitutive model of DP from DBLeaves, DP of Bui and Elastic-Plastic of MC;
	revise the determiniation of time step increment.
2023-04-02: V4.13.4.0
	change the module of particle_into_cell to the Binary sorting method in NNPS.
2023-06-19: V4.13.5.0
	Revise the DP model and Bui-DP model.
	Revise the segmentation fault bug in the CPP_3-NNPS.cpp.
2023-06-26: V4.14.1.0b
	in development
2023-08-18£ºV5.1.1
    start a new version number.
2023-08-20£ºV5.1.2b
	revise a bug in the CPP_3-NNPS.cpp: cellid for the rainfall particle is not calculated,
	but checked in the following code, thus it will throw an exception.
2023-08-21£ºV5.1.2c
	revise a bug in the CPP_1-SPHLoops.cpp: mod(pPPro->l - pPPro->inip, lp) == 0
	------->(pPPro->l - pPPro->inip) % lp == 0
************************************************************************************/

#define _CRT_SECURE_NO_DEPRECATE
#include <stdlib.h>
#include <stdio.h>
#include "string.h"
#include <time.h>
#include "Header_Option.h"

int SPH_Calculation(int argc, char *argv);

int main(int argc, char *argv[])
{

	char *direct_name, temp[100];
	int i;
	time_t start, end;
	int result, run_result;
	FILE *fpch;

	run_result = 0;

	/*SPH simulation process with arguments: argc, argv*/
	if (argc > 2)
	{

		printf("-------------------------------------------------------------------\n");
		printf("--||The simulation has    %3d working project.||----------\n", (argc - 2) < 1 ? 1 : (argc - 2));
		printf("-------------------------------------------------------------------\n");

		if (argv[1][0] == '-' && argv[1][1] == 'p')
		{
			for (i = 2; i < argc; i++)
			{
				printf("-------------------------------------------------------------------\n");
				printf("--||The Working Project No.%3d is now running.||----------\n", i - 1);
				printf("--<%53s>--\n", argv[i]);
				printf("-------------------------------------------------------------------\n");

				memset(temp, 0, 100);
				direct_name = argv[i];

				/*start calculation*/
				start = time(NULL);

				run_result = SPH_Calculation(argc, argv[i]);

				/*consumed time calculation*/
				end = time(NULL);
				result = (int)difftime(end, start);

				//path for windows
				if (win32)
				{
					//file for log output
					strcpy(temp, direct_name);
					strcat(temp, "\\Consumed time.txt");
				}
				//path for Linux
				else
				{
					//file for log output
					strcpy(temp, direct_name);
					strcat(temp, "/Consumed time.txt");
				}

				if ((fpch = fopen(temp, "w")) != NULL)
				{
					fprintf(fpch, "Time used for simulation:%3d : %2d : %2d\n", (int)(result / 3600),
							(int)((result % 3600) / 60), (int)((result % 3600) % 60));
					fclose(fpch);
				}

				if (run_result == 0)
				{
					printf("-------------------------------------------------------------------\n");
					printf("--||The Working Project No.%3d is finished.||-------------\n", i - 1);
					printf("-------------------------------------------------------------------\n");
				}
				else
				{
					printf("-------------------------------------------------------------------\n");
					printf("--||Simulation of Project No.%3d raised an error.|| ------\n", i - 1);
					printf("-------------------------------------------------------------------\n");
				}
			}
		}
		else
		{
			printf("Error 200: option is not -p.\n");
			printf("-------------------------------------------------------------------\n");
			printf("Multi-Phase and Parallelized SPH program\n");
			printf("--Copyright (c) 2018-2023, Weijie ZHANG, GeoHohai, Hohai University.\n");
			printf("-------------------------------------------------------------------\n");
			printf("     Please using no arguments or following arguments.\n");
			printf("-------------------------------------------------------------------\n");
			printf("     Argument list: -p directory1 directory2 ...\n");
			printf("               -p         ---specifying working directory;\n");
			printf("               directory1 ---directory path 1;\n");
			printf("               directory2 ---directory path 2;\n");
			printf("               ...        ---more directory path.\n");
			printf("-------------------------------------------------------------------\n");

			return 200;
		}
	}
	else if (argc == 2)
	{
		printf("Error 100: Argument is not correct, Exit.\n");
		printf("-------------------------------------------------------------------\n");
		printf("Multi-Phase and Parallelized SPH program\n");
		printf("--Copyright (c) 2018-2023, Weijie ZHANG, GeoHohai, Hohai University.\n");
		printf("-------------------------------------------------------------------\n");
		printf("     Please using no arguments or following arguments.\n");
		printf("-------------------------------------------------------------------\n");
		printf("     Argument list: -p directory1 directory2 ...\n");
		printf("               -p         ---specifying working directory;\n");
		printf("               directory1 ---directory path 1;\n");
		printf("               directory2 ---directory path 2;\n");
		printf("               ...        ---more directory path.\n");
		printf("-------------------------------------------------------------------\n");

		return 100;
	}
	else if (argc == 1)
	{

		printf("-------------------------------------------------------------------\n");
		printf("---||No specifying working directory.||--------------------\n");
		printf("---||using the directory of executable file.||-------------\n");
		printf("-------------------------------------------------------------------\n");

		printf("-------------------------------------------------------------------\n");
		printf("---||The Working Project is running.||---------------------\n");
		printf("-------------------------------------------------------------------\n");

		memset(temp, 0, 100);

		/*start calculation*/
		start = time(NULL);

		run_result = SPH_Calculation(argc, (char *)("."));

		/*consumed time calculation*/
		end = time(NULL);
		result = (int)difftime(end, start);

		//path for windows
		if (win32)
		{
			//file for log output
			strcpy(temp, ".");
			strcat(temp, "\\Consumed time.txt");
		}
		//path for Linux
		else
		{
			//file for log output
			strcpy(temp, ".");
			strcat(temp, "/Consumed time.txt");
		}

		if ((fpch = fopen(temp, "w")) != NULL)
		{
			fprintf(fpch, "Time used for simulation:%3d : %2d : %2d\n", (int)(result / 3600),
					(int)((result % 3600) / 60), (int)((result % 3600) % 60));
			fclose(fpch);
		}

		if (run_result == 0)
		{
			printf("-------------------------------------------------------------------\n");
			printf("--||The Working Project No.   1 is finished.||-------------\n");
			printf("-------------------------------------------------------------------\n");
		}
		else
		{
			printf("-------------------------------------------------------------------\n");
			printf("--||Simulation of Project No.  1 raised an error.|| ------\n");
			printf("-------------------------------------------------------------------\n");
		}
	}

	return 0;
}
