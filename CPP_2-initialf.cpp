/************************************************************************************
Multi-Phase and Parallelized SPH program
--parallelized and optimized for SMP (Symmetric multiprocessing ) architecture.
--Copyright (c) 2018-2023, Weijie ZHANG, GeoHohai, Hohai University.
************************************************************************************/

#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <string.h>
#include <time.h>
#include "Class_Problem.h"
#include "Class_Functions.h"
#include "Class_spatial_variability.h"
#include "Header_Option.h"

extern int err_scanf;
double moisturesuct(double satu, double dsatu, const Para_Soil &Soil);

//initializing from input file for problems property
int clFIni_Fun::initial_pro(Para_FL *pFlag, Para_Pro *pPPro, Para_Rain *pPRain, Para_Inter *pPInter, Para_SSInter *pSSInter,
							Para_Boundary *pPBou, Para_Vibra *pPVib, Para_GF *pPGf, Para_EC *pPEc, double *threshold_min, double* threshold_max,
							FILE *flog, FILE *finp)
{
	FILE *pFile;
	char c[2];
	int i;
	time_t rawtime;
	struct tm *timeinfo;

	time(&rawtime);
	timeinfo = localtime(&rawtime);
	printf("Current time: %s", asctime(timeinfo));
	printf("The problem parameters are initializing\n.......\n");
	fprintf(flog, "Current time: %s", asctime(timeinfo));
	fprintf(flog, "The problem parameters are initializing\n.......\n");

	pFile = finp;
	rewind(pFile);

	while (!feof(pFile))
	{
		err_scanf = fscanf(pFile, "%c%c%*[^\n]%*c", &c[0], &c[1]);
		if (c[0] == 'F' && c[1] == 'L')
		{ //flags
			do
			{
				err_scanf = fscanf(pFile, "%c%*[^\n]%*c", &c[0]);
				if (c[0] == 'A')
				{

					err_scanf = fscanf(pFile, "%d%d%d%d%d%d%d%d%d", &pFlag->fla, &pFlag->fls, &pFlag->flr,
									   &pFlag->fli, &pFlag->flv, &pFlag->fleob, &pFlag->flts,
									   &pFlag->flspv, &pFlag->flknl);
					err_scanf = fscanf(pFile, "%d%d", &pFlag->ntd, &pFlag->flbndy);
					break;
				}
			} while (c[0] != 'E');
			break;
		}
	}

	rewind(pFile);
	while (!feof(pFile))
	{
		err_scanf = fscanf(pFile, "%c%c%*[^\n]%*c", &c[0], &c[1]);
		if (c[0] == 'P' && c[1] == 'R')
		{ //problem parameters
			do
			{
				err_scanf = fscanf(pFile, "%c%*[^\n]%*c", &c[0]);
				if (c[0] == 'A')
				{
					err_scanf = fscanf(pFile, "%d%lf%lf", &pPPro->ndim, threshold_min, threshold_max);
					err_scanf = fscanf(pFile, "%d%d%d%d%d", &pPPro->fop, &pPPro->inip, &pPPro->loop, &pPPro->ntotal, &pPPro->step_regu);
					err_scanf = fscanf(pFile, "%d%lf%lf", &pPPro->dttype, &pPPro->dr, &pPPro->dt);
					break;
				}
			} while (c[0] != 'E');
			break;
		}
	}

	rewind(pFile);
	while (!feof(pFile))
	{
		err_scanf = fscanf(pFile, "%c%c%*[^\n]%*c", &c[0], &c[1]);
		if (c[0] == 'X' && c[1] == 'P')
		{ //xsph parameters
			do
			{
				err_scanf = fscanf(pFile, "%c%*[^\n]%*c", &c[0]);
				if (c[0] == 'A')
				{
					err_scanf = fscanf(pFile, "%lf%lf%lf", &pPPro->xsph[0], &pPPro->xsph[1], &pPPro->xsph[2]);
					break;
				}
			} while (c[0] != 'E');
			break;
		}
	}

	rewind(pFile);
	while (!feof(pFile))
	{
		err_scanf = fscanf(pFile, "%c%c%*[^\n]%*c", &c[0], &c[1]);
		if (c[0] == 'A' && c[1] == 'T')
		{ //artificial viscosity parameters
			do
			{
				err_scanf = fscanf(pFile, "%c%*[^\n]%*c", &c[0]);
				if (c[0] == 'A')
				{
					err_scanf = fscanf(pFile, "%lf%lf", &pPPro->art_vis[0][0], &pPPro->art_vis[0][1]);							//water
					err_scanf = fscanf(pFile, "%lf%lf%lf", &pPPro->art_vis[1][0], &pPPro->art_vis[1][1], &pPPro->art_stre_coe); //soil
					//artificial viscosity and artificial stress
					err_scanf = fscanf(pFile, "%lf%lf", &pPPro->art_vis[2][0], &pPPro->art_vis[2][1]); //air
					break;
				}
			} while (c[0] != 'E');
			break;
		}
	}

	rewind(pFile);
	while (!feof(pFile))
	{
		err_scanf = fscanf(pFile, "%c%c%*[^\n]%*c", &c[0], &c[1]);
		if (c[0] == 'E' && c[1] == 'C' && pFlag->fleob == 1)
		{ //excavation or backfill parameters
			do
			{
				err_scanf = fscanf(pFile, "%c%*[^\n]%*c", &c[0]);
				if (c[0] == 'A')
				{
					err_scanf = fscanf(pFile, "%d", &pPEc->ecount);

					//initializing of excavation or backfill parameters
					if (pPEc->ecount >= 1)
					{
						pPEc->memloads();
						for (i = 0; i < pPEc->ecount; i++)
						{

							err_scanf = fscanf(pFile, "%d%d", &pPEc->ecp[i][0], &pPEc->ecp[i][1]);
						}
					}
					else
					{
						return 8; //excavation or backfill steps error
					}

					break;
				}
			} while (c[0] != 'E');
			break;
		}
	}

	rewind(pFile);
	while (!feof(pFile))
	{
		err_scanf = fscanf(pFile, "%c%c%*[^\n]%*c", &c[0], &c[1]);
		if (c[0] == 'G' && c[1] == 'F')
		{ //gravity field
			do
			{
				err_scanf = fscanf(pFile, "%c%*[^\n]%*c", &c[0]);
				if (c[0] == 'A')
				{
					err_scanf = fscanf(pFile, "%lf%lf%lf", &pPGf->gx, &pPGf->gy, &pPGf->gz);
					break;
				}
			} while (c[0] != 'E');
			break;
		}
	}

	rewind(pFile);
	while (!feof(pFile))
	{
		err_scanf = fscanf(pFile, "%c%c%*[^\n]%*c", &c[0], &c[1]);
		if (c[0] == 'R' && c[1] == 'P' && pFlag->flr == 1)
		{ //rainfall parameters
			do
			{
				err_scanf = fscanf(pFile, "%c%*[^\n]%*c", &c[0]);
				if (c[0] == 'A')
				{
					err_scanf = fscanf(pFile, "%d%d%d%lf%lf",
									   &pPRain->rtype, &pPRain->nrain, &pPRain->nhor, &pPRain->press_factor, &pPRain->drop_velocity);
					//2016-6-30
					break;
				}
			} while (c[0] != 'E');
			break;
		}
	}

	if (pPRain->nhor == 0)
		pPRain->nhor = 35;

	rewind(pFile);
	while (!feof(pFile))
	{
		err_scanf = fscanf(pFile, "%c%c%*[^\n]%*c", &c[0], &c[1]);
		if (c[0] == 'I' && c[1] == 'P')
		{ //interaction parameters
			do
			{
				err_scanf = fscanf(pFile, "%c%*[^\n]%*c", &c[0]);
				if (c[0] == 'A')
				{
					err_scanf = fscanf(pFile, "%lf", &pPInter->interc);
					err_scanf = fscanf(pFile, "%d%lf%lf%lf", &pPInter->permtype, &pPInter->pers, &pPInter->perl, &pPInter->perm);
					err_scanf = fscanf(pFile, "%d%lf%lf%d", &pSSInter->type, &pSSInter->myu, &pSSInter->zeta, &pSSInter->acc_type);
					break;
				}
			} while (c[0] != 'E');
			break;
		}
	}

	rewind(pFile);
	while (!feof(pFile))
	{
		err_scanf = fscanf(pFile, "%c%c%*[^\n]%*c", &c[0], &c[1]);
		if (c[0] == 'B' && c[1] == 'P')
		{ //boundary parameters
			do
			{
				err_scanf = fscanf(pFile, "%c%*[^\n]%*c", &c[0]);
				if (c[0] == 'A')
				{
					err_scanf = fscanf(pFile, "%d%d%lf", &pPBou->if_cal, &pPBou->if_move, &pPBou->coe_bndy);
					if (pPBou->if_move >= 1)
						err_scanf = fscanf(pFile, "%lf%lf%lf", &pPBou->move_vx[0], &pPBou->move_vx[1], &pPBou->move_vx[2]);
					break;
				}
			} while (c[0] != 'E');
			break;
		}
	}

	rewind(pFile);
	while (!feof(pFile))
	{
		err_scanf = fscanf(pFile, "%c%c%*[^\n]%*c", &c[0], &c[1]);
		if (c[0] == 'V' && c[1] == 'L' && pFlag->flv >= 1)
		{ //vibration parameters
			do
			{
				err_scanf = fscanf(pFile, "%c%*[^\n]%*c", &c[0]);
				if (c[0] == 'A')
				{
					err_scanf = fscanf(pFile, "%lf%d", &pPVib->ttime, &pPVib->nsteps);
					err_scanf = fscanf(pFile, "%lf%lf", &pPVib->sttime, &pPVib->edtime);

					//initializing of vibrations
					if (pPVib->nsteps >= 1)
					{
						pPVib->memloads();
					}
					else
					{
						return 7; //vibrations steps error
					}

					break;
				}
			} while (c[0] != 'E');
			break;
		}
	}
	
	pPPro->tt = pPPro->dt * pPPro->loop;

	time(&rawtime);
	timeinfo = localtime(&rawtime);
	printf("Current time: %s", asctime(timeinfo));
	printf("The problem parameters have been initialized.\n");
	printf("-------------------------------------------------------------------\n");

	fprintf(flog, "Current time: %s", asctime(timeinfo));
	fprintf(flog, "The problem parameters have been initialized.\n");
	fprintf(flog, "-------------------------------------------------------------------\n");

	return 0;
}

//initialization of vibration loading from vibra.dat
int clFIni_Fun::initial_vib(Para_Vibra *pPVib, FILE *flog, FILE *fvib)
{
	int i;
	time_t rawtime;
	struct tm *timeinfo;

	time(&rawtime);
	timeinfo = localtime(&rawtime);
	printf("Current time: %s", asctime(timeinfo));
	printf("The vibration loadings are initializing\n.......\n");
	fprintf(flog, "Current time: %s", asctime(timeinfo));
	fprintf(flog, "The vibration loadings are initializing\n.......\n");

	for (i = 0; i < pPVib->nsteps; i++)
	{
		err_scanf = fscanf(fvib, "%lf%lf%lf", &(pPVib->loads[i][0]), &(pPVib->loads[i][1]), &(pPVib->loads[i][2]));
	}

	time(&rawtime);
	timeinfo = localtime(&rawtime);
	printf("Current time: %s", asctime(timeinfo));
	printf("The vibration loadings have been initialized.\n");
	printf("-------------------------------------------------------------------\n");

	fprintf(flog, "Current time: %s", asctime(timeinfo));
	fprintf(flog, "The vibration loadings have been initialized.\n");
	fprintf(flog, "-------------------------------------------------------------------\n");

	return 0;
}

//initializing from input file for material properties
int clFIni_Fun::initial_mat(Para_Model *pMatter, Para_Soil *pSoil, Para_Fluid *pWater,
							Para_Fluid *pAir, Para_Structure *pStruct, const Para_FL &pFlag, Flag_SpatialVariability *pFSV,
							Para_SpatialVariability *pPSV, FILE *flog, FILE *finp)
{
	FILE *pFile;
	char c[2];
	double maxp;
	time_t rawtime;
	struct tm *timeinfo;
	int n_layer, n_count, index, i;
	double damp_coe = 0.0;

	time(&rawtime);
	timeinfo = localtime(&rawtime);
	printf("Current time: %s", asctime(timeinfo));
	printf("The material parameters are initializing\n.......\n");

	fprintf(flog, "Current time: %s", asctime(timeinfo));
	fprintf(flog, "The material parameters are initializing\n.......\n");

	pFile = finp;

	rewind(pFile);
	while (!feof(pFile))
	{
		err_scanf = fscanf(pFile, "%c%c %*[^\n]%*c", &c[0], &c[1]);
		if (c[0] == 'S' && c[1] == 'F' && pFlag.flspv >= 1)
		{ //Spatial Variable Flags
			do
			{
				err_scanf = fscanf(pFile, "%c%*[^\n]%*c", &c[0]);
				if (c[0] == 'A')
				{
					err_scanf = fscanf(pFile, "%d %d %d %d", &pFSV->flag_fai, &pFSV->flag_c, &pFSV->flag_cop, &pFSV->flag_ds);
					err_scanf = fscanf(pFile, "%d %d", &pFSV->flag_faiw, &pFSV->flag_cw);
					err_scanf = fscanf(pFile, "%lf %lf %d %lf %lf", &pFSV->cell_ndx, &pFSV->coe_frict_cohesion, &pFSV->flag_auto, 
						&pFSV->beta, &pFSV->beta_1);
				}
			} while (c[0] != 'E');
			break;
		}
	}

	if (pFSV->flag_auto != 1 && pFSV->flag_auto != 2)
		pFSV->flag_auto = 1;  // default auto-correlation function is the single exponential funciton

	rewind(pFile);
	while (!feof(pFile))
	{
		err_scanf = fscanf(pFile, "%c%c %*[^\n]%*c", &c[0], &c[1]);
		if (c[0] == 'S' && c[1] == 'V' && pFlag.flspv >= 1)
		{ //Spatial Variable Parameters
			do
			{
				err_scanf = fscanf(pFile, "%c%*[^\n]%*c", &c[0]);
				if (c[0] == 'A')
				{

					//folid phase: fai, c, cop, ds
					err_scanf = fscanf(pFile, "%d%d%lf%lf%lf%lf%lf%lf", &pPSV[0].no, &pPSV[0].type,
									   &pPSV[0].mean_value, &pPSV[0].sd_value, &pPSV[0].rl_value[0], &pPSV[0].rl_value[1], &pPSV[0].rl_value[2], &pPSV[0].coe_depth);
					err_scanf = fscanf(pFile, "%d%d%lf%lf%lf%lf%lf%lf", &pPSV[1].no, &pPSV[1].type,
									   &pPSV[1].mean_value, &pPSV[1].sd_value, &pPSV[1].rl_value[0], &pPSV[1].rl_value[1], &pPSV[1].rl_value[2], &pPSV[1].coe_depth);
					err_scanf = fscanf(pFile, "%d%d%lf%lf%lf%lf%lf%lf", &pPSV[2].no, &pPSV[2].type, 
									   &pPSV[2].mean_value, &pPSV[2].sd_value, &pPSV[2].rl_value[0], &pPSV[2].rl_value[1], &pPSV[2].rl_value[2], &pPSV[2].coe_depth);
					err_scanf = fscanf(pFile, "%d%d%lf%lf%lf%lf%lf%lf", &pPSV[3].no, &pPSV[3].type, 
									   &pPSV[3].mean_value, &pPSV[3].sd_value, &pPSV[3].rl_value[0], &pPSV[3].rl_value[1], &pPSV[3].rl_value[2], &pPSV[3].coe_depth);

					//fluid phase: faiw, cw
					err_scanf = fscanf(pFile, "%d%d%lf%lf%lf%lf%lf%lf", &pPSV[4].no, &pPSV[4].type, 
									   &pPSV[4].mean_value, &pPSV[4].sd_value, &pPSV[4].rl_value[0], &pPSV[4].rl_value[1], &pPSV[4].rl_value[2], &pPSV[4].coe_depth);
					err_scanf = fscanf(pFile, "%d%d%lf%lf%lf%lf%lf%lf", &pPSV[5].no, &pPSV[5].type, 
									   &pPSV[5].mean_value, &pPSV[5].sd_value, &pPSV[5].rl_value[0], &pPSV[5].rl_value[1], &pPSV[5].rl_value[2], &pPSV[5].coe_depth);
				}
			} while (c[0] != 'E');
			break;
		}
	}

	rewind(pFile);
	while (!feof(pFile))
	{
		err_scanf = fscanf(pFile, "%c%c %*[^\n]%*c", &c[0], &c[1]);
		if (c[0] == 'W' && c[1] == 'P')
		{ //water phase
			do
			{
				err_scanf = fscanf(pFile, "%c%*[^\n]%*c", &c[0]);
				if (c[0] == 'A')
				{
					//add the input of density
					//different materials in fluid
					err_scanf = fscanf(pFile, "%lf%lf%lf%lf%lf%d%d%d%lf%lf", &pWater[0].dens, &pWater[0].vis, &pWater[0].fai, &pWater[0].c, &pWater[0].vpu,
									   &pWater[0].press_type, &pWater[0].acc_type, &pWater[0].bhtype, &pWater[0].pre_coe, &pWater[0].strain_min);
					err_scanf = fscanf(pFile, "%lf%lf%lf%lf%lf", &pWater[1].dens, &pWater[1].vis, &pWater[1].fai, &pWater[1].c, &pWater[1].vpu);
					err_scanf = fscanf(pFile, "%lf%lf%lf%lf%lf", &pWater[2].dens, &pWater[2].vis, &pWater[2].fai, &pWater[2].c, &pWater[2].vpu);
					break;
				}
			} while (c[0] != 'E');
			break;
		}
	}

	//different materials in fluid
	pWater[1].press_type = pWater[0].press_type;
	pWater[2].press_type = pWater[0].press_type;

	pWater[1].bhtype = pWater[0].bhtype;
	pWater[2].bhtype = pWater[0].bhtype;

	pWater[1].pre_coe = pWater[0].pre_coe;
	pWater[2].pre_coe = pWater[0].pre_coe;

	pWater[1].strain_min = pWater[0].strain_min;
	pWater[2].strain_min = pWater[0].strain_min;

	rewind(pFile);
	while (!feof(pFile))
	{
		err_scanf = fscanf(pFile, "%c%c %*[^\n]%*c", &c[0], &c[1]);
		if (c[0] == 'A' && c[1] == 'P' && pFlag.fla == 1)
		{ //air phase
			do
			{
				err_scanf = fscanf(pFile, "%c%*[^\n]%*c", &c[0]);
				if (c[0] == 'A')
				{
					//add the input of density
					err_scanf = fscanf(pFile, "%lf%lf%lf%lf%d%d%lf", &pAir->dens, &pAir->paa0, &pAir->vis,
									   &pAir->vpu, &pAir->press_type, &pAir->acc_type, &pAir->pre_coe);
					break;
				}
			} while (c[0] != 'E');
			break;
		}
	}

	rewind(pFile);
	while (!feof(pFile))
	{
		err_scanf = fscanf(pFile, "%c%c%*[^\n]%*c", &c[0], &c[1]);
		if (c[0] == 'M' && c[1] == 'S' && pFlag.fls == 1)
		{ //flags for soil model
			do
			{
				err_scanf = fscanf(pFile, "%c%*[^\n]%*c", &c[0]);
				if (c[0] == 'A')
				{
					//adding the elastic and perfect plastic model by M-C criterion
					err_scanf = fscanf(pFile, "%d%d%d%d%d%d%d%d%d%d", &pMatter->flagem, &pMatter->flagsb, &pMatter->flagdp, &pMatter->flagum,
									   &pMatter->flagam, &pMatter->flagbui, &pMatter->flagpl, &pMatter->flagve, &pMatter->flagfocc, &pMatter->flagep);
					//flagtsam 2017-7-8
					if (pFlag.flts >= 1)
						err_scanf = fscanf(pFile, "%lf%lf", &pMatter->alphats, &damp_coe);
					break;
				}
			} while (c[0] != 'E');
			break;
		}
	}

	rewind(pFile);
	while (!feof(pFile))
	{
		err_scanf = fscanf(pFile, "%c%c%*[^\n]%*c", &c[0], &c[1]);
		if (c[0] == 'E' && c[1] == 'M' && pFlag.fls == 1 && pMatter->flagem >= 1)
		{ //elastic
			do
			{
				err_scanf = fscanf(pFile, "%c%*[^\n]%*c", &c[0]);
				if (c[0] == 'A')
				{
					n_layer = pMatter->flagem;
					n_count = 0;
					do
					{
						index = n_count * 10 + 0;
						err_scanf = fscanf(pFile, "%lf%lf%lf", &pSoil[index].kapa, &pSoil[index].poi, &pSoil[index].dens);
						err_scanf = fscanf(pFile, "%lf%lf%lf%lf", &pSoil[index].damp1, &pSoil[index].damp2, &pSoil[index].cop, &pSoil[index].eps0);
						err_scanf = fscanf(pFile, "%lf%lf", &pSoil[index].e, &pSoil[index].ds);
						n_count++;
					} while (n_layer > n_count);
					break;
				}
			} while (c[0] != 'E');
			break;
		}
	}

	rewind(pFile);
	while (!feof(pFile))
	{
		err_scanf = fscanf(pFile, "%c%c%*[^\n]%*c", &c[0], &c[1]);
		if (c[0] == 'S' && c[1] == 'B' && pFlag.fls == 1 && pMatter->flagsb >= 1)
		{ //sub loading cam-clay model
			do
			{
				err_scanf = fscanf(pFile, "%c%*[^\n]%*c", &c[0]);
				if (c[0] == 'A')
				{
					n_layer = pMatter->flagsb;
					n_count = 0;
					do
					{
						index = n_count * 10 + 1;
						err_scanf = fscanf(pFile, "%lf%lf%lf%lf%lf%lf", &pSoil[index].poi, &pSoil[index].dens, &pSoil[index].eps0, &pSoil[index].zmf, &pSoil[index].sme0, &pSoil[index].fs);
						err_scanf = fscanf(pFile, "%lf%lf%lf%lf%lf", &pSoil[index].ramda, &pSoil[index].kapa, &pSoil[index].damp1, &pSoil[index].damp2, &pSoil[index].cop);
						err_scanf = fscanf(pFile, "%lf%lf", &pSoil[index].ann, &pSoil[index].ocr);
						err_scanf = fscanf(pFile, "%lf%lf", &pSoil[index].e, &pSoil[index].ds);
						n_count++;
					} while (n_layer > n_count);
					break;
				}
			} while (c[0] != 'E');
			break;
		}
	}

	rewind(pFile);
	while (!feof(pFile))
	{
		err_scanf = fscanf(pFile, "%c%c%*[^\n]%*c", &c[0], &c[1]);
		if (c[0] == 'D' && c[1] == 'P' && pFlag.fls == 1 && pMatter->flagdp >= 1)
		{ //drucker prager model of associated flow rule
			do
			{
				err_scanf = fscanf(pFile, "%c%*[^\n]%*c", &c[0]);
				if (c[0] == 'A')
				{
					n_layer = pMatter->flagdp;
					n_count = 0;
					do
					{
						index = n_count * 10 + 2;
						err_scanf = fscanf(pFile, "%lf%lf%lf%lf%lf%lf", &pSoil[index].kapa, &pSoil[index].poi, &pSoil[index].dens, &pSoil[index].c, &pSoil[index].fai, &pSoil[index].fs);
						err_scanf = fscanf(pFile, "%lf%lf%lf%lf", &pSoil[index].damp1, &pSoil[index].damp2, &pSoil[index].cop, &pSoil[index].eps0);
						err_scanf = fscanf(pFile, "%lf%lf", &pSoil[index].e, &pSoil[index].ds);
						n_count++;
					} while (n_layer > n_count);
					break;
				}
			} while (c[0] != 'E');
			break;
		}
	}

	rewind(pFile);
	while (!feof(pFile))
	{
		err_scanf = fscanf(pFile, "%c%c%*[^\n]%*c", &c[0], &c[1]);
		if (c[0] == 'U' && c[1] == 'M' && pFlag.fls == 1 && pMatter->flagum >= 1)
		{ //unsaturated cam-clay model
			do
			{
				err_scanf = fscanf(pFile, "%c%*[^\n]%*c", &c[0]);
				if (c[0] == 'A')
				{
					n_layer = pMatter->flagum;
					n_count = 0;
					do
					{
						index = n_count * 10 + 3;
						err_scanf = fscanf(pFile, "%lf%lf%lf%lf%lf%lf", &pSoil[index].poi, &pSoil[index].dens, &pSoil[index].eps0, &pSoil[index].zmf, &pSoil[index].sme0, &pSoil[index].fs);
						err_scanf = fscanf(pFile, "%lf%lf%lf%lf%lf", &pSoil[index].ramda, &pSoil[index].kapa, &pSoil[index].damp1, &pSoil[index].damp2, &pSoil[index].cop);
						err_scanf = fscanf(pFile, "%lf%lf%lf%lf%lf%lf%lf", &pSoil[index].ann, &pSoil[index].bnn, &pSoil[index].beta,
										   &pSoil[index].epss, &pSoil[index].epsr, &pSoil[index].sitas, &pSoil[index].sitar);
						err_scanf = fscanf(pFile, "%lf%lf%lf%lf%lf%lf", &pSoil[index].sd, &pSoil[index].sw, &pSoil[index].kse,
										   &pSoil[index].c1, &pSoil[index].c2, &pSoil[index].c3);
						err_scanf = fscanf(pFile, "%lf%lf", &pSoil[index].e, &pSoil[index].ds);
						n_count++;
					} while (n_layer > n_count);
					break;
				}
			} while (c[0] != 'E');
			break;
		}
	}

	rewind(pFile);
	while (!feof(pFile))
	{
		err_scanf = fscanf(pFile, "%c%c%*[^\n]%*c", &c[0], &c[1]);
		if (c[0] == 'A' && c[1] == 'M' && pFlag.fls == 1 && pMatter->flagam >= 1)
		{
			//anisotropic cam-clay model and revised at 2017-7-8
			do
			{
				err_scanf = fscanf(pFile, "%c%*[^\n]%*c", &c[0]);
				if (c[0] == 'A')
				{
					n_layer = pMatter->flagam;
					n_count = 0;
					do
					{
						index = n_count * 10 + 4;
						err_scanf = fscanf(pFile, "%lf%lf%lf%lf%lf%lf", &pSoil[index].poi, &pSoil[index].dens, &pSoil[index].eps0, &pSoil[index].zmf, &pSoil[index].sme0, &pSoil[index].fs);
						err_scanf = fscanf(pFile, "%lf%lf%lf%lf%lf", &pSoil[index].ramda, &pSoil[index].kapa, &pSoil[index].damp1, &pSoil[index].damp2, &pSoil[index].cop);
						err_scanf = fscanf(pFile, "%lf%lf%lf%lf%lf%lf", &pSoil[index].mzm, &pSoil[index].azm, &pSoil[index].r0, &pSoil[index].bl, &pSoil[index].bz, &pSoil[index].zet);
						err_scanf = fscanf(pFile, "%lf%lf", &pSoil[index].e, &pSoil[index].ds);
						n_count++;
					} while (n_layer > n_count);
					break;
				}
			} while (c[0] != 'E');
			break;
		}
	}

	rewind(pFile);
	while (!feof(pFile))
	{
		err_scanf = fscanf(pFile, "%c%c%*[^\n]%*c", &c[0], &c[1]);
		if (c[0] == 'B' && c[1] == 'U' && pFlag.fls == 1 && pMatter->flagbui >= 1)
		{ //drucker prager model proposed by Bui H.H.
			do
			{
				err_scanf = fscanf(pFile, "%c%*[^\n]%*c", &c[0]);
				if (c[0] == 'A')
				{
					n_layer = pMatter->flagbui;
					n_count = 0;
					do
					{
						index = n_count * 10 + 5;
						err_scanf = fscanf(pFile, "%lf%lf%lf%lf%lf%lf%lf", &pSoil[index].kapa,
										   &pSoil[index].poi, &pSoil[index].dens, &pSoil[index].c, &pSoil[index].fai, &pSoil[index].beta, &pSoil[index].fs);
						err_scanf = fscanf(pFile, "%lf%lf%lf%lf", &pSoil[index].damp1, &pSoil[index].damp2, &pSoil[index].cop, &pSoil[index].eps0);
						err_scanf = fscanf(pFile, "%lf%lf", &pSoil[index].e, &pSoil[index].ds);
						n_count++;
					} while (n_layer > n_count);
					break;
				}
			} while (c[0] != 'E');
			break;
		}
	}

	rewind(pFile);
	while (!feof(pFile))
	{
		err_scanf = fscanf(pFile, "%c%c%*[^\n]%*c", &c[0], &c[1]);
		if (c[0] == 'P' && c[1] == 'P' && pFlag.fls == 1 && pMatter->flagpl >= 1)
		{ //elastic model for pile
			do
			{
				err_scanf = fscanf(pFile, "%c%*[^\n]%*c", &c[0]);
				if (c[0] == 'A')
				{
					n_layer = pMatter->flagpl;
					n_count = 0;
					do
					{
						index = n_count * 10 + 6;
						err_scanf = fscanf(pFile, "%lf%lf%lf", &pSoil[index].kapa, &pSoil[index].poi, &pSoil[index].dens);
						err_scanf = fscanf(pFile, "%lf%lf%lf%lf", &pSoil[index].damp1, &pSoil[index].damp2, &pSoil[index].cop, &pSoil[index].eps0);
						err_scanf = fscanf(pFile, "%lf%lf", &pSoil[index].e, &pSoil[index].ds);
						n_count++;
					} while (n_layer > n_count);
					break;
				}
			} while (c[0] != 'E');
			break;
		}
	}

	rewind(pFile);
	while (!feof(pFile))
	{
		err_scanf = fscanf(pFile, "%c%c%*[^\n]%*c", &c[0], &c[1]);
		if (c[0] == 'V' && c[1] == 'E' && pFlag.fls == 1 && pMatter->flagve >= 1)
		{ //viscous elastic model for soil
			do
			{
				err_scanf = fscanf(pFile, "%c%*[^\n]%*c", &c[0]);
				if (c[0] == 'A')
				{
					n_layer = pMatter->flagve;
					n_count = 0;
					do
					{
						index = n_count * 10 + 7;
						err_scanf = fscanf(pFile, "%lf%lf%lf%lf", &pSoil[index].kapa, &pSoil[index].poi, &pSoil[index].dens, &pSoil[index].c);
						err_scanf = fscanf(pFile, "%lf%lf%lf%lf", &pSoil[index].damp1, &pSoil[index].damp2, &pSoil[index].cop, &pSoil[index].eps0);
						err_scanf = fscanf(pFile, "%lf%lf", &pSoil[index].e, &pSoil[index].ds);
						n_count++;
					} while (n_layer > n_count);
					break;
				}
			} while (c[0] != 'E');
			break;
		}
	}

	rewind(pFile);
	while (!feof(pFile))
	{
		err_scanf = fscanf(pFile, "%c%c%*[^\n]%*c", &c[0], &c[1]);
		if (c[0] == 'F' && c[1] == 'C' && pFlag.fls == 1 && pMatter->flagfocc >= 1)
		{ //fractional order model for MCC
			do
			{
				err_scanf = fscanf(pFile, "%c%*[^\n]%*c", &c[0]);
				if (c[0] == 'A')
				{
					n_layer = pMatter->flagfocc;
					n_count = 0;
					do
					{
						index = n_count * 10 + 8;
						err_scanf = fscanf(pFile, "%lf%lf%lf%lf%lf", &pSoil[index].dens, &pSoil[index].c, &pSoil[index].eps0, &pSoil[index].cop, &pSoil[index].sme0);
						err_scanf = fscanf(pFile, "%lf%lf%lf%lf%lf", &pSoil[index].g0, &pSoil[index].poi, &pSoil[index].zmf, &pSoil[index].ramda, &pSoil[index].etao);
						err_scanf = fscanf(pFile, "%lf%lf%lf%lf%lf", &pSoil[index].zet, &pSoil[index].ann, &pSoil[index].bnn, &pSoil[index].damp1, &pSoil[index].damp2);
						err_scanf = fscanf(pFile, "%lf", &pSoil[index].ds);
						pSoil[index].e = 5.0e6;
						n_count++;
					} while (n_layer > n_count);
					break;
				}
			} while (c[0] != 'E');
			break;
		}
	}

	rewind(pFile);
	//structure material parameter for Soil-Structure interaction calculation
	while (!feof(pFile))
	{
		err_scanf = fscanf(pFile, "%c%c%*[^\n]%*c", &c[0], &c[1]);
		if (c[0] == 'S' && c[1] == 'M')
		{ //elastic model for structure material
			do
			{
				err_scanf = fscanf(pFile, "%c%*[^\n]%*c", &c[0]);
				if (c[0] == 'A')
				{
					err_scanf = fscanf(pFile, "%lf%lf%lf", &pStruct->dens, &pStruct->eps0, &pStruct->cop);
					err_scanf = fscanf(pFile, "%lf%lf", &pStruct->e, &pStruct->poi);
					break;
				}
			} while (c[0] != 'E');
			break;
		}
	}

	rewind(pFile);
	while (!feof(pFile))
	{
		err_scanf = fscanf(pFile, "%c%c%*[^\n]%*c", &c[0], &c[1]);
		if (c[0] == 'E' && c[1] == 'P' && pFlag.fls == 1 && pMatter->flagep >= 1)
		{ //M-C failure criteria
			do
			{
				err_scanf = fscanf(pFile, "%c%*[^\n]%*c", &c[0]);
				if (c[0] == 'A')
				{
					n_layer = pMatter->flagep;
					n_count = 0;
					do
					{
						index = n_count * 10 + 9;
						err_scanf = fscanf(pFile, "%lf%lf%lf%lf%lf%lf", &pSoil[index].kapa, &pSoil[index].poi, &pSoil[index].dens, &pSoil[index].c, &pSoil[index].fai, &pSoil[index].fs);
						err_scanf = fscanf(pFile, "%lf%lf%lf%lf", &pSoil[index].damp1, &pSoil[index].damp2, &pSoil[index].cop, &pSoil[index].eps0);
						err_scanf = fscanf(pFile, "%lf%lf", &pSoil[index].e, &pSoil[index].ds);
						n_count++;
					} while (n_layer > n_count);
					break;
				}
			} while (c[0] != 'E');
			break;
		}
	}

	for (i = 0; i < 30; i++)
	{
		pSoil[i].vpu = 0600.0;
		pSoil[i].damp_coe = damp_coe;
	}
	pStruct->vpu = 0600.0;

	if (pFlag.fls == 0)
	{
		pAir->porosity = 1.0;
		pWater[0].porosity = 1.0;
		pWater[1].porosity = 1.0;
		pWater[2].porosity = 1.0;
	}
	else
	{

		for (i = 0; i < 30; i++)
		{
			pSoil[i].porosity = pSoil[i].eps0 / (1 + pSoil[i].eps0);
		}

		maxp = 0.0;
		for (i = 0; i < 30; i++)
		{
			maxp = fmax(maxp, pSoil[i].porosity);
		}

		pAir->porosity = maxp;
		pWater[0].porosity = maxp;
		pWater[1].porosity = maxp;
		pWater[2].porosity = maxp;
	}
	pStruct->porosity = pStruct->eps0 / (1 + pStruct->eps0);

	if (pWater[0].dens <= 0.0)
		pWater[0].dens = 1000.0;
	if (pWater[1].dens <= 0.0)
		pWater[1].dens = 1000.0;
	if (pWater[2].dens <= 0.0)
		pWater[2].dens = 1000.0;
	if (pAir->dens <= 0.0)
		pAir->dens = 1.293;

	time(&rawtime);
	timeinfo = localtime(&rawtime);
	printf("Current time: %s", asctime(timeinfo));
	printf("The material parameters have been initialized.\n");
	printf("-------------------------------------------------------------------\n");

	fprintf(flog, "Current time: %s", asctime(timeinfo));
	fprintf(flog, "The material parameters have been initialized.\n");
	fprintf(flog, "-------------------------------------------------------------------\n");

	return 0;
}

//setting particle information from file
void clFIni_Fun::input_par(Particle *pPar, Para_Pro *pPPro, Cell_Con *pCellc, Para_Rain *pPRain, FILE *flog, FILE *finp)
{
	double bmax[3], bmin[3];
	double max_nb[3], min_nb[3];
	double water_min[3], water_max[3]; /*problem domain of water phase*/
	double soil_min[3], soil_max[3]; /*problem domain of soil phase*/
	double air_min[3], air_max[3]; /*problem domain of air phase*/
	double stru_min[3], stru_max[3]; /*problem domain of structure phase*/
	int i, j, dim;
	FILE *pFile;
	char c[2];
	time_t rawtime;
	struct tm *timeinfo;
	int nrain = 0;

	dim = pPPro->ndim;

	time(&rawtime);
	timeinfo = localtime(&rawtime);
	printf("Current time: %s", asctime(timeinfo));
	printf("Particle information are initializing\n.......\n");

	fprintf(flog, "Current time: %s", asctime(timeinfo));
	fprintf(flog, "Particle information are initializing\n.......\n");

	pFile = finp;
	rewind(pFile);

	//file inputting
	while (!feof(pFile))
	{
		err_scanf = fscanf(pFile, "%c%c%*[^\n]%*c", &c[0], &c[1]);
		if (c[0] == 'P' && c[1] == 'A')
		{
			do
			{
				err_scanf = fscanf(pFile, "%c%*[^\n]%*c", &c[0]);
				if (c[0] == 'A')
				{
					if (pPPro->ndim == 2)
					{
						for (i = 0; i < pPPro->ntotal; i++)

							err_scanf = fscanf(pFile, "%lf%lf%lf%lf%d%d%d%d",
											   &pPar[i].xp[0], &pPar[i].xp[1], &pPar[i].vxp[0], &pPar[i].vxp[1],
											   &pPar[i].type, &pPar[i].matype, &pPar[i].permtype, &pPar[i].etype);
					}
					else
					{
						for (i = 0; i < pPPro->ntotal; i++)

							err_scanf = fscanf(pFile, "%lf%lf%lf%lf%lf%lf%d%d%d%d",
											   &pPar[i].xp[0], &pPar[i].xp[1], &pPar[i].xp[2], &pPar[i].vxp[0], &pPar[i].vxp[1],
											   &pPar[i].vxp[2], &pPar[i].type, &pPar[i].matype, &pPar[i].permtype, &pPar[i].etype);
					}
					break;
				}
			} while (c[0] != 'E');
			break;
		}
	}

	for (i = 0; i < 3; i++)
	{
		bmax[i] = -10000.0;
		bmin[i] = 10000.0;
		max_nb[i] = -10000.0;
		min_nb[i] = 10000.0;
		water_min[i] = 100000.000;
		water_max[i] = -100000.000;
		soil_min[i] = 100000.000;
		soil_max[i] = -100000.000;
		air_min[i] = 100000.000;
		air_max[i] = -100000.000;
		stru_min[i] = 100000.000;
		stru_max[i] = -100000.000;
	}

	for (i = 0; i < pPPro->ntotal; i++)
	{
		//search the domain for all particles
		for (j = 0; j < 3; j++)
		{
			bmax[j] = fmax(bmax[j], pPar[i].xp[j]);
			bmin[j] = fmin(bmin[j], pPar[i].xp[j]);
		}
		//search the domain for non boundary particles
		if (pPar[i].type != 0 && pPar[i].type != 7)
		{
			for (j = 0; j < 3; j++)
			{
				max_nb[j] = fmax(max_nb[j], pPar[i].xp[j]);
				min_nb[j] = fmin(min_nb[j], pPar[i].xp[j]);
			}
			//for each phase
			if (pPar[i].type == 1) {
				for (j = 0; j < 3; j++)
				{
					water_min[j] = fmin(water_min[j], pPar[i].xp[j]);
					water_max[j] = fmax(water_min[j], pPar[i].xp[j]);
				}
			}
			else if (pPar[i].type == 2) {
				for (j = 0; j < 3; j++)
				{
					soil_min[j] = fmin(soil_min[j], pPar[i].xp[j]);
					soil_max[j] = fmax(soil_min[j], pPar[i].xp[j]);
				}
			}
			else if (pPar[i].type == 3) {
				for (j = 0; j < 3; j++)
				{
					air_min[j] = fmin(air_min[j], pPar[i].xp[j]);
					air_max[j] = fmax(air_min[j], pPar[i].xp[j]);
				}
			}
			else if (pPar[i].type == 4) {
				for (j = 0; j < 3; j++)
				{
					stru_min[j] = fmin(stru_min[j], pPar[i].xp[j]);
					stru_max[j] = fmax(stru_min[j], pPar[i].xp[j]);
				}
			}
		}
		//search the number of rain particles
		if (pPar[i].type == 7 && pPar[i].etype == 1)
			nrain++;
	}

	for (j = 0; j < dim; j++)
	{
		bmax[j] = bmax[j] + 2.0 * pPPro->dr;
		bmin[j] = bmin[j] - 2.0 * pPPro->dr;
	}

	pCellc->xmax = bmax[0];
	pCellc->xmin = bmin[0];
	pCellc->ymax = bmax[1];
	pCellc->ymin = bmin[1];
	pCellc->zmax = bmax[2];
	pCellc->zmin = bmin[2];

	pCellc->xmax_nb = max_nb[0];
	pCellc->xmin_nb = min_nb[0];
	pCellc->ymax_nb = max_nb[1];
	pCellc->ymin_nb = min_nb[1];
	pCellc->zmax_nb = max_nb[2];
	pCellc->zmin_nb = min_nb[2];

	//reflecting data to pCellc variable
	for (j = 0; j < 3; j++)
	{
		pCellc->water_min[j] = water_min[j];
		pCellc->water_max[j] = water_max[j];
		pCellc->soil_min[j] = soil_min[j];
		pCellc->soil_max[j] = soil_max[j];
		pCellc->air_min[j] = air_min[j];
		pCellc->air_max[j] = air_max[j];
		pCellc->stru_min[j] = stru_min[j];
		pCellc->stru_max[j] = stru_max[j];
	}

	pPRain->nrain = nrain;
}

//setting viscosity and constitutive model parameter for each particle
void clFIni_Fun::setting_vis_cons(Particle *pPar, Para_Pro &pPPro, Para_Fluid *pVFluid, Para_Fluid *pWater, Para_Fluid &pAir,
								  Para_Soil *pParti_ConsPara, Para_Soil *pSoil, Para_Structure &pStruct)
{
	int i, index, index_1;

	for (i = 0; i < pPPro.ntotal; i++)
	{
		if (pPar[i].type == 1)
		{
			index = pPar[i].matype;
			pVFluid[i] = pWater[index];
		}
		if (pPar[i].type == 2)
		{
			index = pPar[i].matype;
			if (index < 30)
				pParti_ConsPara[i] = pSoil[index];
			else if (index >= 50) {
				index_1 = index % 10;
				pParti_ConsPara[i] = pSoil[index_1];
			}	
		}
		if (pPar[i].type == 3)
		{
			pVFluid[i] = pAir;
		}
		if (pPar[i].type == 7 && pPar[i].etype == 1)
		{
			index = pPar[i].matype;
			pVFluid[i] = pWater[index];
		}
		if (pPar[i].type == 7 && pPar[i].etype == 2)
		{
			index = pPar[i].matype;
			pParti_ConsPara[i] = pSoil[index];
		}
	}
}

//initializing from input file for particles
void clFIni_Fun::setting_par(Particle *pPar, Para_Pro *pPPro, Para_Fluid *pVFluid, Para_Soil *pParti_ConsPara,
							 Para_Soil *pSoil, Para_GF pPGf, const Para_Structure &pStruct, Para_Fluid *pWater, Para_Rain *pPRain, FILE *flog)
{

	double hmaxw, hmaxs, hminw, hmins;
	double porosity, sitar, sitas, dens, e0, pre_coe;
	int i, dim;
	time_t rawtime;
	struct tm *timeinfo;
	double g; //gravity acceleration

	dim = pPPro->ndim;
	hmaxw = -100000.0;
	hmaxs = -100000.0;
	hminw = 100000.0;
	hmins = 100000.0;
	porosity = 0.0;
	sitas = 0.0;
	sitar = 0.0;

	for (i = 0; i < 30; i++)
	{
		sitas = fmax(sitas, pSoil[i].sitas);
		sitar = fmax(sitar, pSoil[i].sitar);
	}

	/*initial setting*/
	for (i = 0; i < pPPro->ntotal; i++)
	{
		/*initial setting of parameters*/
		if (pPar[i].type == 0) /*boundary*/
		{
			pPar[i].rhop = 1500.0;
			pPar[i].satu = sitas;
			pPar[i].hl = pPPro->dr;
			pPar[i].massini(dim);
			pPPro->nboun = pPPro->nboun + 1;
			//initialization of free field particles gravititional acceleration
			if (pPar[i].matype == 5)
			{
				pPar[i].axp[0] = pPGf.gx;
				pPar[i].axp[1] = pPGf.gy;
				pPar[i].axp[2] = pPGf.gz;
			}
		}
		else if (pPar[i].type == 1) /*water*/
		{
			pPar[i].rhop = pVFluid[i].porosity * pVFluid[i].dens;
			pPar[i].satu = sitas;
			pPar[i].hl = pPPro->dr;
			pPar[i].massini(dim);
			pPPro->nwater = pPPro->nwater + 1;
			pPar[i].axp[0] = pPGf.gx;
			pPar[i].axp[1] = pPGf.gy;
			pPar[i].axp[2] = pPGf.gz;
			pPar[i].porosity = pVFluid[i].porosity;

			hmaxw = fmax(hmaxw, pPar[i].xp[dim - 1]);
			hminw = fmin(hminw, pPar[i].xp[dim - 1]);
		}
		else if (pPar[i].type == 7 && pPar[i].etype == 1) /*water*/
		{
			pPar[i].rhop = pVFluid[i].porosity * pVFluid[i].dens;
			pPar[i].satu = sitas;
			pPar[i].hl = pPPro->dr;
			pPar[i].massini(dim);
			pPar[i].axp[0] = pPGf.gx;
			pPar[i].axp[1] = pPGf.gy;
			pPar[i].axp[2] = pPGf.gz;
			pPar[i].porosity = pVFluid[i].porosity;
		}
		else if (pPar[i].type == 2) /*soil*/
		{

			porosity = pParti_ConsPara[i].porosity;
			e0 = pParti_ConsPara[i].eps0;
			pPar[i].rhop = (1 - porosity) * pParti_ConsPara[i].dens;

			if (pPar[i].matype == 0 || pPar[i].matype == 10 || pPar[i].matype == 20 || pPar[i].matype == 50)
			{
				pPar[i].satu = 1.0;
			}
			else if (pPar[i].matype == 1 || pPar[i].matype == 11 || pPar[i].matype == 21 || pPar[i].matype == 51)
			{
				pPar[i].satu = 1.0;
				pPar[i].roue = (pParti_ConsPara[i].ramda - pParti_ConsPara[i].kapa) * log(pParti_ConsPara[i].ocr);
			}
			else if (pPar[i].matype == 2 || pPar[i].matype == 12 || pPar[i].matype == 22 || pPar[i].matype == 52)
			{
				pPar[i].satu = 1.0;
			}
			else if (pPar[i].matype == 3 || pPar[i].matype == 13 || pPar[i].matype == 23 || pPar[i].matype == 53)
			{
				pPar[i].satu = pParti_ConsPara[i].sitar;
			}
			else if (pPar[i].matype == 4 || pPar[i].matype == 14 || pPar[i].matype == 24 || pPar[i].matype == 54)
			{
				pPar[i].satu = 1.0;
			}
			else if (pPar[i].matype == 5 || pPar[i].matype == 15 || pPar[i].matype == 25 || pPar[i].matype == 55)
			{
				pPar[i].satu = 1.0;
			}
			else if (pPar[i].matype == 6 || pPar[i].matype == 16 || pPar[i].matype == 26 || pPar[i].matype == 56)
			{
				pPar[i].satu = 1.0;
			}
			else if (pPar[i].matype == 7 || pPar[i].matype == 17 || pPar[i].matype == 27 || pPar[i].matype == 57)
			{
				pPar[i].satu = 1.0;
			}
			else if (pPar[i].matype == 8 || pPar[i].matype == 18 || pPar[i].matype == 28 || pPar[i].matype == 58)
			{ //focc
				pPar[i].satu = 1.0;
			}
			else if (pPar[i].matype == 9 || pPar[i].matype == 19 || pPar[i].matype == 29 || pPar[i].matype == 59)
			{ //MC
				pPar[i].satu = 1.0;
			}

			//2017-7-8
			hmaxs = fmax(hmaxs, pPar[i].xp[dim - 1]);
			hmins = fmin(hmins, pPar[i].xp[dim - 1]);
			pPar[i].e = e0;
			pPar[i].porosity = porosity;
			pPar[i].hl = pPPro->dr;
			pPar[i].massini(dim);
			pPPro->nsoil = pPPro->nsoil + 1;
			pPar[i].axp[0] = pPGf.gx;
			pPar[i].axp[1] = pPGf.gy;
			pPar[i].axp[2] = pPGf.gz;
			pPar[i].cd = pParti_ConsPara[i].damp_coe * sqrt(pParti_ConsPara[i].e / pPar[i].rhop / pPar[i].hl / pPar[i].hl);
		}
		else if (pPar[i].type == 7 && pPar[i].etype == 2) /*soil*/
		{

			porosity = pParti_ConsPara[i].porosity;
			e0 = pParti_ConsPara[i].eps0;
			pPar[i].rhop = (1 - porosity) * pParti_ConsPara[i].dens;

			if (pPar[i].matype == 0 || pPar[i].matype == 10 || pPar[i].matype == 20 || pPar[i].matype == 50)
			{
				pPar[i].satu = 1.0;
			}
			else if (pPar[i].matype == 1 || pPar[i].matype == 11 || pPar[i].matype == 21 || pPar[i].matype == 51)
			{
				pPar[i].satu = 1.0;
				pPar[i].roue = (pParti_ConsPara[i].ramda - pParti_ConsPara[i].kapa) * log(pParti_ConsPara[i].ocr);
			}
			else if (pPar[i].matype == 2 || pPar[i].matype == 12 || pPar[i].matype == 22 || pPar[i].matype == 52)
			{
				pPar[i].satu = 1.0;
			}
			else if (pPar[i].matype == 3 || pPar[i].matype == 13 || pPar[i].matype == 23 || pPar[i].matype == 53)
			{
				pPar[i].satu = pParti_ConsPara[i].sitar;
			}
			else if (pPar[i].matype == 4 || pPar[i].matype == 14 || pPar[i].matype == 24 || pPar[i].matype == 54)
			{
				pPar[i].satu = 1.0;
			}
			else if (pPar[i].matype == 5 || pPar[i].matype == 15 || pPar[i].matype == 25 || pPar[i].matype == 55)
			{
				pPar[i].satu = 1.0;
			}
			else if (pPar[i].matype == 6 || pPar[i].matype == 16 || pPar[i].matype == 26 || pPar[i].matype == 56)
			{
				pPar[i].satu = 1.0;
			}
			else if (pPar[i].matype == 7 || pPar[i].matype == 17 || pPar[i].matype == 27 || pPar[i].matype == 57)
			{
				pPar[i].satu = 1.0;
			}
			else if (pPar[i].matype == 8 || pPar[i].matype == 18 || pPar[i].matype == 28 || pPar[i].matype == 58)
			{ //focc
				pPar[i].satu = 1.0;
			}
			else if (pPar[i].matype == 9 || pPar[i].matype == 19 || pPar[i].matype == 29 || pPar[i].matype == 59)
			{ //MC
				pPar[i].satu = 1.0;
			}

			pPar[i].axp[0] = pPGf.gx;
			pPar[i].axp[1] = pPGf.gy;
			pPar[i].axp[2] = pPGf.gz;

			pPar[i].cop = pParti_ConsPara[i].cop;
			pPar[i].e = e0;
			pPar[i].porosity = porosity;
			pPar[i].hl = pPPro->dr;
			pPar[i].massini(dim);
			pPar[i].cd = pParti_ConsPara[i].damp_coe * sqrt(pParti_ConsPara[i].e / pPar[i].rhop / pPar[i].hl / pPar[i].hl);
		}
		else if (pPar[i].type == 3) /*air*/
		{
			pPar[i].rhop = pVFluid[i].porosity * pVFluid[i].dens;
			pPar[i].satu = 0.0;
			pPar[i].hl = pPPro->dr;
			pPar[i].massini(dim);
			pPPro->nair = pPPro->nair + 1;

			pPar[i].axp[0] = pPGf.gx;
			pPar[i].axp[1] = pPGf.gy;
			pPar[i].axp[2] = pPGf.gz;
		}
		else if (pPar[i].type == 4) /*structure*/
		{
			porosity = pStruct.porosity;
			pPar[i].rhop = (1 - porosity) * pStruct.dens;
			pPar[i].satu = 1.0;
			pPar[i].e = pStruct.eps0;
			pPar[i].hl = pPPro->dr;
			pPar[i].massini(dim);
			pPPro->nstruct = pPPro->nstruct + 1;
			pPar[i].axp[0] = pPGf.gx;
			pPar[i].axp[1] = pPGf.gy;
			pPar[i].axp[2] = pPGf.gz;
			pPar[i].cd = pParti_ConsPara[i].damp_coe * sqrt(pStruct.e / pPar[i].rhop / pPar[i].hl / pPar[i].hl);
		}

		pPar[i].rho = pPar[i].rhop;
	}
	dens = (pWater[0].dens + pWater[1].dens + pWater[2].dens) * 0.333333333333333;

	//setting gravity acceleration
	if (dim == 2) g = pPGf.gy;
	else if (dim == 3) g = pPGf.gz;
	else g = 0.0;

	/*initial pressure and suction*/
	for (i = 0; i < pPPro->ntotal; i++)
	{
		if (pPar[i].type == 0)
		{ //initializing of solid boundary pressure
			pPar[i].prep = 1500.0 * (-g) * fabs(hmaxw - hminw) * 0.5;
			pPar[i].pre = pPar[i].prep;
		}
		if (pPar[i].type == 1)
		{ //initializing of liquid pressure
			pre_coe = pVFluid[i].pre_coe;
			if (pVFluid[i].press_type == 0)
				pPar[i].prep = pPar[i].rhop * (-g) * pPPro->dr * pre_coe;
			else
				pPar[i].prep = pPar[i].rhop * (-g) * fabs(pPar[i].xp[dim - 1] - hmaxw - pre_coe * pPPro->dr);
			pPar[i].pre = pPar[i].prep;

			if (pPar[i].etype == 7)
			{
				pPRain->x_max = fmax(pPRain->x_max, pPar[i].xp[dim - 2]);
				pPRain->x_min = fmin(pPRain->x_min, pPar[i].xp[dim - 2]);
			}
		}
		if (pPar[i].type == 7 && pPar[i].etype == 1)
		{ //initializing of liquid pressure
			pre_coe = pVFluid[i].pre_coe;
			pPar[i].prep = pPar[i].rhop * (-g) * pPPro->dr * pre_coe;
			pPar[i].pre = pPar[i].prep;
		}
		if (pPar[i].type == 2)
		{ //initializing of soil stress
			if (fabs(pParti_ConsPara[i].sme0) > 1.0e-3)
			{
				pPar[i].sig[0] = -pParti_ConsPara[i].sme0;
				pPar[i].sig[1] = -pParti_ConsPara[i].sme0;
				pPar[i].sig[2] = -pParti_ConsPara[i].sme0;
			}
			else {
				if (dim == 2) {
					/*pPar[i].sig[1] = pPar[i].rho * (-g) * (pPar[i].xp[1] - hmaxs);
					pPar[i].sig[0] = pParti_ConsPara[i].poi / (1.0 - pParti_ConsPara[i].poi) * pPar[i].sig[1];
					pPar[i].sig[2] = pParti_ConsPara[i].poi / (1.0 - pParti_ConsPara[i].poi) * pPar[i].sig[1];*/
					pPar[i].sig[0] = pPar[i].rho * g * (hmaxs - hmins) * 0.5;
					pPar[i].sig[1] = pPar[i].rho * g * (hmaxs - hmins) * 0.5;
					pPar[i].sig[2] = pPar[i].rho * g * (hmaxs - hmins) * 0.5;
				}
				else if(dim == 3) {
					/*pPar[i].sig[2] = pPar[i].rho * (-g) * (pPar[i].xp[2] - hmaxs);
					pPar[i].sig[0] = pParti_ConsPara[i].poi / (1.0 - pParti_ConsPara[i].poi) * pPar[i].sig[2];
					pPar[i].sig[1] = pParti_ConsPara[i].poi / (1.0 - pParti_ConsPara[i].poi) * pPar[i].sig[2];*/
					pPar[i].sig[0] = pPar[i].rho * g * (hmaxs - hmins) * 0.5;
					pPar[i].sig[1] = pPar[i].rho * g * (hmaxs - hmins) * 0.5;
					pPar[i].sig[2] = pPar[i].rho * g * (hmaxs - hmins) * 0.5;
				}

			}

			if (pPar[i].matype == 3 || pPar[i].matype == 13 || pPar[i].matype == 23)
				pPar[i].suct = moisturesuct(pPar[i].satu, -0.0001, pParti_ConsPara[i]);
		}
		if (pPar[i].type == 7 && pPar[i].etype == 2)
		{ //initializing of soil stress
			if (fabs(pParti_ConsPara[i].sme0) > 1.0e-3)
			{
				pPar[i].sig[0] = -pParti_ConsPara[i].sme0;
				pPar[i].sig[1] = -pParti_ConsPara[i].sme0;
				pPar[i].sig[2] = -pParti_ConsPara[i].sme0;
			}
			if (pPar[i].matype == 3 || pPar[i].matype == 13 || pPar[i].matype == 23)
				pPar[i].suct = moisturesuct(pPar[i].satu, -0.0001, pParti_ConsPara[i]);
		}
		if (pPar[i].type == 3)
		{ //initializing of air pressure
			if (pVFluid[i].press_type == 0)
				pPar[i].prep = pVFluid[i].paa0;
			else
				pPar[i].prep = pVFluid[i].dens * (-g) * fabs(pPar[i].xp[dim - 1] - hmaxw - pVFluid[i].pre_coe * pPPro->dr);
			pPar[i].pre = pPar[i].prep;
		}
	}

	time(&rawtime);
	timeinfo = localtime(&rawtime);
	printf("Current time: %s", asctime(timeinfo));
	printf("Particles information has been initialized.\n");
	printf("-------------------------------------------------------------------\n");

	fprintf(flog, "Current time: %s", asctime(timeinfo));
	fprintf(flog, "Particles information has been initialized.\n");
	fprintf(flog, "-------------------------------------------------------------------\n");
}