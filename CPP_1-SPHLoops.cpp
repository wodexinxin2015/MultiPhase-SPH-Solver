/************************************************************************************
Multi-Phase and Parallelized SPH program
--parallelized and optimized for SMP (Symmetric multiprocessing ) architecture.
--Copyright (c) 2018-2023, Weijie ZHANG, GeoHohai, Hohai University.
************************************************************************************/
#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <omp.h>
#include <time.h>
#include "string.h"
#include "Header_Option.h"
#include "Class_spatial_variability.h"
#include "Class_Functions.h"

int err_scanf;

void statement(FILE *flog)
{
	/*program information*/
	/*the version of program must be reflected in this function*/
	printf("-------------------------------------------------------------------\n");
	printf("Multi-Phase and Parallelized SPH program 5.1.1\n");
	printf("--Copyright (c) 2018-2023 Weijie ZHANG, GeoHohai, Hohai University.\n");
	printf("-------------------------------------------------------------------\n");
	printf("2023-03-20: V4.13.1.0\n");
	printf("	reduce the memory consumption of NNPS variables.\n");
	printf("2023-03-24: V4.13.2.0\n");
	printf("	add the correction algorithm of CSPM into the density, strain and XSPH calculation.\n");
	printf("2023-03-28: V4.13.3.0\n");
	printf("	revise the constitutive model of DP from DBLeaves, DP of Bui and Elastic-Plastic of MC;\n");
	printf("	revise the determiniation of time step increment.\n");
	printf("2023-04-02: V4.13.4.0\n");
	printf("	change the module of particle_into_cell to the Binary sorting method in NNPS.\n");
	printf("2023-06-19: V4.13.5.0\n");
	printf("	Revise the DP model and Bui-DP model.\n");
	printf("	Revise the segmentation fault bug in the CPP_3-NNPS.cpp.\n");
	printf("2023-06-26: V4.14.1.0b\n");
	printf("	in development.\n");
	printf("2023-08-18£ºV5.1.1\n");
	printf("	start a new version number.\n");
	printf("-------------------------------------------------------------------\n");
	fprintf(flog, "-------------------------------------------------------------------\n");
	fprintf(flog, "	Multi-Phase and Parallelized SPH program 5.1.1\n");
	fprintf(flog, "	--Copyright (c) 2018-2023 Weijie ZHANG, GeoHohai, Hohai University.\n");
	fprintf(flog, "-------------------------------------------------------------------\n");
	fprintf(flog, "2023-03-20: V4.13.1.0\n");
	fprintf(flog, "	reduce the memory consumption of NNPS variables.\n");
	fprintf(flog, "2023-03-24: V4.13.2.0\n");
	fprintf(flog, "	add the correction algorithm of CSPM into the density, strain and XSPH calculation.\n");
	fprintf(flog, "2023-03-24: V4.13.3.0\n");
	fprintf(flog, "	revise the constitutive model of DP from DBLeaves, DP of Bui and Elastic-Plastic of MC;\n");
	fprintf(flog, "	revise the determiniation of time step increment.\n");
	fprintf(flog, "2023-04-02: V4.13.4.0\n");
	fprintf(flog, "	change the module of particle_into_cell to the Binary sorting method in NNPS.\n");
	fprintf(flog, "2023-06-19: V4.13.5.0\n");
	fprintf(flog, "	Revise the DP model and Bui-DP model.\n");
	fprintf(flog, "	Revise the segmentation fault bug in the CPP_3-NNPS.cpp.\n");
	fprintf(flog, "2023-06-26: V4.14.1.0b\n");
	fprintf(flog, "	in development.\n");
	fprintf(flog, "2023-08-18£ºV5.1.1\n");
	fprintf(flog, "	start a new version number.\n");
	fprintf(flog, "-------------------------------------------------------------------\n");
}

int SPH_Calculation(int argc, char *argv)
{

	// file pointer for MGF file, time log, flowing distance, total impact force,
	// input and vibration input file
	FILE *fMgf, *ftime, *flog, *fFld, *fTotf, *finp, *fvib;
	int lp, i, err_t, st_count;
	double vrain, cement_base[3];
	time_t rawtime;
	struct tm *timeinfo;
	int nstage;
	int k, size_m; //variables for the resize of parti_cell, size_m=2^k.

	// for free surface check
	double threshold_min, threshold_max;
	double(*vx_temp)[3] = NULL;
	threshold_min = 0.1;
	threshold_max = 0.12;

	// Function instances
	clFIni_Fun ItFini;
	clErr_Fun ItErrD;
	clNNPS_Fun ItNNPS;
	clDensity_Fun ItDens;
	clStraStre_Fun ItStaSts;
	clInterFor_Fun ItInterF;
	clAcce_Fun ItAcce;
	clVel_Fun ItVel;
	clUpdate_Fun ItUpdate;
	clBndy_Fun ItBndy;
	clOutput_Fun ItOutput;
	clExcBac_Fun ItExcBac;
	clRain_Fun ItRain;
	clSoilStruct_Fun ItSSInter;
	clFSpatialVariability_Fun ItSpatialFun;
	clFParaOutput_Fun ItParaOutput;
	clFreeFieldBndy_Fun ItFreeFieldParti;
	clFreeFieldBndy_Fun_soil ItFreeFieldParti_soil;

	// working directory initializing
	char temp[100];
	memset(temp, 0, 100);

	// path for windows
	if (win32)
	{
		// file for log output
		strcpy(temp, argv);
		strcat(temp, "\\svlog.txt");
		if ((flog = fopen(temp, "w")) == NULL)
		{
			ItErrD.Err_Deal(10, flog);
			return 10;
		}
		memset(temp, 0, 100);
	}
	// path for Linux
	else
	{
		// file for log output
		strcpy(temp, argv);
		strcat(temp, "/svlog.txt");
		if ((flog = fopen(temp, "w")) == NULL)
		{
			ItErrD.Err_Deal(10, flog);
			return 10;
		}
		memset(temp, 0, 100);
	}
	/*statement*/
	statement(flog);

	// path for windows
	if (win32)
	{
		// file for flowing distance output
		strcpy(temp, argv);
		strcat(temp, "\\Flowing distance.dat");
		if ((fFld = fopen(temp, "w")) == NULL)
		{
			ItErrD.Err_Deal(11, flog);
			return 11;
		}
		memset(temp, 0, 100);
		// file for time step
		strcpy(temp, argv);
		strcat(temp, "\\step_time.log");
		if ((ftime = fopen(temp, "w")) == NULL)
		{
			ItErrD.Err_Deal(11, flog);
			return 11;
		}
		memset(temp, 0, 100);
		// Initializing of MGF file
		strcpy(temp, argv);
		strcat(temp, "\\AVS_MGF.mgf");
		if ((fMgf = fopen(temp, "w")) == NULL)
		{
			ItErrD.Err_Deal(11, flog);
			return 11;
		}
		memset(temp, 0, 100);
		// Initializing of total impact force file
		strcpy(temp, argv);
		strcat(temp, "\\Total_Force.dat");
		if ((fTotf = fopen(temp, "w")) == NULL)
		{
			ItErrD.Err_Deal(11, flog);
			return 11;
		}
		memset(temp, 0, 100);
		// Initializing of input file
		strcpy(temp, argv);
		strcat(temp, "\\input.dat");
		if ((finp = fopen(temp, "r")) == NULL)
		{
			ItErrD.Err_Deal(1, flog);
			return 1;
		}
		memset(temp, 0, 100);
	}
	// path for Linux
	else
	{
		// file for flowing distance output
		strcpy(temp, argv);
		strcat(temp, "/Flowing distance.dat");
		if ((fFld = fopen(temp, "w")) == NULL)
		{
			ItErrD.Err_Deal(11, flog);
			return 11;
		}
		memset(temp, 0, 100);
		// file for time step
		strcpy(temp, argv);
		strcat(temp, "/step_time.log");
		if ((ftime = fopen(temp, "w")) == NULL)
		{
			ItErrD.Err_Deal(11, flog);
			return 11;
		}
		memset(temp, 0, 100);
		// Initializing of MGF file
		strcpy(temp, argv);
		strcat(temp, "/AVS_MGF.mgf");
		if ((fMgf = fopen(temp, "w")) == NULL)
		{
			ItErrD.Err_Deal(11, flog);
			return 11;
		}
		memset(temp, 0, 100);
		// Initializing of total impact force file
		strcpy(temp, argv);
		strcat(temp, "/Total_Force.dat");
		if ((fTotf = fopen(temp, "w")) == NULL)
		{
			ItErrD.Err_Deal(11, flog);
			return 11;
		}
		memset(temp, 0, 100);
		// Initializing of input file
		strcpy(temp, argv);
		strcat(temp, "/input.dat");
		if ((finp = fopen(temp, "r")) == NULL)
		{
			ItErrD.Err_Deal(1, flog);
			return 1;
		}
		memset(temp, 0, 100);
	}

	// problem global parameters
	Para_FL *pFlag = new Para_FL();
	Para_Pro *pPPro = new Para_Pro();
	Para_Rain *pPRain = new Para_Rain();
	Para_Inter *pPInter = new Para_Inter();
	Para_Boundary *pPBou = new Para_Boundary();
	Para_Vibra *pPVib = new Para_Vibra();
	Para_GF *pPGf = new Para_GF();
	Para_EC *pPEc = new Para_EC();
	Para_SSInter *pSSInter = new Para_SSInter();
	Flag_SpatialVariability *pFSV = new Flag_SpatialVariability();

	// material property parameters
	Para_Model *pMatter = new Para_Model();
	Para_Soil *pSoil = new Para_Soil[30];
	Para_Soil *pFFBound = new Para_Soil();
	Para_Fluid *pWater = new Para_Fluid[3];
	Para_Fluid *pAir = new Para_Fluid();
	Para_Structure *pStruct = new Para_Structure();

	// parameters of spatial variability
	Para_SpatialVariability *pPSV = new Para_SpatialVariability[6];

	// NNPS variables
	Cell_Con *pCellc = new Cell_Con();
	cell_info *pCell_Info = NULL;
	cell_link *pCell_Link = NULL;
	Par_Cell *pParCell = NULL;
	parti_cellid *pParti_Cell_Sorted = NULL;
	parti_cellid *pParti_Cell_Temp = NULL;

	// variables for particles
	Particle *pPar = NULL;
	StiffMat *pParStiff = NULL;

	// variables for soil-structure interaction
	clPair_SS *pPair_SS = NULL;
	clInterf_SS **pInterf_SS = NULL;
	clParti_Pair *pPartiPair_SS = NULL;

	// variables for the random field simulation
	Para_Fluid *pVFluid = NULL;
	Para_Soil *pParti_ConsPara = NULL;

	// variable for boundary particle velocity
	cl_bndy_pair *pBndy_pair = NULL;

	// Raij for artificial stress calculatio
	clRaij *Rij = NULL;

	// variable for boundary treatment
	clVar_Boundary *pParti_VariBndy = NULL;

	// initializing from input file for problems and material property
	err_t = ItFini.initial_pro(pFlag, pPPro, pPRain, pPInter, pSSInter, pPBou, pPVib, pPGf, pPEc,
							   &threshold_min, &threshold_max, flog, finp);
	if (err_t > 0)
	{
		ItErrD.Err_Deal(err_t, flog);
		return err_t;
	}

	// update time increment
	if (pPPro->dttype == 2)
	{
		ItUpdate.timestep(pPPro);
	}

	// initialization of vibrations from file vibra.dat
	if (pFlag->flv >= 1)
	{
		if (win32)
		{
			strcpy(temp, argv);
			strcat(temp, "\\vibra.dat");
			if ((fvib = fopen(temp, "r")) == NULL)
			{
				ItErrD.Err_Deal(1, flog);
				return 1;
			}
			memset(temp, 0, 100);
		}
		else
		{
			strcpy(temp, argv);
			strcat(temp, "/vibra.dat");
			if ((fvib = fopen(temp, "r")) == NULL)
			{
				ItErrD.Err_Deal(1, flog);
				return 1;
			}
			memset(temp, 0, 100);
		}
		err_t = ItFini.initial_vib(pPVib, flog, fvib);
		if (err_t > 0)
		{
			ItErrD.Err_Deal(err_t, flog);
			return err_t;
		}

		if (fvib != NULL)
			fclose(fvib);
	}

	// setting openmp environment
	omp_set_num_threads(pFlag->ntd);

	time(&rawtime);
	timeinfo = localtime(&rawtime);
	printf("Current time: %s", asctime(timeinfo));
	printf("OpenMP environment is ready for calculating with %d threads.\n", pFlag->ntd);
	printf("-------------------------------------------------------------------\n");
	fprintf(flog, "Current time: %s", asctime(timeinfo));
	fprintf(flog, "OpenMP environment is ready for calculating with %d threads.\n", pFlag->ntd);
	fprintf(flog, "-------------------------------------------------------------------\n");

	// allocate memory
	k = (int)log2(pPPro->ntotal) + 1;
	size_m = 1 << k;
	if (pPPro->ntotal > 0)
	{
		pPar = new Particle[pPPro->ntotal];
		Rij = new clRaij[pPPro->ntotal];
		pParCell = new Par_Cell[pPPro->ntotal];
		pParStiff = new StiffMat[pPPro->ntotal];

		pInterf_SS = new clInterf_SS *[pPPro->ntotal];
		for (i = 0; i < pPPro->ntotal; i++)
		{
			pInterf_SS[i] = new clInterf_SS[pFlag->ntd];
		}
		pPair_SS = new clPair_SS[pPPro->ntotal];
		pPartiPair_SS = new clParti_Pair[pPPro->ntotal];
		pVFluid = new Para_Fluid[pPPro->ntotal];
		pParti_ConsPara = new Para_Soil[pPPro->ntotal];
		pBndy_pair = new cl_bndy_pair[pPPro->ntotal];
		pParti_VariBndy = new clVar_Boundary[pPPro->ntotal];
		pParti_Cell_Sorted = new parti_cellid[size_m]; // array
		pParti_Cell_Temp = new parti_cellid[size_m];   // temp array
		vx_temp = new double[pPPro->ntotal][3];
	}
	else
	{
		ItErrD.Err_Deal(2, flog);
		return 2;
	}

	// initialization of material property
	err_t = ItFini.initial_mat(pMatter, pSoil, pWater, pAir, pStruct, *pFlag, pFSV, pPSV, flog, finp);
	if (err_t > 0)
	{
		ItErrD.Err_Deal(err_t, flog);
		return err_t;
	}

	// initializing from input file for particle information
	ItFini.input_par(pPar, pPPro, pCellc, pPRain, flog, finp);

	// close file pointer finp
	fclose(finp);

	// setting constitutive model parameters for each particle
	ItFini.setting_vis_cons(pPar, *pPPro, pVFluid, pWater, *pAir, pParti_ConsPara, pSoil, *pStruct);

	// generation of random field variables
	if (pFlag->flspv == 1 && pPPro->ndim == 2) // K-L expansion of non-stationary random field
		err_t = ItSpatialFun.SpatialVariables_Generate_2D(pPar, pVFluid, pParti_ConsPara, *pPPro, *pFSV, pPSV, *pCellc);
	else if (pFlag->flspv == 1 && pPPro->ndim == 3)
		err_t = ItSpatialFun.SpatialVariables_Generate_3D(pPar, pVFluid, pParti_ConsPara, *pPPro, *pFSV, pPSV, *pCellc);

	if (pFlag->flspv == 2 && pPPro->ndim == 2) // K-L expansion of non-stationary random field with correlation of c and fai
		err_t = ItSpatialFun.SpatialVariables_Generate_2D_corr(pPar, pVFluid, pParti_ConsPara, *pPPro, *pFSV, pPSV, *pCellc);
	else if (pFlag->flspv == 2 && pPPro->ndim == 3)
		err_t = ItSpatialFun.SpatialVariables_Generate_3D_corr(pPar, pVFluid, pParti_ConsPara, *pPPro, *pFSV, pPSV, *pCellc);

	if (pFlag->flspv == 3 && pPPro->ndim == 2) // correlateion matrix method and nonstationary random field
		err_t = ItSpatialFun.nonstationary_random_2D(pPar, pVFluid, pParti_ConsPara, *pPPro, *pFSV, pPSV, *pCellc);
	else if (pFlag->flspv == 3 && pPPro->ndim == 3)
		err_t = ItSpatialFun.nonstationary_random_3D(pPar, pVFluid, pParti_ConsPara, *pPPro, *pFSV, pPSV, *pCellc);

	if (pFlag->flspv == 4 && pPPro->ndim == 2) // correlateion matrix method, nonstationary and correlated random field between friction angle and cohesion
		err_t = ItSpatialFun.nonstationary_correlated_random_2D(pPar, pVFluid, pParti_ConsPara, *pPPro, *pFSV, pPSV, *pCellc);
	else if (pFlag->flspv == 4 && pPPro->ndim == 3)
		err_t = ItSpatialFun.nonstationary_correlated_random_3D(pPar, pVFluid, pParti_ConsPara, *pPPro, *pFSV, pPSV, *pCellc);

	if (pFlag->flspv == 5)
	{ // 5--using the correlated random varible for the homogeneous field
		err_t = ItSpatialFun.random_variable_generate(pPar, pVFluid, pParti_ConsPara, *pPPro, *pFSV, pPSV);
	}

	if (pFlag->flspv == -1)
	{ // using the random field from the parameters.txt
		err_t = ItSpatialFun.random_field_input(pPar, pVFluid, pParti_ConsPara, *pPPro, argv);
	}

	if (err_t > 0)
	{
		ItErrD.Err_Deal(err_t, flog);
		return err_t;
	}

	// setting the initial density, pressure and stress for each particle
	ItFini.setting_par(pPar, pPPro, pVFluid, pParti_ConsPara, pSoil, *pPGf, *pStruct, pWater, pPRain, flog);

	// output for the parameters
	ItParaOutput.paramters_output(pPar, pVFluid, pParti_ConsPara, *pPPro, argv);

	// initializing for Cell information
	ItNNPS.celllist_ini1(pCellc, *pPPro);
	pCell_Info = new cell_info[pCellc->ctotal + 4];
	pCell_Link = new cell_link[pCellc->ctotal + 4];

	ItNNPS.celllist_ini2(pCell_Link, *pCellc, pPPro->ndim, pFlag->ntd);

	fprintf(fMgf, "# Micro AVS Geom:2.10\n# SPH simulation\n");
	fprintf(fMgf, "%d\n", (int)(pPPro->loop / pPPro->fop));

	// initializing of rainfall parameters
	if (pFlag->flr == 1 && pPRain->rtype <= 7 && pPRain->rtype >= 1)
	{
		ItRain.rainini(*pPRain, *pPPro, &lp, &vrain);
		ItRain.center_cement(pPar, *pPPro, cement_base);
	}
	else if (pFlag->flr == 1)
	{
		ItErrD.Err_Deal(6, flog);
		return 6;
	}

	time(&rawtime);
	timeinfo = localtime(&rawtime);
	printf("Current time: %s", asctime(timeinfo));
	printf("Main loops of SPH are running.\n");
	fprintf(flog, "Current time: %s", asctime(timeinfo));
	fprintf(flog, "Main loops of SPH are running.\n");

	st_count = 1;
	nstage = 0;

	// Main SPH loop
	for (pPPro->l = 1; pPPro->l <= pPPro->loop; pPPro->l++)
	{
		// print information on screen
		if (pPPro->t < pPPro->tt)
		{
			pPPro->t = pPPro->t + pPPro->dt;
			if (mod(pPPro->l, pPPro->fop) == 0)
			{
				st_count += 1;
				printf("It is %8d step of total %8d steps, %6.3e s for output.  \n", pPPro->l, pPPro->loop, pPPro->t);
				fprintf(ftime, "Time of %8d step is %6.3e s\n", pPPro->l, pPPro->t);
				fprintf(flog, "It is %8d step of total %8d steps, %6.3e s for output.  \n", pPPro->l, pPPro->loop, pPPro->t);
			}
			else
			{
				ItOutput.output_status(pPPro->l, (st_count - 1) * pPPro->fop, st_count * pPPro->fop);
			}

			// produce rainfall
			if (pFlag->flr == 1 && (pPPro->l - pPPro->inip) > 0 && mod(pPPro->l - pPPro->inip, lp) == 0)
			{
				// rainfall producing
				if (pPRain->rtype == 1)
					ItRain.rainpro(pPar, pPRain, pPPro, pVFluid, *pPGf, vrain);
				// simulation of decreasing water level in the left side
				else if (pPRain->rtype == 2)
					ItRain.waterdel_min(pPar, pPRain, pPPro);
				// simulation of decreasing water level in the right side
				else if (pPRain->rtype == 3)
					ItRain.waterdel_max(pPar, pPRain, pPPro);
				// cement injection particles producing
				else if (pPRain->rtype == 4)
					ItRain.cementpro(pPar, pPRain, pPPro, pVFluid, *pPGf, vrain, cement_base);
			}

			// excavation or backfill
			if (pFlag->fleob == 1)
			{
				if (pPPro->l == pPEc->ecp[nstage][2] && nstage < pPEc->ecount)
				{
					ItExcBac.exca_or_back(pPar, pPPro, *pPEc, pFlag->ntd, nstage);
					nstage += 1;
				}
			}

			// Nearest neighboring particles searching
			err_t = ItNNPS.parttocell(pParCell, pCell_Info, *pCellc, pParti_Cell_Sorted, pParti_Cell_Temp,
				*pPPro, pPar, pPPro->ndim, pFlag->ntd, size_m, k);
			if (err_t > 0)
			{
				ItErrD.Err_Debug(pPPro->l, -1, flog);
				ItErrD.Err_Deal(err_t, flog);
				// output to the inp file
				ItOutput.output_total(pPar, *pPPro, argv);
				return err_t;
			}
			err_t = ItNNPS.partisearching(pParCell, pCell_Info, pCell_Link, pParti_Cell_Sorted, *pCellc, *pPPro, pPar, pFlag->ntd, pFlag->flknl);
			if (err_t > 0)
			{
				ItErrD.Err_Debug(pPPro->l, -1, flog);
				ItErrD.Err_Deal(err_t, flog);
				// output to the inp file
				ItOutput.output_total(pPar, *pPPro, argv);
				return err_t;
			}

			// kernel gradient correction
			if (pFlag->flknl > 2)
				ItNNPS.kernel_gradient_correction(pPar, pParCell, *pPPro);

			// density for velocity particles and stress particles
			ItDens.density(pPar, pParCell, pVFluid, pBndy_pair, pParti_VariBndy, *pPPro, pFlag->flbndy);

			// strain rate tensor and spin rate tensor
			ItStaSts.strain_noreg(pPar, pParCell, pBndy_pair, pParti_VariBndy, *pPPro, pFlag->flbndy);

			// strain regularization
			if (mod(pPPro->l, pPPro->step_regu) == 0)
			{
				ItStaSts.strain_regu(pPar, pParCell, *pPPro);
			}

			// saturation
			if (pMatter->flagum == 1)
				ItInterF.saturation(pPar, pParCell, pParti_ConsPara, *pPPro, pFlag->ntd);

			// calculate permeability
			if (pFlag->fli == 1)
				ItInterF.cal_permeability(pPar, *pPInter, pParti_ConsPara, *pPPro, pFlag->ntd);

			// pressure calculation from equation of state
			ItStaSts.eos(pPar, pParCell, pVFluid, *pPPro, pFlag->ntd);

			// strain->stress using soil model
			ItStaSts.strain_to_stress(pPar, pParti_ConsPara, pParStiff, pVFluid, *pStruct, *pPPro, pFlag->ntd, pFlag->flts, pMatter->alphats);

			// stress regularization
			if (mod(pPPro->l, pPPro->step_regu) == 0)
			{
				ItStaSts.stress_regu(pPar, pParCell, *pPPro);
			}

			// free surface detection and stress correction
			ItStaSts.free_surface_check(pPar, pParCell, *pPPro, threshold_max, threshold_min);

			// principle stress calculation
			ItStaSts.Stress_eigValue(pPar, *pPPro, pFlag->ntd);

			// Accelerations
			if (pWater[0].acc_type == 0)
				ItAcce.fluid_acceleration_water(pPar, pParCell, pVFluid, *pPPro, *pPGf, pFlag->ntd);
			else
				ItAcce.fluid_acceleration_minorwater(pPar, pParCell, pVFluid, *pPPro, *pPGf, pFlag->ntd);

			if (pAir->acc_type == 0)
				ItAcce.fluid_acceleration_air(pPar, pParCell, pVFluid, *pPPro, *pPGf, pFlag->ntd);
			else
				ItAcce.fluid_acceleration_minorair(pPar, pParCell, pVFluid, *pPPro, *pPGf, pFlag->ntd);

			ItAcce.soil_acceleration_noreg(pPar, pParCell, pParti_ConsPara, *pPPro, *pPGf, pFlag->ntd);

			// acceleration of structure particles
			ItAcce.structure_acceleration_noreg(pPar, pParCell, *pPPro, *pPGf, pFlag->ntd);

			// artifical viscosity
			ItAcce.artificial_viscosity(pPar, pParCell, pParti_ConsPara, pVFluid, *pPPro, pFlag->ntd);

			// artifical stress
			if (pPPro->ndim == 2 && pPPro->art_stre_coe > 0.0)
				ItAcce.artificial_stress_2d(pPar, pParCell, *pPPro, Rij, pFlag->ntd);
			else if (pPPro->ndim == 3 && pPPro->art_stre_coe > 0.0)
				ItAcce.artificial_stress_3d(pPar, pParCell, *pPPro, Rij, pFlag->ntd);

			// boundary treatment of Bui HH and Chong Peng
			if (pFlag->flv == 0)
			{
				// action of free field particles on the moving particles
				if (pFlag->flbndy == 1 || pFlag->flbndy == 2)
					ItAcce.boundary_moving_stress_effect(pPar, pParCell, *pPPro, pFlag->ntd);
				else if (pFlag->flbndy == 3 || pFlag->flbndy == 4 || pFlag->flbndy == 13 || pFlag->flbndy == 14)
					ItAcce.boundary_stress_effect_tran(pPar, pParCell, pParti_VariBndy, *pPPro, pFlag->ntd);
				else
					ItAcce.boundary_effect_free_struct(pPar, pParCell, *pPPro, pFlag->ntd);
			}
			else
			{
				// action of free field particles on the moving particles
				if (pFlag->flbndy == 1 || pFlag->flbndy == 2)
					ItAcce.boundary_moving_stress_effect(pPar, pParCell, *pPPro, pFlag->ntd);
				else if (pFlag->flbndy == 3 || pFlag->flbndy == 4 || pFlag->flbndy == 13 || pFlag->flbndy == 14)
					ItAcce.boundary_stress_effect_tran(pPar, pParCell, pParti_VariBndy, *pPPro, pFlag->ntd);
				else
					ItAcce.boundary_effect_free_struct(pPar, pParCell, *pPPro, pFlag->ntd);

				ItFreeFieldParti.free_field_solid(pPar, pBndy_pair, *pFFBound, pParCell,
												  *pPPro, pFlag->ntd);
				ItFreeFieldParti_soil.free_field_solid(pPar, pParti_ConsPara, pParCell, *pPPro, pFlag->ntd);
			}

			// acceleration of soil, water and structure
			if (pSSInter->type >= 1 && pSSInter->type <= 3)
			{
				if (pSSInter->acc_type == 2)
					ItAcce.coupled_Water_Soil_Structure(pPar, pParCell, pVFluid, pParti_ConsPara, *pPPro, *pPGf, pFlag->ntd);
			}

			// vibration acceleration
			if (pFlag->flv == 1)
				ItAcce.vibration_acceleration(pPar, *pPVib, *pPPro, pFlag->ntd);

			// Interaction force
			if (pFlag->fli == 1)
			{ // interaction force of Darcy's law
				ItInterF.interact1(pPar, pParCell, *pPInter, pVFluid, *pPPro, *pPGf, pFlag->ntd);
			}
			else if (pFlag->fli == 2)
			{ // interaction force inlcuding laminar flow and turbulent flow from Tsuji
				ItInterF.interact2(pPar, pParCell, *pPInter, pParti_ConsPara, pVFluid, *pPPro, *pPGf, pFlag->ntd);
			}
			else if (pFlag->fli == 3)
			{ // interaction force inlcuding laminar flow and turbulent flow from Peng C. et al.
				ItInterF.interact3(pPar, pParCell, *pPInter, pParti_ConsPara, pVFluid, *pPPro, *pPGf, pFlag->ntd);
			}

			// Interaction force of soil and structure
			if (pSSInter->type == 1)
			{
				err_t = ItSSInter.interaction_ss_rigid(pPar, pParCell, pPair_SS, pInterf_SS, pPartiPair_SS,
													   *pSSInter, *pPPro, pFlag->ntd);
				if (err_t > 0)
				{
					ItErrD.Err_Deal(err_t, flog);
					return err_t;
				}
			}
			else if (pSSInter->type == 2)
			{
				err_t = ItSSInter.interaction_ss_deformable(pPar, pParCell, pPair_SS, pInterf_SS, pPartiPair_SS,
															*pSSInter, *pPPro, pFlag->ntd);
				if (err_t > 0)
				{
					ItErrD.Err_Deal(err_t, flog);
					return err_t;
				}
			}
			else if (pSSInter->type == 3)
			{
				err_t = ItSSInter.interaction_ss_penalty(pPar, pParCell, pPair_SS, pInterf_SS, pPartiPair_SS,
														 *pSSInter, *pPPro, pFlag->ntd);
				if (err_t > 0)
				{
					ItErrD.Err_Deal(err_t, flog);
					return err_t;
				}
			}

			// damping effect in initialization of geostress
			if (pPPro->l <= pPPro->inip)
				ItAcce.soil_damping_geostress(pPar, pParCell, *pPPro, pFlag->ntd);
			if (pFlag->flv > 0)
				ItAcce.soil_damping_vibration(pPar, pParCell, pParti_ConsPara, pParStiff, *pPPro);

			// boundary effect using the free particle-solid wall interaction
			if (pFlag->flbndy == 5)
			{
				ItBndy.boundary_free_solid(pPar, pParCell, *pPPro, *pSSInter, pFlag->ntd);
			}

			// velocity update
			ItVel.veocity_update(pPar, *pPPro, *pPBou, pFlag->flv);

			// vibration velocity
			if (pFlag->flv == 2)
				ItAcce.vibration_velocity(pPar, *pPVib, *pPPro, pFlag->ntd);

			// XSPH correction
			if (pPBou->if_move != 2)
			{
				ItVel.xsph(pPar, pParCell, pBndy_pair, pParti_VariBndy, *pPPro, vx_temp, pFlag->flbndy);
			}

			// Boundary condition
			if (pFlag->flbndy == 1)
			{
				ItBndy.tmboundary(pPar, pParCell, pBndy_pair, *pPPro);
			}
			else if (pFlag->flbndy == 2)
			{
				ItBndy.velocity_set_domain(pPar, pParCell, pBndy_pair, *pPPro, pPBou->coe_bndy);
				if (pPPro->ndim == 2)
					ItBndy.tm_boundary_2d(pPar, pParCell, pBndy_pair, pParti_VariBndy, *pPPro, pFlag->ntd);
				else if (pPPro->ndim == 3)
					ItBndy.tm_boundary_3d(pPar, pParCell, pBndy_pair, pParti_VariBndy, *pPPro, pFlag->ntd);
			}
			else if (pFlag->flbndy == 3)
			{
				ItBndy.no_slip_tran(pPar, pParCell, pParti_VariBndy, *pPPro, pFlag->ntd);
			}
			else if (pFlag->flbndy == 4)
			{
				ItBndy.free_slip_tran(pPar, pParCell, pParti_VariBndy, *pPPro, pFlag->ntd);
			}
			else if (pFlag->flbndy == 13)
			{
				ItBndy.velocity_set_domain(pPar, pParCell, pBndy_pair, *pPPro, pPBou->coe_bndy);
				ItBndy.no_slip_tran(pPar, pParCell, pParti_VariBndy, *pPPro, pFlag->ntd);
			}
			else if (pFlag->flbndy == 14)
			{
				ItBndy.velocity_set_domain(pPar, pParCell, pBndy_pair, *pPPro, pPBou->coe_bndy);
				ItBndy.free_slip_tran(pPar, pParCell, pParti_VariBndy, *pPPro, pFlag->ntd);
			}
			ItBndy.check_domain(pPar, *pCellc, *pPBou, *pPPro);

			// information update for all particles
			ItUpdate.update(pPar, *pPPro, *pPBou, pFlag->ntd);

			// particles delete
			if (pFlag->flr == 1)
			{
				ItRain.partidel(pPar, pPRain, pPPro, *pCellc, pFlag->ntd);
			}

			// output to file
			if (mod(pPPro->l, pPPro->fop) == 0 || pPPro->l == 1)
			{
				// output to the inp file
				ItOutput.output_total(pPar, *pPPro, argv);
				ItOutput.output_flowdist(fFld, pPar, *pPPro);
				// output the impact forces
				if (pPBou->if_cal == 1)
					ItOutput.output_impact(fTotf, pPar, *pPPro, argv);
				// output to the mgf file
				if (pPPro->l != 1)
					ItOutput.mgf(pPar, *pPPro, *pCellc, fMgf);
			}
		}
		else
		{
			// output to the inp file
			ItOutput.output_total(pPar, *pPPro, argv);
			ItOutput.output_flowdist(fFld, pPar, *pPPro);
			break;
		}
	}

	time(&rawtime);
	timeinfo = localtime(&rawtime);
	printf("Current time: %s", asctime(timeinfo));
	printf("Main loops of SPH have been finished.\n");
	printf("-------------------------------------------------------------------\n");

	fprintf(flog, "Current time: %s", asctime(timeinfo));
	fprintf(flog, "Main loops of SPH have been finished.\n");
	fprintf(flog, "-------------------------------------------------------------------\n");
	fclose(ftime);

	// Finishing MGF file
	fclose(fMgf);
	fclose(fFld);
	fclose(flog);
	fclose(fTotf);

	// delete memory
	delete pFlag;
	delete pPRain;
	delete pPInter;
	delete pPBou;
	delete pPVib;
	delete pPGf;
	delete pPEc;
	delete pSSInter;
	delete pMatter;
	delete[] pSoil;
	delete pFFBound;
	delete[] pWater;
	delete pAir;
	delete pCellc;
	delete[] pParCell;
	delete[] pCell_Info;
	delete[] pCell_Link;
	delete[] pPar;
	delete[] Rij;
	delete[] pParStiff;
	delete[] pBndy_pair;
	delete[] pParti_VariBndy;
	delete[] pParti_Cell_Temp;
	delete[] pParti_Cell_Sorted;

	delete pStruct;
	delete[] pPair_SS;
	delete[] pPartiPair_SS;
	for (i = 0; i < pPPro->ntotal; i++)
	{
		delete[] pInterf_SS[i];
	}
	delete[] pInterf_SS;
	delete pPPro;
	delete[] pVFluid;
	delete[] pParti_ConsPara;

	delete pFSV;
	delete[] pPSV;

	delete[] vx_temp;

	return 0;
}