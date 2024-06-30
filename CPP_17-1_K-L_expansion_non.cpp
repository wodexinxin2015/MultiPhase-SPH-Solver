/************************************************************************************
Multi-Phase and Parallelized SPH program
--parallelized and optimized for SMP (Symmetric multiprocessing ) architecture.
--Copyright (c) 2018-2023, Weijie ZHANG and Jian JI, GeoHohai, Hohai University.
************************************************************************************/

#define _CRT_SECURE_NO_WARNINGS
#include <random>
#include <chrono>
#include <iostream>
#include <time.h>
#include <string.h>
#include "Header_Parameters.h"
#include "Class_particle.h"
#include "Class_ParticleProperty.h"
#include "Class_Problem.h"
#include "Header_Option.h"
#include "Class_spatial_variability.h"

#define klterm 200  //terms for K-L expansions
#define precision 1.0e-6  //the precision for the equation solving

//generating a random variable that corresponds with Standard Normal Distribution
inline void generate_SND(double* kexi, int n)
{
	// construct a trivial random generator engine from a time-based seed:
	unsigned seed = (unsigned)std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);

	std::normal_distribution<double> distribution(0.0, 1.0);

	for (int i = 0; i < n; i++)
		kexi[i] = distribution(generator);
}

//solve the value ramdai
inline double equation_ramdai(double dirac, double wi)
{
	double ramdai;
	ramdai = 4.0 * dirac / (4.0 + wi * wi * dirac * dirac);
	return ramdai;
}

// solve the faii
inline double solve_faii(double x, double a, double wi, int i)
{
	double alphai, faii;

	if ((i + 1) % 2 == 1)
	{ // odd
		alphai = 1.0 / sqrt(sin(2.0 * wi * a) / 2.0 / wi + a);
		faii = alphai * cos(wi * x);
	}
	else
	{ // even
		alphai = 1.0 / sqrt(-sin(2.0 * wi * a) / 2.0 / wi + a);
		faii = alphai * sin(wi * x);
	}
	return faii;
}

//generate the spatial-distributed variables
int clFSpatialVariability_Fun::SpatialVariables_Generate_2D(Particle* pPar, Para_Fluid* pVFluid, Para_Soil* pParti_ConsPara,
	const Para_Pro& pPPro, const Flag_SpatialVariability fsv, Para_SpatialVariability* psv, const Cell_Con Cellc)
{
	double xmin[3], xmax[3], x_trans[3];
	int i, j, order_id[200][3], flag[6];
	int err_type = 0;
	double uu[klterm], wi[klterm][3], eigval[klterm];
	double water_max, soil_max;
	water_max = Cellc.water_max[1];
	soil_max = Cellc.soil_max[1];

	//problem domain a
	xmin[0] = Cellc.xmin_nb - 0.2 * (Cellc.xmax_nb - Cellc.xmin_nb);
	xmax[0] = Cellc.xmax_nb + 0.2 * (Cellc.xmax_nb - Cellc.xmin_nb);
	xmin[1] = Cellc.ymin_nb - 0.2 * (Cellc.ymax_nb - Cellc.ymin_nb);
	xmax[1] = Cellc.ymax_nb + 0.2 * (Cellc.ymax_nb - Cellc.ymin_nb);
	xmin[2] = 0.0;
	xmax[2] = 0.0;

	for (i = 0; i < 2; i++)
	{
		x_trans[i] = (xmax[i] + xmin[i]) / 2.0;
		for (j = 0; j < 6; j++)
			psv[j].a[i] = (xmax[i] - xmin[i]) / 2.0;
	}

	//check if a is zero
	for (i = 0; i < 6; i++)
	{
		if (psv[i].a[0] == 0.0)
			return 13;
		if (psv[i].a[1] == 0.0)
			return 13;
	}

	//SPV flag transformation
	if (fsv.flag_fai == 1)
	{
		flag[0] = 1;
	}
	else flag[0] = 0;

	if (fsv.flag_c == 1)
	{
		flag[1] = 1;
	}
	else flag[1] = 0;

	if (fsv.flag_cop == 1)
	{
		flag[2] = 1;
	}
	else flag[2] = 0;

	if (fsv.flag_ds == 1)
	{
		flag[3] = 1;
	}
	else flag[3] = 0;

	if (fsv.flag_faiw == 1)
	{
		flag[4] = 1;
	}
	else flag[4] = 0;

	if (fsv.flag_cw == 1)
	{
		flag[5] = 1;
	}
	else flag[5] = 0;

	//SPV flag loop
	for (i = 0; i < 6; i++) {
		if (flag[i] == 1) {
			//generate standard normal distribution variables
			generate_SND(uu, klterm);

			//solve wi
			solve_wi(wi, psv[i].a, psv[i].rl_value, 2);

			if (err_type != 0)
				return err_type;

			//solve ramda
			err_type = solve_ramdai(eigval, wi, psv[i].rl_value, order_id, 2);
			if (err_type != 0)
				return err_type;

			//solve H
			err_type = solve_h_2D(pPar, pVFluid, pParti_ConsPara, pPPro, wi, eigval,uu, x_trans, order_id, psv[i], soil_max, water_max, i);
			if (err_type != 0)
				return err_type;
		}
	}

	//output random variables
	output_Para_SpatialVariability(fsv, psv);

	return 0;
}

//generate the spatial-distributed variables
int clFSpatialVariability_Fun::SpatialVariables_Generate_3D(Particle* pPar, Para_Fluid* pVFluid, Para_Soil* pParti_ConsPara,
	const Para_Pro& pPPro, const Flag_SpatialVariability fsv, Para_SpatialVariability* psv, const Cell_Con Cellc)
{
	double xmin[3], xmax[3], x_trans[3];
	int i, j, order_id[200][3], flag[6];
	int err_type = 0;
	double uu[klterm], wi[klterm][3], eigval[klterm];
	double water_max, soil_max;
	water_max = Cellc.water_max[2];
	soil_max = Cellc.soil_max[2];

	//problem domain a
	xmin[0] = Cellc.xmin_nb - 0.2 * (Cellc.xmax_nb - Cellc.xmin_nb);
	xmax[0] = Cellc.xmax_nb + 0.2 * (Cellc.xmax_nb - Cellc.xmin_nb);
	xmin[1] = Cellc.ymin_nb - 0.2 * (Cellc.ymax_nb - Cellc.ymin_nb);
	xmax[1] = Cellc.ymax_nb + 0.2 * (Cellc.ymax_nb - Cellc.ymin_nb);
	xmin[2] = Cellc.zmin_nb - 0.2 * (Cellc.zmax_nb - Cellc.zmin_nb);
	xmax[2] = Cellc.zmax_nb + 0.2 * (Cellc.zmax_nb - Cellc.zmin_nb);
	for (i = 0; i < 3; i++)
	{
		x_trans[i] = (xmax[i] + xmin[i]) / 2.0;
		for (j = 0; j < 6; j++)
			psv[j].a[i] = (xmax[i] - xmin[i]) / 2.0;
	}

	//check if a is zero
	for (i = 0; i < 6; i++)
	{
		if (psv[i].a[0] == 0.0)
			return 13;
		if (psv[i].a[1] == 0.0)
			return 13;
	}

	//SPV flag transformation
	if (fsv.flag_fai == 1)
	{
		flag[0] = 1;
	}
	else flag[0] = 0;

	if (fsv.flag_c == 1)
	{
		flag[1] = 1;
	}
	else flag[1] = 0;

	if (fsv.flag_cop == 1)
	{
		flag[2] = 1;
	}
	else flag[2] = 0;

	if (fsv.flag_ds == 1)
	{
		flag[3] = 1;
	}
	else flag[3] = 0;

	if (fsv.flag_faiw == 1)
	{
		flag[4] = 1;
	}
	else flag[4] = 0;

	if (fsv.flag_cw == 1)
	{
		flag[5] = 1;
	}
	else flag[5] = 0;

	//SPV flag loop
	for (i = 0; i < 6; i++) {
		if (flag[i] == 1) {
			//generate standard normal distribution variables
			generate_SND(uu, klterm);

			//solve wi
			solve_wi(wi, psv[i].a, psv[i].rl_value, 3);

			if (err_type != 0)
				return err_type;

			//solve ramda
			err_type = solve_ramdai(eigval, wi, psv[i].rl_value, order_id, 3);
			if (err_type != 0)
				return err_type;

			//solve H
			err_type = solve_h_3D(pPar, pVFluid, pParti_ConsPara, pPPro, wi, eigval,uu, x_trans, order_id, psv[i], soil_max, water_max, i);
			if (err_type != 0)
				return err_type;
		}
	}

	//output random variables
	output_Para_SpatialVariability(fsv, psv);

	return 0;
}

//solve the value wi
void clFSpatialVariability_Fun::solve_wi(double (*wi)[3] , double a_x[3], double dirac_r[3], int ndim)
{
	double x1, x2, xm;
	double f1, fm, a, dira;
	double wi_upper[klterm], wi_lower[klterm];

	//upper and lower boundary for the wi
	for (int d = 0; d < ndim; d++)
	{
		a = a_x[d];
		dira = dirac_r[d];

#pragma omp parallel for
		for (int i = 0; i < klterm; i++)
		{
			if ((i + 1) % 2 == 1)
			{ // odd
				wi_upper[i] = (i + 0.5) * pi / a;
				wi_lower[i] = i * pi / a;
			}
			else
			{ // even
				wi_upper[i] = (i + 1.0) * pi / a;
				wi_lower[i] = (i + 0.5) * pi / a;
			}
		}

		// sovel wi
		for (int i = 0; i < klterm; i++)
		{
			if ((i + 1) % 2 == 1)
			{ // odd
				x1 = wi_lower[i] + precision;
				x2 = wi_upper[i] - precision;
				xm = (x1 + x2) * 0.5;
				do
				{
					f1 = 2.0 / dira - x1 * tan(x1 * a);
					fm = 2.0 / dira - xm * tan(xm * a);

					if ((f1 * fm) < 0.0)
					{
						x1 = x1;
						x2 = xm;
						xm = (x1 + x2) * 0.5;
					}
					else
					{
						x1 = xm;
						x2 = x2;
						xm = (x1 + x2) * 0.5;
					}
				} while (fabs(fm) > precision);
			}
			else
			{ // even
				x1 = wi_lower[i] + precision;
				x2 = wi_upper[i] - precision;
				xm = (x1 + x2) / 2.0;
				do
				{
					f1 = 2.0 / dira * tan(x1 * a) + x1;
					fm = 2.0 / dira * tan(xm * a) + xm;

					if ((f1 * fm) < 0.0)
					{
						x1 = x1;
						x2 = xm;
						xm = (x1 + x2) * 0.5;
					}
					else
					{
						x1 = xm;
						x2 = x2;
						xm = (x1 + x2) * 0.5;
					}
				} while (fabs(fm) > precision);
			}
			wi[i][d] = xm;
		}
	}
}

//solve the ramda i
int clFSpatialVariability_Fun::solve_ramdai(double *eigval, double (*wi)[3], double dirac_r[3],
											   int (*order_id)[3], int ndim)
{
	double ramda[klterm][3];
	double (*eigval_2d)[klterm] = new double[klterm][klterm];
	double (*eigval_3d)[klterm][klterm] = new double[klterm][klterm][klterm];

	// initialization of orders
#pragma omp parallel for
	for (int i = 0; i < klterm; i++)
	{
		order_id[i][0] = i;
		order_id[i][1] = i;
		order_id[i][2] = i;
	}

	// calculating ramda i for each idx in klterm
	for (int d = 0; d < ndim; d++)
	{
#pragma omp parallel for
		for (int i = 0; i < klterm; i++)
		{
			ramda[i][d] = equation_ramdai(dirac_r[d], wi[i][d]);
		}
	}

	//sorting eigenvalues by descending order
	if (ndim == 2){
		for (int i = 0; i < klterm; i++)
		{
			for (int j = 0; j < klterm; j++)
			{
				eigval_2d[i][j] = ramda[i][0] * ramda[j][1];
			}
		}

		//sorting the eigval_2d and find the order
		for (int i = 0; i < klterm; i++)
		{
			eigval[i] = -1.0;
			for (int m = 0; m < klterm; m++)
			{
				for (int n = 0; n < klterm; n++)
				{
					if (eigval[i] < eigval_2d[m][n])
					{
						eigval[i] = eigval_2d[m][n];
						order_id[i][0] = m;
						order_id[i][1] = n;
					}
				}
			}
			eigval_2d[order_id[i][0]][order_id[i][1]] = -1.0;
		}
	}
	else if (ndim == 3){
		for (int i = 0; i < klterm; i++)
		{
			for (int j = 0; j < klterm; j++)
			{
				for (int k = 0; k < klterm; k++)
				{
					eigval_3d[i][j][k] = ramda[i][0] * ramda[j][1] * ramda[k][2];
				}
			}
		}

		//sorting the eigval_3d and find the order
		for (int i = 0; i < klterm; i++)
		{
			eigval[i] = -1.0;
			for (int m = 0; m < klterm; m++)
			{
				for (int n = 0; n < klterm; n++)
				{
					for (int k = 0; k < klterm; k++)
					{
						if (eigval[i] < eigval_3d[m][n][k])
						{
							eigval[i] = eigval_3d[m][n][k];
							order_id[i][0] = m;
							order_id[i][1] = n;
							order_id[i][2] = k;
						}
					}
				}
			}
			eigval_3d[order_id[i][0]][order_id[i][1]][order_id[i][2]] = -1.0;
		}
	}

	delete[]eigval_2d;
	delete[]eigval_3d;

	return 0;
}

//output the problem parameters for the spatial variability
void clFSpatialVariability_Fun::output_Para_SpatialVariability(Flag_SpatialVariability fsv,
	Para_SpatialVariability* psv)
{
	int i;
	FILE* fop;
	time_t rawtime;
	struct tm* timeinfo;

	//open and append data at the end of a file
	fop = fopen("Para_Spatial_Dist.txt", "a");

	//output the time information
	time(&rawtime);
	timeinfo = localtime(&rawtime);
	fprintf(fop, "-----------------------------------------\n");
	fprintf(fop, "Current time: %s", asctime(timeinfo));
	fprintf(fop, "-----------------------------------------\n");

	//write the paramters
	fprintf(fop, "--flag_fai--flag_c--flag_cop--flag_ds\
--flag_faiw--flag_cw--\n");
	fprintf(fop, "%4d %4d %4d %4d %4d %4d \n", fsv.flag_fai, fsv.flag_c, fsv.flag_cop, fsv.flag_ds,  fsv.flag_faiw, fsv.flag_cw);
	fprintf(fop, "-----------------------------------------\n");

	fprintf(fop, "----no----type----dd_type----mean----sd_value\
----lx----ly----lz----ax----ay----az----\n");
	for (i = 0; i < 6; i++)
	{
		fprintf(fop, "%2d %2d  %7.3e  %7.3e  %7.3e  %7.3e  %7.3e  %7.3e  %7.3e  %7.3e\n",
			i + 1, psv[i].type,
			psv[i].mean_value, psv[i].sd_value, psv[i].rl_value[0],
			psv[i].rl_value[1], psv[i].rl_value[2], psv[i].a[0],
			psv[i].a[1], psv[i].a[2]);
	}
	fprintf(fop, "-----------------------------------------\n");

	//close the file
	if (fop != NULL)
		fclose(fop);
}

// sum for the random field of 2D
int clFSpatialVariability_Fun::solve_h_2D(Particle* pPar, Para_Fluid* pVFluid, Para_Soil* pParti_ConsPara,
	const Para_Pro& pPPro, double (*wi)[3], double *eigval, double *uu, double x_trans[3], int (*order_id)[3], 
	Para_SpatialVariability psv, double soil_max, double water_max, int type)
{
	int i, j, idx, idy;
	int ntotal = pPPro.ntotal;
	double eigvect[klterm], fai[3], x_a[3];
	double temp_h, temp_mean, temp_sigma;
	double cof_dep;

	for (i = 0; i < ntotal; i++)
	{
		if (pPar[i].type != 0 && pPar[i].type != 7)
		{
			//x value
			x_a[0] = pPar[i].xp[0] - x_trans[0];
			x_a[1] = pPar[i].xp[1] - x_trans[1];

			//sovle faii
#pragma omp parallel for private(idx, idy, fai)
			for (j = 0; j < klterm; j++)
			{
				idx = order_id[j][0];
				idy = order_id[j][1];

				fai[0] = solve_faii(x_a[0], psv.a[0], wi[idx][0], idx);
				fai[1] = solve_faii(x_a[1], psv.a[1], wi[idy][1], idy);

				eigvect[j] = fai[0] * fai[1];
			}
			
			//solve H value
			if (psv.type == 1)
			{ //standard normal distribution
				//initial setting
				cof_dep = psv.coe_depth;
				if (pPar[i].type == 2)
					temp_mean = psv.mean_value + (soil_max - pPar[i].xp[1]) * cof_dep;
				else if (pPar[i].type == 1)
					temp_mean = psv.mean_value + (water_max - pPar[i].xp[1]) * cof_dep;
				else
					temp_mean = psv.mean_value;

				temp_sigma = psv.sd_value;

				//solve H
				temp_h = 0.0;
				for (j = 0; j < klterm; j++)
				{
					temp_h = temp_h + sqrt(eigval[j]) * eigvect[j] * uu[j];
				}

				//reflect the random field to particle information
				if (type == 0 && pPar[i].type == 2)
					pParti_ConsPara[i].fai = temp_mean + temp_h * temp_sigma;
				if (type == 1 && pPar[i].type == 2)
					pParti_ConsPara[i].c = temp_mean + temp_h * temp_sigma;
				if (type == 2 && pPar[i].type == 2)
					pParti_ConsPara[i].cop = temp_mean + temp_h * temp_sigma;
				if (type == 3 && pPar[i].type == 2)
					pParti_ConsPara[i].ds = temp_mean + temp_h * temp_sigma;

				if (type == 4 && pPar[i].type == 1)
					pVFluid[i].fai = temp_mean + temp_h * temp_sigma;
				if (type == 5 && pPar[i].type == 1)
					pVFluid[i].c = temp_mean + temp_h * temp_sigma;
			}
			else if (psv.type == 2)
			{ //log normal distribution
				//initial setting
				cof_dep = psv.coe_depth;
				if (pPar[i].type == 2)
					temp_mean = psv.mean_value + (soil_max - pPar[i].xp[1]) * cof_dep;
				else if (pPar[i].type == 1)
					temp_mean = psv.mean_value + (water_max - pPar[i].xp[1]) * cof_dep;
				else
					temp_mean = psv.mean_value;

				temp_sigma = sqrt(log(1.0 + (psv.sd_value / temp_mean) * (psv.sd_value / temp_mean)));
				temp_mean = log(temp_mean) - 0.5 * temp_sigma * temp_sigma;
				
				//solve H
				temp_h = 0.0;
				for (j = 0; j < klterm; j++)
				{
					temp_h = temp_h + sqrt(eigval[j]) * eigvect[j] * uu[j];
				}

				//reflect the random field to particle information
				if (type == 0 && pPar[i].type == 2)
					pParti_ConsPara[i].fai = exp(temp_mean + temp_h * temp_sigma);
				if (type == 1 && pPar[i].type == 2)
					pParti_ConsPara[i].c = exp(temp_mean + temp_h * temp_sigma);
				if (type == 2 && pPar[i].type == 2)
					pParti_ConsPara[i].cop = exp(temp_mean + temp_h * temp_sigma);
				if (type == 3 && pPar[i].type == 2)
					pParti_ConsPara[i].ds = exp(temp_mean + temp_h * temp_sigma);

				if (type == 4 && pPar[i].type == 1)
					pVFluid[i].fai = exp(temp_mean + temp_h * temp_sigma);
				if (type == 5 && pPar[i].type == 1)
					pVFluid[i].c = exp(temp_mean + temp_h * temp_sigma);
			}
			else
				return 14;
		}
	}

	return 0;
}

// sum for the random field of 3D
int clFSpatialVariability_Fun::solve_h_3D(Particle* pPar, Para_Fluid* pVFluid, Para_Soil* pParti_ConsPara,
	const Para_Pro& pPPro, double (*wi)[3], double *eigval, double *uu, double x_trans[3], int (*order_id)[3], 
	Para_SpatialVariability psv, double soil_max, double water_max, int type)
{
	int i, j, idx, idy, idz;
	int ntotal = pPPro.ntotal;
	double eigvect[klterm], fai[3], x_a[3];
	double temp_h, temp_mean, temp_sigma;
	double cof_dep;

	for (i = 0; i < ntotal; i++)
	{
		if (pPar[i].type != 0 && pPar[i].type != 7)
		{
			//x value
			x_a[0] = pPar[i].xp[0] - x_trans[0];
			x_a[1] = pPar[i].xp[1] - x_trans[1];
			x_a[2] = pPar[i].xp[2] - x_trans[2];

			//sovle faii
#pragma omp parallel for private(idx, idy, fai)
			for (j = 0; j < klterm; j++)
			{
				idx = order_id[j][0];
				idy = order_id[j][1];
				idz = order_id[j][2];

				fai[0] = solve_faii(x_a[0], psv.a[0], wi[idx][0], idx);
				fai[1] = solve_faii(x_a[1], psv.a[1], wi[idy][1], idy);
				fai[2] = solve_faii(x_a[2], psv.a[2], wi[idz][2], idz);

				eigvect[j] = fai[0] * fai[1] * fai[2];
			}
			
			//solve H value
			if (psv.type == 1)
			{ //standard normal distribution
				//initial setting
				cof_dep = psv.coe_depth;
				if (pPar[i].type == 2)
					temp_mean = psv.mean_value + (soil_max - pPar[i].xp[1]) * cof_dep;
				else if (pPar[i].type == 1)
					temp_mean = psv.mean_value + (water_max - pPar[i].xp[1]) * cof_dep;
				else
					temp_mean = psv.mean_value;

				temp_sigma = psv.sd_value;

				//solve H
				temp_h = 0.0;
				for (j = 0; j < klterm; j++)
				{
					temp_h = temp_h + sqrt(eigval[j]) * eigvect[j] * uu[j];
				}

				//reflect the random field to particle information
				if (type == 0 && pPar[i].type == 2)
					pParti_ConsPara[i].fai = temp_mean + temp_h * temp_sigma;
				if (type == 1 && pPar[i].type == 2)
					pParti_ConsPara[i].c = temp_mean + temp_h * temp_sigma;
				if (type == 2 && pPar[i].type == 2)
					pParti_ConsPara[i].cop = temp_mean + temp_h * temp_sigma;
				if (type == 3 && pPar[i].type == 2)
					pParti_ConsPara[i].ds = temp_mean + temp_h * temp_sigma;

				if (type == 4 && pPar[i].type == 1)
					pVFluid[i].fai = temp_mean + temp_h * temp_sigma;
				if (type == 5 && pPar[i].type == 1)
					pVFluid[i].c = temp_mean + temp_h * temp_sigma;
			}
			else if (psv.type == 2)
			{ //log normal distribution
				//initial setting
				cof_dep = psv.coe_depth;
				if (pPar[i].type == 2)
					temp_mean = psv.mean_value + (soil_max - pPar[i].xp[1]) * cof_dep;
				else if (pPar[i].type == 1)
					temp_mean = psv.mean_value + (water_max - pPar[i].xp[1]) * cof_dep;
				else
					temp_mean = psv.mean_value;

				temp_sigma = sqrt(log(1.0 + (psv.sd_value / temp_mean) * (psv.sd_value / temp_mean)));
				temp_mean = log(temp_mean) - 0.5 * temp_sigma * temp_sigma;
				
				//solve H
				temp_h = 0.0;
				for (j = 0; j < klterm; j++)
				{
					temp_h = temp_h + sqrt(eigval[j]) * eigvect[j] * uu[j];
				}

				//reflect the random field to particle information
				if (type == 0 && pPar[i].type == 2)
					pParti_ConsPara[i].fai = exp(temp_mean + temp_h * temp_sigma);
				if (type == 1 && pPar[i].type == 2)
					pParti_ConsPara[i].c = exp(temp_mean + temp_h * temp_sigma);
				if (type == 2 && pPar[i].type == 2)
					pParti_ConsPara[i].cop = exp(temp_mean + temp_h * temp_sigma);
				if (type == 3 && pPar[i].type == 2)
					pParti_ConsPara[i].ds = exp(temp_mean + temp_h * temp_sigma);

				if (type == 4 && pPar[i].type == 1)
					pVFluid[i].fai = exp(temp_mean + temp_h * temp_sigma);
				if (type == 5 && pPar[i].type == 1)
					pVFluid[i].c = exp(temp_mean + temp_h * temp_sigma);
			}
			else
				return 14;
		}
	}

	return 0;
}

int clFSpatialVariability_Fun::random_field_input(Particle* pPar, Para_Fluid* pVFluid, Para_Soil* pParti_ConsPara,
	const Para_Pro& pPPro, char* argv){
	char file_name[100];
	int i, err, id, type, matype;
	FILE* fpt;
	char c[2];
	double x[3], temp[2];

	//open file
	memset(file_name, 0, 100);
	strcpy(file_name, argv);
	if (win32)
		strcat(file_name, "\\Parameters");
	else
		strcat(file_name, "/Parameters");
	strcat(file_name, ".txt");
	fpt = fopen(file_name, "r");

	//input constitutive parameters from Parameters.txt
	err = fscanf(fpt, "%c%c%*[^\n]%*c", &c[0], &c[1]);
	if (fpt != NULL) {
		for (i = 0; i < pPPro.ntotal; i++)
		{
			if (pPar[i].type == 1) {   //for water particles
				err = fscanf(fpt, "%d%lf%lf%lf%lf%lf%lf%lf%d%d", &id, &x[0], &x[1], &x[2],
					&pVFluid[i].fai, &pVFluid[i].c, &temp[0], &temp[1], &type, &matype);
			}
			else if (pPar[i].type == 2) {   //for soil particles
				err = fscanf(fpt, "%d%lf%lf%lf%lf%lf%lf%lf%d%d", &id, &x[0], &x[1], &x[2],
					&pParti_ConsPara[i].fai, &pParti_ConsPara[i].c, &pParti_ConsPara[i].cop,
					&pParti_ConsPara[i].ds, &type, &matype);
			}
			else if (pPar[i].type == 3) {   //for air particles
				err = fscanf(fpt, "%d%lf%lf%lf%lf%lf%lf%lf%d%d", &id, &x[0], &x[1], &x[2],
					&pVFluid[i].fai, &pVFluid[i].c, &temp[0], &temp[1], &type, &matype);
			}
			else if (pPar[i].type == 4) {   //for structure particles
				err = fscanf(fpt, "%d%lf%lf%lf%lf%lf%lf%lf%d%d", &id, &x[0], &x[1], &x[2],
					&pParti_ConsPara[i].fai, &pParti_ConsPara[i].c, &pParti_ConsPara[i].cop,
					&pParti_ConsPara[i].ds, &type, &matype);
			}
			else {  //other particles: slip the input
				err = fscanf(fpt, "%c%c%*[^\n]%*c", &c[0], &c[1]);
			}
		}

		/*Close files*/
		if (fpt != NULL)
			fclose(fpt);

		return 0;
	}
	else return 1;
	
}