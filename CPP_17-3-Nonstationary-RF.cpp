/************************************************************************************
Multi-Phase and Parallelized SPH program
--parallelized and optimized for SMP (Symmetric multiprocessing ) architecture.
--Copyright (c) 2018-2023, Weijie ZHANG, GeoHohai, Hohai University.
************************************************************************************/

#define _CRT_SECURE_NO_WARNINGS
#include <cstring>
#include <random>
#include <chrono>
#include <iostream>
#include <time.h>
#include <cmath>
#include "Header_Parameters.h"
#include "Class_particle.h"
#include "Class_ParticleProperty.h"
#include "Class_Problem.h"
#include "Class_spatial_variability.h"

inline void matrix_multi(double* mat_inp, double* vec_inp, double* vec_out, int n) {
	int i, j;
	double temp;
#pragma omp parallel for schedule(static) private(j,temp)
	for (i = 0; i < n; i++) {
		temp = 0.0;
		for (j = 0; j < n; j++)
			temp = temp + mat_inp[i * n + j] * vec_inp[j];
		vec_out[i] = temp;
	}
}

inline void generate_SND(double* kexi, int n)
{
	// construct a trivial random generator engine from a time-based seed:
	unsigned seed = (unsigned)std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);

	std::normal_distribution<double> distribution(0.0, 1.0);

	for (int i = 0; i < n; i++)
		kexi[i] = distribution(generator);
}

//Cholesky decomposition for matrix A
inline void Cholesky(double* A, double* L, int n)
{
	int k, i, j;
	double sum;

	for (k = 0; k < n; k++)
	{
		//solve gii
		sum = 0;
#pragma omp parallel for reduction(+:sum)
		for (i = 0; i < k; i++)
			sum += L[k * n + i] * L[k * n + i];
		sum = A[k * n + k] - sum;
		L[k * n + k] = sqrt(sum > 1.0e-6 ? sum : 1.0e-6);

		//sovle gij
#pragma omp parallel for private(sum, j)
		for (i = k + 1; i < n; i++)
		{
			sum = 0;
			for (j = 0; j < k; j++)
				sum += L[i * n + j] * L[k * n + j];
			L[i * n + k] = (A[i * n + k] - sum) / L[k * n + k];
		}

#pragma omp parallel for
		for (j = 0; j < k; j++)
			L[j * n + k] = 0;
	}
}

int clFSpatialVariability_Fun::nonstationary_random_2D(Particle* pPar, Para_Fluid* pVFluid, Para_Soil* pParti_ConsPara,
	const Para_Pro& pPPro, const Flag_SpatialVariability fsv, Para_SpatialVariability* psv, const Cell_Con Cellc) {
	double xmin[3], xmax[3], cell_dx, base_x[3], dx[3], acf, lx, ly, n_dx;
	int nx, ny, np_x, np_y, flag[6], dist_type;
	int cell_total, cell_id, temp;
	int i, j, k, pati_num;
	int err_type = 0;
	spv_cell* gSPVC = NULL;
	double* mat_C = NULL;
	double* mat_L = NULL;
	double* vec_zeta = NULL;
	double* vec_E = NULL;
	double dr = pPPro.dr;
	int ntotal = pPPro.ntotal;
	Para_SpatialVariability psv_temp;

	double mean, std_dev, coe_depth;
	double std_temp, temp_rand_var, temp_mean, beta_temp;
	double water_max = Cellc.water_max[1];
	double soil_max = Cellc.soil_max[1];

	//problem domain a
	xmin[0] = Cellc.xmin_nb - 0.1 * (Cellc.xmax_nb - Cellc.xmin_nb);
	xmax[0] = Cellc.xmax_nb + 0.1 * (Cellc.xmax_nb - Cellc.xmin_nb);
	xmin[1] = Cellc.ymin_nb - 0.1 * (Cellc.ymax_nb - Cellc.ymin_nb);
	xmax[1] = Cellc.ymax_nb + 0.1 * (Cellc.ymax_nb - Cellc.ymin_nb);
	xmin[2] = 0.0;
	xmax[2] = 0.0;

	//cell partition
	n_dx = fsv.cell_ndx;
	cell_dx = n_dx * dr;
	nx = (int)((xmax[0] - xmin[0]) / cell_dx) + 1;
	ny = (int)((xmax[1] - xmin[1]) / cell_dx) + 1;
	cell_total = nx * ny;
	temp = cell_total * cell_total;

	if (cell_total > 0 && cell_total < 15000) {
		gSPVC = new spv_cell[cell_total];
		mat_C = new double[temp];
		mat_L = new double[temp];
		vec_zeta = new double[cell_total];
		vec_E = new double[cell_total];
	}
	else {
		err_type = 41;
		return err_type;
	}

	memset(mat_C, 0, sizeof(double) * temp);
	memset(mat_L, 0, sizeof(double) * temp);
	memset(vec_zeta, 0, sizeof(double) * cell_total);
	memset(vec_E, 0, sizeof(double) * cell_total);

	//cell coordinate calculation and particles caclulation
	base_x[0] = -(nx * cell_dx - (xmax[0] - xmin[0])) / 2.0 + 0.01 * dr + xmin[0];
	base_x[1] = -(ny * cell_dx - (xmax[1] - xmin[1])) / 2.0 + 0.01 * dr + xmin[1];
	base_x[2] = 0.0;

	for (i = 0; i < ntotal; i++) {
		if (pPar[i].type == 1 || pPar[i].type == 2) {
			np_x = (int)(fabs(pPar[i].xp[0] - base_x[0]) / cell_dx);
			np_y = (int)(fabs(pPar[i].xp[1] - base_x[1]) / cell_dx);
			cell_id = np_y * nx + np_x;
			if (gSPVC[cell_id].ntotal < 39) {
				gSPVC[cell_id].ntotal = gSPVC[cell_id].ntotal + 1;
				gSPVC[cell_id].parti_id[gSPVC[cell_id].ntotal] = i;
			}
			else {
				err_type = 42;
				return err_type;
			}
		}
	}

	for (j = 0; j < ny; j++) {
		for (i = 0; i < nx; i++) {
			cell_id = j * nx + i;
			gSPVC[cell_id].cent_x[0] = base_x[0] + ((double)i + 0.5) * cell_dx;
			gSPVC[cell_id].cent_x[1] = base_x[1] + ((double)j + 0.5) * cell_dx;
		}
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
			//variables initialization
			psv_temp = psv[i];
			dist_type = psv_temp.type; //1--normal or Guassian distribution; 2--log normal distribution
			mean = psv_temp.mean_value;
			std_dev = psv_temp.sd_value;
			lx = psv_temp.rl_value[0];
			ly = psv_temp.rl_value[1];
			coe_depth = psv_temp.coe_depth;
			beta_temp = fsv.beta * pi / 180.0;

			//autocorrelation function matrix
#pragma omp parallel for schedule(static) private(j,dx,acf)
			for (k = 0; k < cell_total; k++) {
				for (j = 0; j < cell_total; j++) {
					//calculating distance difference
					dx[0] = fabs(gSPVC[k].cent_x[0] - gSPVC[j].cent_x[0]);
					dx[1] = fabs(gSPVC[k].cent_x[1] - gSPVC[j].cent_x[1]);
					dx[2] = 0.0;

					if (fsv.flag_auto == 1) acf = exp(-2.0 * (dx[0] / lx + dx[1] / ly));
					else if (fsv.flag_auto == 2) acf = exp(-pi * (pow(dx[0] / lx, 2) + pow(dx[1] / ly, 2)));
					else acf = exp(-2.0 * ((dx[0] * cos(beta_temp) + dx[1] * sin(beta_temp)) / lx
						+ (-dx[0] * sin(beta_temp) + dx[1] * cos(beta_temp)) / ly));

					mat_C[k * cell_total + j] = acf;
				}
			}

			//calculating Matrix L by cholesky decomposition
			Cholesky(mat_C, mat_L, cell_total);

			//generating random number
			generate_SND(vec_zeta, cell_total);

			//cacluating random field
			matrix_multi(mat_L, vec_zeta, vec_E, cell_total);

			//calculating ramdom variable for each particle
			if (i == 0) { //friction angle for soil
				if (dist_type == 1) { //normal
#pragma omp parallel for schedule(static) private(k,pati_num,temp_rand_var)
					for (j = 0; j < cell_total; j++) {
						pati_num = gSPVC[j].ntotal;
						if (soil_max > gSPVC[j].cent_x[1])
							temp_rand_var = mean + coe_depth * (soil_max - gSPVC[j].cent_x[1]) + std_dev * vec_E[j];
						else
							temp_rand_var = mean + std_dev * vec_E[j];
						for (k = 1; k <= pati_num; k++) {
							if (pPar[gSPVC[j].parti_id[k]].type == 2) pParti_ConsPara[gSPVC[j].parti_id[k]].fai = temp_rand_var;
						}
					}
				}
				else if (dist_type == 2) { //log normal
#pragma omp parallel for schedule(static) private(k,pati_num,temp_rand_var,temp_mean,std_temp)
					for (j = 0; j < cell_total; j++) {
						if (soil_max > gSPVC[j].cent_x[1]) //nonstationary effect
							temp_mean = mean + coe_depth * (soil_max - gSPVC[j].cent_x[1]);
						else 
							temp_mean = mean;
						//lognormal std_dev and mean value
						std_temp = sqrt(log(1.0 + (std_dev / temp_mean) * (std_dev / temp_mean)));
						temp_mean = log(temp_mean) - 0.5 * std_temp * std_temp;
						//random field
						temp_rand_var = temp_mean + std_temp * vec_E[j];
						//put random field information to particle
						pati_num = gSPVC[j].ntotal;
						for (k = 1; k <= pati_num; k++) {
							if (pPar[gSPVC[j].parti_id[k]].type == 2) pParti_ConsPara[gSPVC[j].parti_id[k]].fai = pow(exp(1), temp_rand_var);
						}
					}
				}
			}

			if (i == 1) { //cohesion for soil
				if (dist_type == 1) { //normal
#pragma omp parallel for schedule(static) private(k,pati_num,temp_rand_var)
					for (j = 0; j < cell_total; j++) {
						pati_num = gSPVC[j].ntotal;
						if (soil_max > gSPVC[j].cent_x[1])
							temp_rand_var = mean + coe_depth * (soil_max - gSPVC[j].cent_x[1]) + std_dev * vec_E[j];
						else
							temp_rand_var = mean + std_dev * vec_E[j];
						for (k = 1; k <= pati_num; k++) {
							if (pPar[gSPVC[j].parti_id[k]].type == 2) pParti_ConsPara[gSPVC[j].parti_id[k]].c = temp_rand_var;
						}
					}
				}
				else if (dist_type == 2) { //log normal
#pragma omp parallel for schedule(static) private(k,pati_num,temp_rand_var,temp_mean,std_temp)
					for (j = 0; j < cell_total; j++) {
						if (soil_max > gSPVC[j].cent_x[1]) //nonstationary effect
							temp_mean = mean + coe_depth * (soil_max - gSPVC[j].cent_x[1]);
						else 
							temp_mean = mean;
						//lognormal std_dev and mean value
						std_temp = sqrt(log(1.0 + (std_dev / temp_mean) * (std_dev / temp_mean)));
						temp_mean = log(temp_mean) - 0.5 * std_temp * std_temp;
						//random field
						temp_rand_var = temp_mean + std_temp * vec_E[j];
						//put random field information to particle
						pati_num = gSPVC[j].ntotal;
						for (k = 1; k <= pati_num; k++) {
							if (pPar[gSPVC[j].parti_id[k]].type == 2) pParti_ConsPara[gSPVC[j].parti_id[k]].c = pow(exp(1), temp_rand_var);
						}
					}
				}
			}

			if (i == 2) { //coefficient of permeability for soil
				if (dist_type == 1) { //normal
#pragma omp parallel for schedule(static) private(k,pati_num,temp_rand_var)
					for (j = 0; j < cell_total; j++) {
						pati_num = gSPVC[j].ntotal;
						if (soil_max > gSPVC[j].cent_x[1])
							temp_rand_var = mean + coe_depth * (soil_max - gSPVC[j].cent_x[1]) + std_dev * vec_E[j];
						else
							temp_rand_var = mean + std_dev * vec_E[j];
						for (k = 1; k <= pati_num; k++) {
							if (pPar[gSPVC[j].parti_id[k]].type == 2) pParti_ConsPara[gSPVC[j].parti_id[k]].cop = temp_rand_var;
						}
					}
				}
				else if (dist_type == 2) { //log normal
#pragma omp parallel for schedule(static) private(k,pati_num,temp_rand_var,temp_mean,std_temp)
					for (j = 0; j < cell_total; j++) {
						if (soil_max > gSPVC[j].cent_x[1]) //nonstationary effect
							temp_mean = mean + coe_depth * (soil_max - gSPVC[j].cent_x[1]);
						else 
							temp_mean = mean;
						//lognormal std_dev and mean value
						std_temp = sqrt(log(1.0 + (std_dev / temp_mean) * (std_dev / temp_mean)));
						temp_mean = log(temp_mean) - 0.5 * std_temp * std_temp;
						//random field
						temp_rand_var = temp_mean + std_temp * vec_E[j];
						//put random field information to particle
						pati_num = gSPVC[j].ntotal;
						for (k = 1; k <= pati_num; k++) {
							if (pPar[gSPVC[j].parti_id[k]].type == 2) pParti_ConsPara[gSPVC[j].parti_id[k]].cop = pow(exp(1), temp_rand_var);
						}
					}
				}
			}

			if (i == 3) { //soil grain size for soil
				if (dist_type == 1) { //normal
#pragma omp parallel for schedule(static) private(k,pati_num,temp_rand_var)
					for (j = 0; j < cell_total; j++) {
						pati_num = gSPVC[j].ntotal;
						if (soil_max > gSPVC[j].cent_x[1])
							temp_rand_var = mean + coe_depth * (soil_max - gSPVC[j].cent_x[1]) + std_dev * vec_E[j];
						else
							temp_rand_var = mean + std_dev * vec_E[j];
						for (k = 1; k <= pati_num; k++) {
							if (pPar[gSPVC[j].parti_id[k]].type == 2) pParti_ConsPara[gSPVC[j].parti_id[k]].ds = temp_rand_var;
						}
					}
				}
				else if (dist_type == 2) { //log normal
#pragma omp parallel for schedule(static) private(k,pati_num,temp_rand_var,temp_mean,std_temp)
					for (j = 0; j < cell_total; j++) {
						if (soil_max > gSPVC[j].cent_x[1]) //nonstationary effect
							temp_mean = mean + coe_depth * (soil_max - gSPVC[j].cent_x[1]);
						else 
							temp_mean = mean;
						//lognormal std_dev and mean value
						std_temp = sqrt(log(1.0 + (std_dev / temp_mean) * (std_dev / temp_mean)));
						temp_mean = log(temp_mean) - 0.5 * std_temp * std_temp;
						//random field
						temp_rand_var = temp_mean + std_temp * vec_E[j];
						//put random field information to particle
						pati_num = gSPVC[j].ntotal;
						for (k = 1; k <= pati_num; k++) {
							if (pPar[gSPVC[j].parti_id[k]].type == 2) pParti_ConsPara[gSPVC[j].parti_id[k]].ds = pow(exp(1), temp_rand_var);
						}
					}
				}
			}

			
			if (i == 4) { //friction angle for fluid
				if (dist_type == 1) { //normal
#pragma omp parallel for schedule(static) private(k,pati_num,temp_rand_var)
					for (j = 0; j < cell_total; j++) {
						pati_num = gSPVC[j].ntotal;
						if (water_max > gSPVC[j].cent_x[1])
							temp_rand_var = mean + coe_depth * (water_max - gSPVC[j].cent_x[1]) + std_dev * vec_E[j];
						else
							temp_rand_var = mean + std_dev * vec_E[j];
						for (k = 1; k <= pati_num; k++) {
							if (pPar[gSPVC[j].parti_id[k]].type == 1) pVFluid[gSPVC[j].parti_id[k]].fai = temp_rand_var;
						}
					}
				}
				else if (dist_type == 2) { //log normal
#pragma omp parallel for schedule(static) private(k,pati_num,temp_rand_var,temp_mean,std_temp)
					for (j = 0; j < cell_total; j++) {
						if (soil_max > gSPVC[j].cent_x[1]) //nonstationary effect
							temp_mean = mean + coe_depth * (soil_max - gSPVC[j].cent_x[1]);
						else 
							temp_mean = mean;
						//lognormal std_dev and mean value
						std_temp = sqrt(log(1.0 + (std_dev / temp_mean) * (std_dev / temp_mean)));
						temp_mean = log(temp_mean) - 0.5 * std_temp * std_temp;
						//random field
						temp_rand_var = temp_mean + std_temp * vec_E[j];
						//put random field information to particle
						pati_num = gSPVC[j].ntotal;
						for (k = 1; k <= pati_num; k++) {
							if (pPar[gSPVC[j].parti_id[k]].type == 1) pVFluid[gSPVC[j].parti_id[k]].fai = pow(exp(1), temp_rand_var);
						}
					}
				}
			}

			if (i == 5) { //cohesion for fluid
				if (dist_type == 1) { //normal
#pragma omp parallel for schedule(static) private(k,pati_num,temp_rand_var)
					for (j = 0; j < cell_total; j++) {
						pati_num = gSPVC[j].ntotal;
						if (water_max > gSPVC[j].cent_x[1])
							temp_rand_var = mean + coe_depth * (water_max - gSPVC[j].cent_x[1]) + std_dev * vec_E[j];
						else
							temp_rand_var = mean + std_dev * vec_E[j];
						for (k = 1; k <= pati_num; k++) {
							if (pPar[gSPVC[j].parti_id[k]].type == 1) pVFluid[gSPVC[j].parti_id[k]].c = temp_rand_var;
						}
					}
				}
				else if (dist_type == 2) { //log normal
#pragma omp parallel for schedule(static) private(k,pati_num,temp_rand_var,temp_mean,std_temp)
					for (j = 0; j < cell_total; j++) {
						if (soil_max > gSPVC[j].cent_x[1]) //nonstationary effect
							temp_mean = mean + coe_depth * (soil_max - gSPVC[j].cent_x[1]);
						else 
							temp_mean = mean;
						//lognormal std_dev and mean value
						std_temp = sqrt(log(1.0 + (std_dev / temp_mean) * (std_dev / temp_mean)));
						temp_mean = log(temp_mean) - 0.5 * std_temp * std_temp;
						//random field
						temp_rand_var = temp_mean + std_temp * vec_E[j];
						//put random field information to particle
						pati_num = gSPVC[j].ntotal;
						for (k = 1; k <= pati_num; k++) {
							if (pPar[gSPVC[j].parti_id[k]].type == 1) pVFluid[gSPVC[j].parti_id[k]].c = pow(exp(1), temp_rand_var);
						}
					}
				}
			}
		}
	}

	//free memory
	if (gSPVC != NULL) delete[]gSPVC;
	if (mat_C != NULL) delete[]mat_C;
	if (mat_L != NULL) delete[]mat_L;
	if (vec_zeta != NULL) delete[]vec_zeta;
	if (vec_E != NULL) delete[]vec_E;

	return 0;
}

int clFSpatialVariability_Fun::nonstationary_random_3D(Particle* pPar, Para_Fluid* pVFluid, Para_Soil* pParti_ConsPara,
	const Para_Pro& pPPro, const Flag_SpatialVariability fsv, Para_SpatialVariability* psv, const Cell_Con Cellc) {
	double xmin[3], xmax[3], cell_dx, base_x[3], dx[3], acf, lx, ly, lz, n_dx, temp_dx[3];
	int nx, ny, nz, np_x, np_y, np_z, flag[6], dist_type;
	int cell_total, cell_id, temp;
	int i, j, k, pati_num;
	int err_type = 0;
	spv_cell* gSPVC = NULL;
	double* mat_C = NULL;
	double* mat_L = NULL;
	double* vec_zeta = NULL;
	double* vec_E = NULL;
	double dr = pPPro.dr;
	int ntotal = pPPro.ntotal;
	Para_SpatialVariability psv_temp;

	double mean, std_dev, coe_depth;
	double std_temp, temp_rand_var, temp_mean, beta_temp, beta1_temp;
	double water_max = Cellc.water_max[2];
	double soil_max = Cellc.soil_max[2];

	//problem domain a
	xmin[0] = Cellc.xmin_nb - 0.1 * (Cellc.xmax_nb - Cellc.xmin_nb);
	xmax[0] = Cellc.xmax_nb + 0.1 * (Cellc.xmax_nb - Cellc.xmin_nb);
	xmin[1] = Cellc.ymin_nb - 0.1 * (Cellc.ymax_nb - Cellc.ymin_nb);
	xmax[1] = Cellc.ymax_nb + 0.1 * (Cellc.ymax_nb - Cellc.ymin_nb);
	xmin[2] = Cellc.zmin_nb - 0.1 * (Cellc.zmax_nb - Cellc.zmin_nb);
	xmax[2] = Cellc.zmax_nb + 0.1 * (Cellc.zmax_nb - Cellc.zmin_nb);

	//cell partition
	n_dx = fsv.cell_ndx;
	cell_dx = n_dx * dr;
	nx = (int)((xmax[0] - xmin[0]) / cell_dx) + 1;
	ny = (int)((xmax[1] - xmin[1]) / cell_dx) + 1;
	nz = (int)((xmax[2] - xmin[2]) / cell_dx) + 1;
	cell_total = nx * ny * nz;
	temp = cell_total * cell_total;

	if (cell_total > 0 && cell_total < 15000) {
		gSPVC = new spv_cell[cell_total];
		mat_C = new double[temp];
		mat_L = new double[temp];
		vec_zeta = new double[cell_total];
		vec_E = new double[cell_total];
	}
	else {
		err_type = 41;
		return err_type;
	}

	memset(mat_C, 0, sizeof(double) * temp);
	memset(mat_L, 0, sizeof(double) * temp);
	memset(vec_zeta, 0, sizeof(double) * cell_total);
	memset(vec_E, 0, sizeof(double) * cell_total);

	//cell coordinate calculation and particles caclulation
	base_x[0] = -(nx * cell_dx - (xmax[0] - xmin[0])) / 2.0 + 0.01 * dr + xmin[0];
	base_x[1] = -(ny * cell_dx - (xmax[1] - xmin[1])) / 2.0 + 0.01 * dr + xmin[1];
	base_x[2] = -(nz * cell_dx - (xmax[2] - xmin[2])) / 2.0 + 0.01 * dr + xmin[2];

	for (i = 0; i < ntotal; i++) {
		np_x = (int)(fabs(pPar[i].xp[0] - base_x[0]) / cell_dx);
		np_y = (int)(fabs(pPar[i].xp[1] - base_x[1]) / cell_dx);
		np_z = (int)(fabs(pPar[i].xp[2] - base_x[2]) / cell_dx);

		cell_id = np_z * nx * ny + np_y * nx + np_x;
		if (gSPVC[cell_id].ntotal < 39) {
			gSPVC[cell_id].ntotal = gSPVC[cell_id].ntotal + 1;
			gSPVC[cell_id].parti_id[gSPVC[cell_id].ntotal] = i;
		}
		else {
			err_type = 42;
			return err_type;
		}
	}
	for (k = 0; k < nz; k++) {
		for (j = 0; j < ny; j++) {
			for (i = 0; i < nx; i++) {
				cell_id = k * nx * ny + j * nx + i;
				gSPVC[cell_id].cent_x[0] = base_x[0] + ((double)i + 0.5) * cell_dx;
				gSPVC[cell_id].cent_x[1] = base_x[1] + ((double)j + 0.5) * cell_dx;
				gSPVC[cell_id].cent_x[2] = base_x[2] + ((double)k + 0.5) * cell_dx;
			}
		}
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
			//variables initialization
			psv_temp = psv[i];
			dist_type = psv_temp.type; //1--normal or Guassian distribution; 2--log normal distribution
			mean = psv_temp.mean_value;
			std_dev = psv_temp.sd_value;
			lx = psv_temp.rl_value[0];
			ly = psv_temp.rl_value[1];
			lz = psv_temp.rl_value[2];
			coe_depth = psv_temp.coe_depth;
			beta_temp = fsv.beta * pi / 180.0;
			beta1_temp = fsv.beta_1 * pi / 180.0;

			//autocorrelation function matrix
#pragma omp parallel for schedule(static) private(j,dx,acf,temp_dx)
			for (k = 0; k < cell_total; k++) {
				for (j = 0; j < cell_total; j++) {
					//calculating distance difference
					dx[0] = fabs(gSPVC[k].cent_x[0] - gSPVC[j].cent_x[0]);
					dx[1] = fabs(gSPVC[k].cent_x[1] - gSPVC[j].cent_x[1]);
					dx[2] = fabs(gSPVC[k].cent_x[2] - gSPVC[j].cent_x[2]);

					temp_dx[0] = dx[0] * cos(beta1_temp) + dx[1] * sin(beta1_temp);
					temp_dx[1] = (-dx[0] * sin(beta1_temp) + dx[1] * cos(beta1_temp)) * cos(beta_temp) + dx[2] * sin(beta_temp);
					temp_dx[2] = -(-dx[0] * sin(beta1_temp) + dx[1] * cos(beta1_temp)) * sin(beta_temp) + dx[2] * cos(beta_temp);

					if (fsv.flag_auto == 1) acf = exp(-2.0 * (dx[0] / lx + dx[1] / ly + dx[2] / lz));
					else if (fsv.flag_auto == 2) acf = exp(-pi * (pow(dx[0] / lx, 2) + pow(dx[1] / ly, 2) + pow(dx[2] / lz, 2)));
					else acf = exp(-2.0 * (temp_dx[0] / lx + temp_dx[1] / ly + temp_dx[2] / lz));

					mat_C[k * cell_total + j] = acf;
				}
			}

			//calculating Matrix L by cholesky decomposition
			Cholesky(mat_C, mat_L, cell_total);

			//generating random number
			generate_SND(vec_zeta, cell_total);

			//cacluating random field
			matrix_multi(mat_L, vec_zeta, vec_E, cell_total);

			//calculating ramdom variable for each particle
			if (i == 0) { //friction angle for soil
				if (dist_type == 1) { //normal
#pragma omp parallel for schedule(static) private(k,pati_num,temp_rand_var)
					for (j = 0; j < cell_total; j++) {
						pati_num = gSPVC[j].ntotal;
						if (soil_max > gSPVC[j].cent_x[2])
							temp_rand_var = mean + coe_depth * (soil_max - gSPVC[j].cent_x[2]) + std_dev * vec_E[j];
						else
							temp_rand_var = mean + std_dev * vec_E[j];
						for (k = 1; k <= pati_num; k++) {
							if (pPar[gSPVC[j].parti_id[k]].type == 2) pParti_ConsPara[gSPVC[j].parti_id[k]].fai = temp_rand_var;
						}
					}
				}
				else if (dist_type == 2) { //log normal
#pragma omp parallel for schedule(static) private(k,pati_num,temp_rand_var,temp_mean,std_temp)
					for (j = 0; j < cell_total; j++) {
						if (soil_max > gSPVC[j].cent_x[2]) //nonstationary effect
							temp_mean = mean + coe_depth * (soil_max - gSPVC[j].cent_x[2]);
						else 
							temp_mean = mean;
						//lognormal std_dev and mean value
						std_temp = sqrt(log(1.0 + (std_dev / temp_mean) * (std_dev / temp_mean)));
						temp_mean = log(temp_mean) - 0.5 * std_temp * std_temp;
						//random field
						temp_rand_var = temp_mean + std_temp * vec_E[j];
						//put random field information to particle
						pati_num = gSPVC[j].ntotal;
						for (k = 1; k <= pati_num; k++) {
							if (pPar[gSPVC[j].parti_id[k]].type == 2) pParti_ConsPara[gSPVC[j].parti_id[k]].fai = pow(exp(1), temp_rand_var);
						}
					}
				}
			}

			if (i == 1) { //cohesion for soil
				if (dist_type == 1) { //normal
#pragma omp parallel for schedule(static) private(k,pati_num,temp_rand_var)
					for (j = 0; j < cell_total; j++) {
						pati_num = gSPVC[j].ntotal;
						if (soil_max > gSPVC[j].cent_x[2])
							temp_rand_var = mean + coe_depth * (soil_max - gSPVC[j].cent_x[2]) + std_dev * vec_E[j];
						else
							temp_rand_var = mean + std_dev * vec_E[j];
						for (k = 1; k <= pati_num; k++) {
							if (pPar[gSPVC[j].parti_id[k]].type == 2) pParti_ConsPara[gSPVC[j].parti_id[k]].c = temp_rand_var;
						}
					}
				}
				else if (dist_type == 2) { //log normal
#pragma omp parallel for schedule(static) private(k,pati_num,temp_rand_var,temp_mean,std_temp)
					for (j = 0; j < cell_total; j++) {
						if (soil_max > gSPVC[j].cent_x[2]) //nonstationary effect
							temp_mean = mean + coe_depth * (soil_max - gSPVC[j].cent_x[2]);
						else 
							temp_mean = mean;
						//lognormal std_dev and mean value
						std_temp = sqrt(log(1.0 + (std_dev / temp_mean) * (std_dev / temp_mean)));
						temp_mean = log(temp_mean) - 0.5 * std_temp * std_temp;
						//random field
						temp_rand_var = temp_mean + std_temp * vec_E[j];
						//put random field information to particle
						pati_num = gSPVC[j].ntotal;
						for (k = 1; k <= pati_num; k++) {
							if (pPar[gSPVC[j].parti_id[k]].type == 2) pParti_ConsPara[gSPVC[j].parti_id[k]].c = pow(exp(1), temp_rand_var);
						}
					}
				}
			}

			if (i == 2) { //coefficient of permeability for soil
				if (dist_type == 1) { //normal
#pragma omp parallel for schedule(static) private(k,pati_num,temp_rand_var)
					for (j = 0; j < cell_total; j++) {
						pati_num = gSPVC[j].ntotal;
						if (soil_max > gSPVC[j].cent_x[2])
							temp_rand_var = mean + coe_depth * (soil_max - gSPVC[j].cent_x[2]) + std_dev * vec_E[j];
						else
							temp_rand_var = mean + std_dev * vec_E[j];
						for (k = 1; k <= pati_num; k++) {
							if (pPar[gSPVC[j].parti_id[k]].type == 2) pParti_ConsPara[gSPVC[j].parti_id[k]].cop = temp_rand_var;
						}
					}
				}
				else if (dist_type == 2) { //log normal
#pragma omp parallel for schedule(static) private(k,pati_num,temp_rand_var,temp_mean,std_temp)
					for (j = 0; j < cell_total; j++) {
						if (soil_max > gSPVC[j].cent_x[2]) //nonstationary effect
							temp_mean = mean + coe_depth * (soil_max - gSPVC[j].cent_x[2]);
						else 
							temp_mean = mean;
						//lognormal std_dev and mean value
						std_temp = sqrt(log(1.0 + (std_dev / temp_mean) * (std_dev / temp_mean)));
						temp_mean = log(temp_mean) - 0.5 * std_temp * std_temp;
						//random field
						temp_rand_var = temp_mean + std_temp * vec_E[j];
						//put random field information to particle
						pati_num = gSPVC[j].ntotal;
						for (k = 1; k <= pati_num; k++) {
							if (pPar[gSPVC[j].parti_id[k]].type == 2) pParti_ConsPara[gSPVC[j].parti_id[k]].cop = pow(exp(1), temp_rand_var);
						}
					}
				}
			}

			if (i == 3) { //soil grain size for soil
				if (dist_type == 1) { //normal
#pragma omp parallel for schedule(static) private(k,pati_num,temp_rand_var)
					for (j = 0; j < cell_total; j++) {
						pati_num = gSPVC[j].ntotal;
						if (soil_max > gSPVC[j].cent_x[2])
							temp_rand_var = mean + coe_depth * (soil_max - gSPVC[j].cent_x[2]) + std_dev * vec_E[j];
						else
							temp_rand_var = mean + std_dev * vec_E[j];
						for (k = 1; k <= pati_num; k++) {
							if (pPar[gSPVC[j].parti_id[k]].type == 2) pParti_ConsPara[gSPVC[j].parti_id[k]].ds = temp_rand_var;
						}
					}
				}
				else if (dist_type == 2) { //log normal
#pragma omp parallel for schedule(static) private(k,pati_num,temp_rand_var,temp_mean,std_temp)
					for (j = 0; j < cell_total; j++) {
						if (soil_max > gSPVC[j].cent_x[2]) //nonstationary effect
							temp_mean = mean + coe_depth * (soil_max - gSPVC[j].cent_x[2]);
						else 
							temp_mean = mean;
						//lognormal std_dev and mean value
						std_temp = sqrt(log(1.0 + (std_dev / temp_mean) * (std_dev / temp_mean)));
						temp_mean = log(temp_mean) - 0.5 * std_temp * std_temp;
						//random field
						temp_rand_var = temp_mean + std_temp * vec_E[j];
						//put random field information to particle
						pati_num = gSPVC[j].ntotal;
						for (k = 1; k <= pati_num; k++) {
							if (pPar[gSPVC[j].parti_id[k]].type == 2) pParti_ConsPara[gSPVC[j].parti_id[k]].ds = pow(exp(1), temp_rand_var);
						}
					}
				}
			}

			
			if (i == 4) { //friction angle for fluid
				if (dist_type == 1) { //normal
#pragma omp parallel for schedule(static) private(k,pati_num,temp_rand_var)
					for (j = 0; j < cell_total; j++) {
						pati_num = gSPVC[j].ntotal;
						if (water_max > gSPVC[j].cent_x[2])
							temp_rand_var = mean + coe_depth * (water_max - gSPVC[j].cent_x[2]) + std_dev * vec_E[j];
						else
							temp_rand_var = mean + std_dev * vec_E[j];
						for (k = 1; k <= pati_num; k++) {
							if (pPar[gSPVC[j].parti_id[k]].type == 1) pVFluid[gSPVC[j].parti_id[k]].fai = temp_rand_var;
						}
					}
				}
				else if (dist_type == 2) { //log normal
#pragma omp parallel for schedule(static) private(k,pati_num,temp_rand_var,temp_mean,std_temp)
					for (j = 0; j < cell_total; j++) {
						if (soil_max > gSPVC[j].cent_x[2]) //nonstationary effect
							temp_mean = mean + coe_depth * (soil_max - gSPVC[j].cent_x[2]);
						else 
							temp_mean = mean;
						//lognormal std_dev and mean value
						std_temp = sqrt(log(1.0 + (std_dev / temp_mean) * (std_dev / temp_mean)));
						temp_mean = log(temp_mean) - 0.5 * std_temp * std_temp;
						//random field
						temp_rand_var = temp_mean + std_temp * vec_E[j];
						//put random field information to particle
						pati_num = gSPVC[j].ntotal;
						for (k = 1; k <= pati_num; k++) {
							if (pPar[gSPVC[j].parti_id[k]].type == 1) pVFluid[gSPVC[j].parti_id[k]].fai = pow(exp(1), temp_rand_var);
						}
					}
				}
			}

			if (i == 5) { //cohesion for fluid
				if (dist_type == 1) { //normal
#pragma omp parallel for schedule(static) private(k,pati_num,temp_rand_var)
					for (j = 0; j < cell_total; j++) {
						pati_num = gSPVC[j].ntotal;
						if (water_max > gSPVC[j].cent_x[2])
							temp_rand_var = mean + coe_depth * (water_max - gSPVC[j].cent_x[2]) + std_dev * vec_E[j];
						else
							temp_rand_var = mean + std_dev * vec_E[j];
						for (k = 1; k <= pati_num; k++) {
							if (pPar[gSPVC[j].parti_id[k]].type == 1) pVFluid[gSPVC[j].parti_id[k]].c = temp_rand_var;
						}
					}
				}
				else if (dist_type == 2) { //log normal
#pragma omp parallel for schedule(static) private(k,pati_num,temp_rand_var,temp_mean,std_temp)
					for (j = 0; j < cell_total; j++) {
						if (soil_max > gSPVC[j].cent_x[2]) //nonstationary effect
							temp_mean = mean + coe_depth * (soil_max - gSPVC[j].cent_x[2]);
						else 
							temp_mean = mean;
						//lognormal std_dev and mean value
						std_temp = sqrt(log(1.0 + (std_dev / temp_mean) * (std_dev / temp_mean)));
						temp_mean = log(temp_mean) - 0.5 * std_temp * std_temp;
						//random field
						temp_rand_var = temp_mean + std_temp * vec_E[j];
						//put random field information to particle
						pati_num = gSPVC[j].ntotal;
						for (k = 1; k <= pati_num; k++) {
							if (pPar[gSPVC[j].parti_id[k]].type == 1) pVFluid[gSPVC[j].parti_id[k]].c = pow(exp(1), temp_rand_var);
						}
					}
				}
			}
		}
	}

	//free memory
	if (gSPVC != NULL) delete[]gSPVC;
	if (mat_C != NULL) delete[]mat_C;
	if (mat_L != NULL) delete[]mat_L;
	if (vec_zeta != NULL) delete[]vec_zeta;
	if (vec_E != NULL) delete[]vec_E;

	return 0;
}