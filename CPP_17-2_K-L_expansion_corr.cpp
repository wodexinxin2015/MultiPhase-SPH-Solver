/************************************************************************************
Multi-Phase and Parallelized SPH program
--parallelized and optimized for SMP (Symmetric multiprocessing ) architecture.
--Copyright (c) 2018-2023, Weijie ZHANG, GeoHohai, Hohai University.
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

//inline code for matrix multiply: [m*n]*[n*p]=[m*p]
inline void matrix_multi(double* mat_inp, double* vec_inp, double* vec_out,
	int m, int n, int p) {
	int i, j, k;
	double temp;
#pragma omp parallel for schedule(static) private(j,k,temp)
	for (i = 0; i < m; i++) {
		for (j = 0; j < p; j++) {
			temp = 0.0;
			for (k = 0; k < n; k++) {
				temp = temp + mat_inp[i * n + k] * vec_inp[k * p + j];
			}
			vec_out[i * p + j] = temp;
		}
	}
}

//generating a random variable that corresponds with Standard Normal Distribution
inline void generate_SND(double *kexi, int n)
{
	// construct a trivial random generator engine from a time-based seed:
	unsigned seed = (unsigned)std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);

	std::normal_distribution<double> distribution(0.0, 1.0);

	for (int i = 0; i < n; i++)
	{
		kexi[i * 2] = distribution(generator);
		kexi[i * 2 + 1] = distribution(generator);
	}
}

//generate the spatial-distributed variables
int clFSpatialVariability_Fun::SpatialVariables_Generate_2D_corr(Particle* pPar, Para_Fluid* pVFluid, Para_Soil* pParti_ConsPara,
	const Para_Pro& pPPro, const Flag_SpatialVariability fsv, Para_SpatialVariability* psv, const Cell_Con Cellc)
{
	double xmin[3], xmax[3], x_trans[3];
	int i, j, order_id[200][3], flag[6];
	int err_type = 0;
	double uu[klterm], wi[klterm][3], eigval[klterm];
	double* uu_cfai = new double[klterm * 2];
	double* zeta = new double[klterm * 2];
	double water_max, soil_max;
	double covariance_mat[4];
	water_max = Cellc.water_max[1];
	soil_max = Cellc.soil_max[1];

	//problem domain a
	xmin[0] = Cellc.xmin_nb - 0.2 * (Cellc.xmax_nb - Cellc.xmin_nb);
	xmax[0] = Cellc.xmax_nb + 0.2 * (Cellc.xmax_nb - Cellc.xmin_nb);
	xmin[1] = Cellc.ymin_nb - 0.2 * (Cellc.ymax_nb - Cellc.ymin_nb);
	xmax[1] = Cellc.ymax_nb + 0.2 * (Cellc.ymax_nb - Cellc.ymin_nb);

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
	if (fsv.flag_fai == 1) //correlated flag
	{
		flag[0] = 1;
	}
	else flag[0] = 0;

	if (fsv.flag_fai == 1) //correlated flag
	{
		flag[1] = 1;
	}
	else flag[1] = 0;

	if (fsv.flag_faiw == 1) //correlated flag
	{
		flag[4] = 1;
	}
	else flag[4] = 0;

	if (fsv.flag_faiw == 1) //correlated flag
	{
		flag[5] = 1;
	}
	else flag[5] = 0;
	
	if (flag[0] == 1) {
		int dist_type = psv[0].type;
		double coefficient = fsv.coe_frict_cohesion;

		// generate the uu for the fai and c of soil
		generate_SND(uu_cfai, klterm);

		//calculating the correlate efficient
		if (dist_type == 2) {  //log normal distribution for soil
			coefficient = log(1.0 +
				fsv.coe_frict_cohesion * (psv[0].sd_value / psv[0].mean_value) * (psv[1].sd_value / psv[1].mean_value));
			coefficient = coefficient / sqrt(log(1.0 + psv[0].sd_value * psv[0].sd_value / psv[0].mean_value / psv[0].mean_value)
				* log(1.0 + psv[1].sd_value * psv[1].sd_value / psv[1].mean_value / psv[1].mean_value));
		}
		else
			coefficient = fsv.coe_frict_cohesion;

		//correlate covariance matrix between friction angle and cohesion
		covariance_mat[0] = 1.0;
		covariance_mat[1] = coefficient;
		covariance_mat[2] = 0.0;
		covariance_mat[3] = sqrt(1.0 - coefficient * coefficient);

		//calculate covariance matrix correlated random number
		matrix_multi(uu_cfai, covariance_mat, zeta, klterm, 2, 2);

		//generate the random field for fai
		for (i = 0; i < klterm; i++) {
			uu[i] = zeta[2 * i];
		}

		//solve wi
		solve_wi(wi, psv[0].a, psv[0].rl_value, 2);

		if (err_type != 0)
			return err_type;

		//solve ramda
		err_type = solve_ramdai(eigval, wi, psv[0].rl_value, order_id, 2);
		if (err_type != 0)
			return err_type;

		//solve H
		err_type = solve_h_2D(pPar, pVFluid, pParti_ConsPara, pPPro, wi, eigval, uu, x_trans, order_id, psv[0], soil_max, water_max, 0);
		if (err_type != 0)
			return err_type;

		//generate the random field for c
		for (i = 0; i < klterm; i++) {
			uu[i] = zeta[2 * i + 1];
		}
		
		//solve wi
		solve_wi(wi, psv[1].a, psv[1].rl_value, 2);

		if (err_type != 0)
			return err_type;

		//solve ramda
		err_type = solve_ramdai(eigval, wi, psv[1].rl_value, order_id, 2);
		if (err_type != 0)
			return err_type;

		//solve H
		err_type = solve_h_2D(pPar, pVFluid, pParti_ConsPara, pPPro, wi, eigval, uu, x_trans, order_id, psv[1], soil_max, water_max, 1);
		if (err_type != 0)
			return err_type;
	}

	if (flag[4] == 1) {
		int dist_type = psv[4].type;
		double coefficient = fsv.coe_frict_cohesion;

		// generate the uu for the fai and c of fluid
		generate_SND(uu_cfai, klterm);

		//calculating the correlate efficient
		if (dist_type == 2) {  //log normal distribution for fluid
			coefficient = log(1.0 +
				fsv.coe_frict_cohesion * (psv[4].sd_value / psv[4].mean_value) * (psv[5].sd_value / psv[5].mean_value));
			coefficient = coefficient / sqrt(log(1.0 + psv[4].sd_value * psv[4].sd_value / psv[4].mean_value / psv[4].mean_value)
				* log(1.0 + psv[5].sd_value * psv[5].sd_value / psv[5].mean_value / psv[5].mean_value));
		}
		else
			coefficient = fsv.coe_frict_cohesion;

		//correlate covariance matrix between friction angle and cohesion
		covariance_mat[0] = 1.0;
		covariance_mat[1] = coefficient;
		covariance_mat[2] = 0.0;
		covariance_mat[3] = sqrt(1.0 - coefficient * coefficient);

		//calculate covariance matrix correlated random number
		matrix_multi(uu_cfai, covariance_mat, zeta, klterm, 2, 2);

		//generate the random field for fai
		for (i = 0; i < klterm; i++) {
			uu[i] = zeta[2 * i];
		}

		//solve wi
		solve_wi(wi, psv[4].a, psv[4].rl_value, 2);

		if (err_type != 0)
			return err_type;

		//solve ramda
		err_type = solve_ramdai(eigval, wi, psv[4].rl_value, order_id, 2);
		if (err_type != 0)
			return err_type;

		//solve H
		err_type = solve_h_2D(pPar, pVFluid, pParti_ConsPara, pPPro, wi, eigval, uu, x_trans, order_id, psv[4], soil_max, water_max, 4);
		if (err_type != 0)
			return err_type;

		//generate the random field for c
		for (i = 0; i < klterm; i++) {
			uu[i] = zeta[2 * i + 1];
		}

		//solve wi
		solve_wi(wi, psv[5].a, psv[5].rl_value, 2);

		if (err_type != 0)
			return err_type;

		//solve ramda
		err_type = solve_ramdai(eigval, wi, psv[5].rl_value, order_id, 2);
		if (err_type != 0)
			return err_type;

		//solve H
		err_type = solve_h_2D(pPar, pVFluid, pParti_ConsPara, pPPro, wi, eigval, uu, x_trans, order_id, psv[5], soil_max, water_max, 5);
		if (err_type != 0)
			return err_type;
	}

	//output random variables
	output_Para_SpatialVariability(fsv, psv);

	//free memory
	delete []uu_cfai;
	delete []zeta;

	return 0;
}

//generate the spatial-distributed variables
int clFSpatialVariability_Fun::SpatialVariables_Generate_3D_corr(Particle* pPar, Para_Fluid* pVFluid, Para_Soil* pParti_ConsPara,
	const Para_Pro& pPPro, const Flag_SpatialVariability fsv, Para_SpatialVariability* psv, const Cell_Con Cellc)
{
	double xmin[3], xmax[3], x_trans[3];
	int i, j, order_id[200][3], flag[6];
	int err_type = 0;
	double uu[klterm], wi[klterm][3], eigval[klterm];
	double* uu_cfai = new double[klterm * 2];
	double* zeta = new double[klterm * 2];
	double water_max, soil_max;
	double covariance_mat[4];
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
	if (fsv.flag_fai == 1) //correlated flag
	{
		flag[0] = 1;
	}
	else flag[0] = 0;

	if (fsv.flag_fai == 1) //correlated flag
	{
		flag[1] = 1;
	}
	else flag[1] = 0;

	if (fsv.flag_faiw == 1) //correlated flag
	{
		flag[4] = 1;
	}
	else flag[4] = 0;

	if (fsv.flag_faiw == 1) //correlated flag
	{
		flag[5] = 1;
	}
	else flag[5] = 0;

	if (flag[0] == 1) {
		int dist_type = psv[0].type;
		double coefficient = fsv.coe_frict_cohesion;

		// generate the uu for the fai and c of soil
		generate_SND(uu_cfai, klterm);

		//calculating the correlate efficient
		if (dist_type == 2) {  //log normal distribution for soil
			coefficient = log(1.0 +
				fsv.coe_frict_cohesion * (psv[0].sd_value / psv[0].mean_value) * (psv[1].sd_value / psv[1].mean_value));
			coefficient = coefficient / sqrt(log(1.0 + psv[0].sd_value * psv[0].sd_value / psv[0].mean_value / psv[0].mean_value)
				* log(1.0 + psv[1].sd_value * psv[1].sd_value / psv[1].mean_value / psv[1].mean_value));
		}
		else
			coefficient = fsv.coe_frict_cohesion;

		//correlate covariance matrix between friction angle and cohesion
		covariance_mat[0] = 1.0;
		covariance_mat[1] = coefficient;
		covariance_mat[2] = 0.0;
		covariance_mat[3] = sqrt(1.0 - coefficient * coefficient);

		//calculate covariance matrix correlated random number
		matrix_multi(uu_cfai, covariance_mat, zeta, klterm, 2, 2);

		//generate the random field for fai
		for (i = 0; i < klterm; i++) {
			uu[i] = zeta[2 * i];
		}

		//solve wi
		solve_wi(wi, psv[0].a, psv[0].rl_value, 3);

		if (err_type != 0)
			return err_type;

		//solve ramda
		err_type = solve_ramdai(eigval, wi, psv[0].rl_value, order_id, 3);
		if (err_type != 0)
			return err_type;

		//solve H
		err_type = solve_h_3D(pPar, pVFluid, pParti_ConsPara, pPPro, wi, eigval, uu, x_trans, order_id, psv[0], soil_max, water_max, 0);
		if (err_type != 0)
			return err_type;

		//generate the random field for c
		for (i = 0; i < klterm; i++) {
			uu[i] = zeta[2 * i + 1];
		}

		//solve wi
		solve_wi(wi, psv[1].a, psv[1].rl_value, 3);

		if (err_type != 0)
			return err_type;

		//solve ramda
		err_type = solve_ramdai(eigval, wi, psv[1].rl_value, order_id, 3);
		if (err_type != 0)
			return err_type;

		//solve H
		err_type = solve_h_3D(pPar, pVFluid, pParti_ConsPara, pPPro, wi, eigval, uu, x_trans, order_id, psv[1], soil_max, water_max, 1);
		if (err_type != 0)
			return err_type;
	}

	if (flag[4] == 1) {
		int dist_type = psv[4].type;
		double coefficient = fsv.coe_frict_cohesion;

		// generate the uu for the fai and c of fluid
		generate_SND(uu_cfai, klterm);

		//calculating the correlate efficient
		if (dist_type == 2) {  //log normal distribution for fluid
			coefficient = log(1.0 +
				fsv.coe_frict_cohesion * (psv[4].sd_value / psv[4].mean_value) * (psv[5].sd_value / psv[5].mean_value));
			coefficient = coefficient / sqrt(log(1.0 + psv[4].sd_value * psv[4].sd_value / psv[4].mean_value / psv[4].mean_value)
				* log(1.0 + psv[5].sd_value * psv[5].sd_value / psv[5].mean_value / psv[5].mean_value));
		}
		else
			coefficient = fsv.coe_frict_cohesion;

		//correlate covariance matrix between friction angle and cohesion
		covariance_mat[0] = 1.0;
		covariance_mat[1] = coefficient;
		covariance_mat[2] = 0.0;
		covariance_mat[3] = sqrt(1.0 - coefficient * coefficient);

		//calculate covariance matrix correlated random number
		matrix_multi(uu_cfai, covariance_mat, zeta, klterm, 2, 2);

		//generate the random field for fai
		for (i = 0; i < klterm; i++) {
			uu[i] = zeta[2 * i];
		}

		//solve wi
		solve_wi(wi, psv[4].a, psv[4].rl_value, 3);

		if (err_type != 0)
			return err_type;

		//solve ramda
		err_type = solve_ramdai(eigval, wi, psv[4].rl_value, order_id, 3);
		if (err_type != 0)
			return err_type;

		//solve H
		err_type = solve_h_3D(pPar, pVFluid, pParti_ConsPara, pPPro, wi, eigval, uu, x_trans, order_id, psv[4], soil_max, water_max, 4);
		if (err_type != 0)
			return err_type;

		//generate the random field for c
		for (i = 0; i < klterm; i++) {
			uu[i] = zeta[2 * i + 1];
		}

		//solve wi
		solve_wi(wi, psv[5].a, psv[5].rl_value, 3);

		if (err_type != 0)
			return err_type;

		//solve ramda
		err_type = solve_ramdai(eigval, wi, psv[5].rl_value, order_id, 3);
		if (err_type != 0)
			return err_type;

		//solve H
		err_type = solve_h_3D(pPar, pVFluid, pParti_ConsPara, pPPro, wi, eigval, uu, x_trans, order_id, psv[5], soil_max, water_max, 5);
		if (err_type != 0)
			return err_type;
	}

	//output random variables
	output_Para_SpatialVariability(fsv, psv);

	//free memory
	delete[]uu_cfai;
	delete[]zeta;

	return 0;
}