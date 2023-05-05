/************************************************************************************
Multi-Phase and Parallelized SPH program
--parallelized and optimized for SMP (Symmetric multiprocessing ) architecture.
--Copyright (c) 2018-2023, Weijie ZHANG and Jian JI, GeoHohai, Hohai University.
************************************************************************************/

#pragma once
#include "Class_particle.h"
#include "Class_ParticleProperty.h"
#include "Class_NNPS.h"

class Flag_SpatialVariability
{
public:
	//common flag for the spatial variability of solid phase
	int flag_fai;	//flag for the internal frictional angle
	int flag_c;		//flag for the cohesion
	int flag_cop;   //flag for the coefficient of permeability
	int flag_ds;	//flag for the soil particle size, not the particle spacing

	//flags  for the spatial variability of fluid phase
	int flag_faiw; //flag for the internal frictional angle
	int flag_cw;   //flag for the cohesion

	double cell_ndx; //set the cell dimension: length = cell_ndx * dr

	//correlation coefficient
	double coe_frict_cohesion; //correlation coefficient between firction angle and cohesion

	//flag for the auto-correlation function: 1--single exponential funciton; 2--Gaussian auto correlation function
	//flag for the auto-correlation function: 3--rotational anisotropic single exponential funciton
	int flag_auto;

	double beta; 
	// the rotational angle of z-axis from the y-axis, unit: degree, for flag_auto=3 and two dimensions
	// the rotational angle of z-axis from the xy-comprehensive axis, unit: degree, for flag_auto=3 and two dimensions
	double beta_1; 
	// the rotational angle of from the x-direction in xy plane, unit: degree, for flag_auto=3 and three dimensions

	// methods
	Flag_SpatialVariability();
	~Flag_SpatialVariability();
};

class Para_SpatialVariability
{
public:
	//No. of variable kinds
	//0-internal friction angle for soil; 1-cohesion for soil; 2-coefficient of permeability
	//3-soil particle size; 4-internal friction angle for fluid; 5-cohesion for fluid
	int no;

	//type of distribution:
	//1--standard normal distribution and 2--log-normal distribution
	int type;

	//parameters
	double mean_value;	//mean value of random
	double sd_value;	//mean value of random variabilities
	double rl_value[3]; //relevant length of problem
	double coe_depth;  //change rate of variable with depth: x=x0+coe_depth*z
	double a[3];		//problem domain [-a, a]

	//methods
	Para_SpatialVariability();
	~Para_SpatialVariability();
};

class spv_cell {//cell class for the randon field generation
public:
	int ntotal;   //total number of particles in the cell
	int parti_id[40];    //particle id in the cell
	double cent_x[3];    //coordinates of cell center

	spv_cell();
	~spv_cell();
};

class clFSpatialVariability_Fun
{
public:
	//generate the spatial-distributed variables
	//function that can be cited in the program
	int SpatialVariables_Generate_2D(Particle *pPar, Para_Fluid *pVFluid, Para_Soil *pParti_ConsPara,
									 const Para_Pro &pPPro, const Flag_SpatialVariability fsv, Para_SpatialVariability *psv, const Cell_Con Cellc);
	int SpatialVariables_Generate_3D(Particle *pPar, Para_Fluid *pVFluid, Para_Soil *pParti_ConsPara,
									 const Para_Pro &pPPro, const Flag_SpatialVariability fsv, Para_SpatialVariability *psv, const Cell_Con Cellc);
	int SpatialVariables_Generate_2D_corr(Particle *pPar, Para_Fluid *pVFluid, Para_Soil *pParti_ConsPara,
									 const Para_Pro &pPPro, const Flag_SpatialVariability fsv, Para_SpatialVariability *psv, const Cell_Con Cellc);
	int SpatialVariables_Generate_3D_corr(Particle *pPar, Para_Fluid *pVFluid, Para_Soil *pParti_ConsPara,
									 const Para_Pro &pPPro, const Flag_SpatialVariability fsv, Para_SpatialVariability *psv, const Cell_Con Cellc);
	int nonstationary_random_2D(Particle* pPar, Para_Fluid* pVFluid, Para_Soil* pParti_ConsPara,
		const Para_Pro& pPPro, const Flag_SpatialVariability fsv, Para_SpatialVariability* psv, const Cell_Con Cellc);
	int nonstationary_random_3D(Particle* pPar, Para_Fluid* pVFluid, Para_Soil* pParti_ConsPara,
		const Para_Pro& pPPro, const Flag_SpatialVariability fsv, Para_SpatialVariability* psv, const Cell_Con Cellc);
	int nonstationary_correlated_random_2D(Particle* pPar, Para_Fluid* pVFluid, Para_Soil* pParti_ConsPara,
		const Para_Pro& pPPro, const Flag_SpatialVariability fsv, Para_SpatialVariability* psv, const Cell_Con Cellc);
	int nonstationary_correlated_random_3D(Particle* pPar, Para_Fluid* pVFluid, Para_Soil* pParti_ConsPara,
		const Para_Pro& pPPro, const Flag_SpatialVariability fsv, Para_SpatialVariability* psv, const Cell_Con Cellc);
	int random_field_input(Particle* pPar, Para_Fluid* pVFluid, Para_Soil* pParti_ConsPara,
		const Para_Pro& pPPro, char* argv);
	int random_variable_generate(Particle* pPar, Para_Fluid* pVFluid, Para_Soil* pParti_ConsPara,
		const Para_Pro& pPPro, const Flag_SpatialVariability fsv, Para_SpatialVariability* psv);
private:
	//solve the value w_i
	void solve_wi(double (*wi)[3] , double a[3], double l[3], int ndim);
	//output the problem parameters for the spatial variability
	void output_Para_SpatialVariability(Flag_SpatialVariability fsv, Para_SpatialVariability *psv);
	//solve the wi for 2D and 3D
	int solve_ramdai(double *eigval, double (*wi)[3], double dirac_r[3],
											   int (*order_id)[3], int ndim);
	//solve H for 2D and 3D
	int solve_h_2D(Particle *pPar, Para_Fluid *pVFluid, Para_Soil *pParti_ConsPara,
				   const Para_Pro &pPPro, double (*wi)[3], double *eigval, double *uu, double x_trans[3], int (*order_id)[3],
				   Para_SpatialVariability psv, double soil_max, double water_max, int type);
	int solve_h_3D(Particle *pPar, Para_Fluid *pVFluid, Para_Soil *pParti_ConsPara,
				   const Para_Pro &pPPro, double (*wi)[3], double *eigval, double *uu, double x_trans[3], int (*order_id)[3],
				   Para_SpatialVariability psv, double soil_max, double water_max, int type);
};

class clFParaOutput_Fun
{
public:
	//output the spatial-distributed variables
	//function that can be cited in the program
	void paramters_output(Particle *pPar, Para_Fluid *pVFluid, Para_Soil *pParti_ConsPara,
						  const Para_Pro &pPPro, char *argv);
};