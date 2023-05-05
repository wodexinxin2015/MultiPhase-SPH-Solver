/************************************************************************************
Multi-Phase and Parallelized SPH program
--parallelized and optimized for SMP (Symmetric multiprocessing ) architecture.
--Copyright (c) 2018-2023, Weijie ZHANG, GeoHohai, Hohai University.
************************************************************************************/

#include "Class_spatial_variability.h"

Flag_SpatialVariability::Flag_SpatialVariability()
{
	flag_fai = 0;
	flag_c = 0;
	flag_cop = 0;
	flag_ds = 0;

	flag_faiw = 0;
	flag_cw = 0;
	
	cell_ndx = 1.0;
	coe_frict_cohesion = 0.0;
	flag_auto = 1;

	beta = 0;
	beta_1 = 0;
}

Flag_SpatialVariability::~Flag_SpatialVariability()
{
}

Para_SpatialVariability::Para_SpatialVariability()
{
	no = -1;
	type = 1;

	mean_value = 0.0;
	sd_value = 0.0;
	rl_value[0] = 0.0;
	rl_value[1] = 0.0;
	rl_value[2] = 0.0;
	coe_depth = 0.0;
	a[0] = 0.0;
	a[1] = 0.0;
	a[2] = 0.0;
}

Para_SpatialVariability::~Para_SpatialVariability()
{
}

spv_cell::spv_cell() {
	ntotal = 0;
	cent_x[0] = 0.0;
	cent_x[1] = 0.0;
	cent_x[2] = 0.0;
}

spv_cell::~spv_cell() {

}