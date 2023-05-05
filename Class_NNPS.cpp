/************************************************************************************
Multi-Phase and Parallelized SPH program
--parallelized and optimized for SMP (Symmetric multiprocessing ) architecture.
--Copyright (c) 2018-2023, Weijie ZHANG, GeoHohai, Hohai University.
************************************************************************************/

#include <stdlib.h>
#include "Class_NNPS.h"

Cell_Con::Cell_Con()
{
	cxm = 1;
	cym = 1;
	czm = 1;
	ctotal = 1;

	xmin = 0.0;
	xmax = 0.0;
	ymin = 0.0;
	ymax = 0.0;
	zmin = 0.0;
	zmax = 0.0;

	xmin_nb = 0.0;
	xmax_nb = 0.0;
	ymin_nb = 0.0;
	ymax_nb = 0.0;
	zmin_nb = 0.0;
	zmax_nb = 0.0;

	for (int i = 0; i < 3; i++)
	{
		water_min[i] = 100000.000;
		water_max[i] = -100000.000;
		soil_min[i] = 100000.000;
		soil_max[i] = -100000.000;
		air_min[i] = 100000.000;
		air_max[i] = -100000.000;
		stru_min[i] = 100000.000;
		stru_max[i] = -100000.000;
	}
}

Cell_Con::~Cell_Con()
{
}

parti_cellid::parti_cellid()
{
	cell_id = 10000000;
	parti_id = 10000000;
}

parti_cellid::~parti_cellid()
{
}

cell_link::cell_link()
{
	for (int i = 0; i < 27; i++)
		nncell[i] = -1;
}

cell_link::~cell_link()
{
}

cell_info::cell_info()
{
	start = 2000000;
	end = -2000000;
}

cell_info::~cell_info()
{
}

Par_Cell::Par_Cell()
{
	cell_id = -1;
	ninflu = 0;
	for (int i = 0; i < 100; i++)
	{
		influ[i] = 0;
		for (int j = 0; j < 4; j++)
			wij[i][j] = 0.0;
	}
}

Par_Cell::~Par_Cell()
{
}

void Par_Cell::initial()
{
	ninflu = 0;
}