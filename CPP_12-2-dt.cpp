/************************************************************************************
Multi-Phase and Parallelized SPH program
--parallelized and optimized for SMP (Symmetric multiprocessing ) architecture.
--Copyright (c) 2018-2023, Weijie ZHANG, GeoHohai, Hohai University.
************************************************************************************/

#include <cmath>
#include <algorithm>
#include "Class_Functions.h"
#include "Header_Option.h"

/* calculation of the time step using the CFL condition */
void clUpdate_Fun::timestep(Para_Pro * pPPro)
{ 
	double c_cour = 0.25;

	pPPro->tt = pPPro->dt;
	pPPro->dt = c_cour * pPPro->dr / 600.0;
	pPPro->loop = int(pPPro->tt / pPPro->dt) + 1;
}
