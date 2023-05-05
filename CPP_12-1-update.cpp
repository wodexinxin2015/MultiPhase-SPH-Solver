/************************************************************************************
Multi-Phase and Parallelized SPH program
--parallelized and optimized for SMP (Symmetric multiprocessing ) architecture.
--Copyright (c) 2018-2023, Weijie ZHANG, GeoHohai, Hohai University.
************************************************************************************/

#include "Class_Functions.h"
#include "Header_Option.h"

/* Update of particle information */
void clUpdate_Fun::update(Particle *pPar, const Para_Pro &pPPro, const Para_Boundary &pPBou,
						  int cn)
{
	int i, k, loop;
	int ntotal = pPPro.ntotal;
	int ndim = pPPro.ndim;
	double dt = pPPro.dt;

	loop = pPPro.l - pPPro.inip;

#pragma omp parallel for schedule(static) private(k)
	for (i = 0; i < ntotal; i++)
	{
		for (k = 0; k < ndim; k++)
		{
			if (pPar[i].type != 0)
			{
				if (loop == 0)
				{
					pPar[i].x[k] = pPar[i].xp[k] + (pPar[i].vx[k] + pPar[i].vxp[k]) * dt  * 0.5;
					pPar[i].ux[k] = 0.0;
				}
				else
				{
					pPar[i].x[k] = pPar[i].xp[k] + (pPar[i].vx[k] + pPar[i].vxp[k]) * dt  * 0.5;
					pPar[i].ux[k] = pPar[i].ux[k] + (pPar[i].vx[k] + pPar[i].vxp[k]) * dt  * 0.5;
				}
				pPar[i].xp[k] = pPar[i].x[k];
				pPar[i].vxp[k] = pPar[i].vx[k];
				pPar[i].axp[k] = pPar[i].ax[k];
				pPar[i].rhop = pPar[i].rho;
			}
			else
			{
				if (pPar[i].matype == 2 || pPar[i].matype == 3 || pPar[i].matype == 4 || pPar[i].matype == 5)
				{
					pPar[i].x[k] = pPar[i].xp[k] + (pPar[i].interfss[k] + pPar[i].interfss[k]) * dt  * 0.5;
					pPar[i].ux[k] = pPar[i].ux[k] + (pPar[i].interfss[k] + pPar[i].interfss[k]) * dt  * 0.5;
					pPar[i].xp[k] = pPar[i].x[k];
					pPar[i].rhop = pPar[i].rho;
					pPar[i].prep = pPar[i].pre;
				}
			}
		}
	}
}