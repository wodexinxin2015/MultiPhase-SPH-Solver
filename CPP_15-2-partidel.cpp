/************************************************************************************
Multi-Phase and Parallelized SPH program
--parallelized and optimized for SMP (Symmetric multiprocessing ) architecture.
--Copyright (c) 2018-2023, Weijie ZHANG, GeoHohai, Hohai University.
************************************************************************************/

#include <math.h>
#include "Class_Functions.h"
#include "Header_Option.h"

inline void setvelocity(double *vx)
{
	vx[0] = 0.0;
	vx[1] = 0.0;
	vx[2] = 0.0;
}

/*delete particles*/
void clRain_Fun::partidel(Particle *pVelPar, Para_Rain *pPRain,
						  Para_Pro *pPPro, const Cell_Con &pCellc, int cn)
{
	int i;
	double dx1, dx2, dy1;

	dx1 = 2.0 * pPPro->dr;
	dx2 = 2.0 * pPPro->dr;

#pragma omp parallel for schedule(static) private(dx1, dx2, dy1)
	for (i = 0; i < pPPro->ntotal; i++)
	{

		dx1 = fabs(pVelPar[i].xp[0] - pCellc.xmax - 7.0 * pPPro->dr);
		//dx1 = fabs(pVelPar[i].xp[0] - pCellc.xmax - 3.0*pPPro->dr);
		dx2 = fabs(pVelPar[i].xp[0] - pCellc.xmin + 7.0 * pPPro->dr);
		//dx2 = fabs(pVelPar[i].xp[0] - pCellc.xmin + 3.0*pPPro->dr);

		dy1 = pVelPar[i].xp[1] - 0.12;

		if (pVelPar[i].type == 1 &&
			(dx1 < pPPro->dr || dx2 < pPPro->dr) && dy1 > 0.0)
		{

			pVelPar[i].type = 7;
			pVelPar[i].vxp[0] = 0.0;
			pVelPar[i].vxp[1] = 0.0;
			pVelPar[i].vxp[2] = 0.0;
			pVelPar[i].vx[0] = 0.0;
			pVelPar[i].vx[1] = 0.0;
			pVelPar[i].vx[2] = 0.0;
			pVelPar[i].ax[0] = 0.0;
			pVelPar[i].ax[1] = 0.0;
			pVelPar[i].ax[2] = 0.0;

#pragma omp critical
			{
				pPPro->nwater = pPPro->nwater - 1;
			}
		}
	}
}

/*delete particles*/
void clRain_Fun::waterdel_min(Particle *pVelPar, Para_Rain *pPRain, Para_Pro *pPPro)
{
	int ndim = pPPro->ndim;
	int i;
	double x_cord, x_min;

	x_min = pPRain->x_min + 3.0 * pPPro->dr;

	for (i = 0; i < pPPro->ntotal; i++)
	{
		if (pVelPar[i].type == 1 && pVelPar[i].etype == 7)
		{
			x_cord = pVelPar[i].xp[ndim - 2];
			if (x_cord <= x_min)
			{
				pVelPar[i].type = 7;
				pVelPar[i].etype = 1;

				pVelPar[i].vxp[0] = 0.0;
				pVelPar[i].vxp[1] = 0.0;
				pVelPar[i].vxp[2] = 0.0;
				pVelPar[i].vx[0] = 0.0;
				pVelPar[i].vx[1] = 0.0;
				pVelPar[i].vx[2] = 0.0;
				pVelPar[i].ax[0] = 0.0;
				pVelPar[i].ax[1] = 0.0;
				pVelPar[i].ax[2] = 0.0;

				pPRain->nrain = pPRain->nrain + 1;
				pPPro->nwater = pPPro->nwater - 1;
			}
		}
	}
}

/*delete particles*/
void clRain_Fun::waterdel_max(Particle *pVelPar, Para_Rain *pPRain, Para_Pro *pPPro)
{
	int ndim = pPPro->ndim;
	int i;
	double x_cord, x_max;

	x_max = pPRain->x_max - 3.0 * pPPro->dr;

	for (i = 0; i < pPPro->ntotal; i++)
	{
		if (pVelPar[i].type == 1 && pVelPar[i].etype == 7)
		{
			x_cord = pVelPar[i].xp[ndim - 2];
			if (x_cord >= x_max)
			{
				pVelPar[i].type = 7;
				pVelPar[i].etype = 1;

				pVelPar[i].vxp[0] = 0.0;
				pVelPar[i].vxp[1] = 0.0;
				pVelPar[i].vxp[2] = 0.0;
				pVelPar[i].vx[0] = 0.0;
				pVelPar[i].vx[1] = 0.0;
				pVelPar[i].vx[2] = 0.0;
				pVelPar[i].ax[0] = 0.0;
				pVelPar[i].ax[1] = 0.0;
				pVelPar[i].ax[2] = 0.0;

				pPRain->nrain = pPRain->nrain + 1;
				pPPro->nwater = pPPro->nwater - 1;
			}
		}
	}
}