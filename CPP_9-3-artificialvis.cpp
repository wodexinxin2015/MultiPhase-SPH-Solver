/************************************************************************************
Multi-Phase and Parallelized SPH program
--parallelized and optimized for SMP (Symmetric multiprocessing ) architecture.
--Copyright (c) 2018-2023, Weijie ZHANG, GeoHohai, Hohai University.
************************************************************************************/

#include <cmath>
#include "Header_Parameters.h"
#include "Header_Option.h"
#include "Class_Functions.h"

void clAcce_Fun::artificial_viscosity(Particle *pPar, Par_Cell *pParCell, const Para_Soil *pParti_ConsPara,
									  const Para_Fluid *pVFluid, const Para_Pro &pPPro, int cn)
{

	int i, j, k, nc, nj;
	int ntotal, ndim;
	double inprod, dvx[3], dx[3], paix, c;
	double tempup[3], dst, myu, hl, rou;
	double a, b;

	ntotal = pPPro.ntotal;
	ndim = pPPro.ndim;

#pragma omp parallel for schedule(static) private(j, k, nc, nj, inprod, dvx, dx, paix, tempup, dst, myu, hl, rou, c, a, b)
	for (i = 0; i < ntotal; i++)
	{

		for (k = 0; k < 3; k++)
		{
			tempup[k] = 0.0;
		}

		nc = pParCell[i].ninflu;
		for (nj = 1; nj <= nc; nj++)
		{
			j = pParCell[i].influ[nj];

			if (pPar[i].type == pPar[j].type && pPar[i].type != 7)
			{
				a = 0.0;
				b = 0.0;
				c = 0.0;
				inprod = 0.0;
				dst = 0.0;
				myu = 0.0;
				paix = 0.0;

				for (k = 0; k < ndim; k++)
				{
					dvx[k] = pPar[i].vxp[k] - pPar[j].vxp[k];
					dx[k] = pPar[i].xp[k] - pPar[j].xp[k];
					dst = dst + dx[k] * dx[k];
					inprod = inprod + dvx[k] * dx[k];
				}

				hl = (pPar[j].hl + pPar[i].hl) * 0.5;
				rou = (pPar[j].rho + pPar[i].rho) * 0.5;

				if (pPar[i].type == 1)
				{
					c = pVFluid[i].vpu;
					a = pPPro.art_vis[0][0];
					b = pPPro.art_vis[0][1];
				}
				if (pPar[i].type == 2)
				{
					c = pParti_ConsPara[i].vpu;
					a = pPPro.art_vis[1][0];
					b = pPPro.art_vis[1][1];
				}
				if (pPar[i].type == 3)
				{
					c = pVFluid[i].vpu;
					a = pPPro.art_vis[2][0];
					b = pPPro.art_vis[2][1];
				}

				if (inprod < 0.0)
				{
					dst = sqrt(dst);
					myu = hl * inprod / (dst * dst + 0.01 * hl * hl);
					paix = (-a * c * myu + b * myu * myu) / rou;
				}
				else
					paix = 0.0;

				for (k = 0; k < ndim; k++)
					tempup[k] = tempup[k] + pPar[j].mass * paix * pParCell[i].wij[nj][k];
			}
		}

		if (pPar[i].type != 7)
		{
			for (k = 0; k < ndim; k++)
				pPar[i].ax[k] = pPar[i].ax[k] - tempup[k];
		}
	}
}