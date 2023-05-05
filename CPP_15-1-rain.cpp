/************************************************************************************
Multi-Phase and Parallelized SPH program
--parallelized and optimized for SMP (Symmetric multiprocessing ) architecture.
--Copyright (c) 2018-2023, Weijie ZHANG, GeoHohai, Hohai University.
************************************************************************************/

#include <math.h>
#include "Class_Functions.h"
#include "Header_Option.h"

void clRain_Fun::rainini(const Para_Rain &pPRain, const Para_Pro &pPPro, int *lp, double *vrain)
{
	double temp;
	int nrain = pPRain.nrain;
	int nhor = pPRain.nhor;
	double vel = pPRain.drop_velocity;

	if (pPRain.rtype == 1)
	{ //customized rainfall
		temp = nrain / nhor + 1.0;
		*lp = (int)((pPPro.loop - pPPro.inip) / temp);
		*vrain = vel;
	}
	if (pPRain.rtype == 2)
	{ //water particle deletion
		temp = fabs(pPRain.x_max - pPRain.x_min) / pPPro.dr + 1.0;
		*lp = (int)((pPPro.loop - pPPro.inip) / temp);
	}
	if (pPRain.rtype == 3)
	{ //water particle deletion
		temp = fabs(pPRain.x_max - pPRain.x_min) / pPPro.dr + 1.0;
		*lp = (int)((pPPro.loop - pPPro.inip) / temp);
	}
	if (pPRain.rtype == 4)
	{ //cement injection of tunnel
		temp = nrain / nhor + 1.0;
		*lp = (int)((pPPro.loop - pPPro.inip) / temp);
		*vrain = 2.0;
	}
}

//produce rainfall particles
void clRain_Fun::rainpro(Particle *pVelPar, Para_Rain *pPRain, Para_Pro *pPPro, const Para_Fluid *pVFluid,
						 Para_GF pgf, double vrain)
{

	int k = pPPro->ntotal - pPRain->nrain;
	double temp_xp, press_factor;
	int ndim = pPPro->ndim;

	press_factor = pPRain->press_factor;

	if (pPRain->nrain > 0)
	{
		do
		{
			pPRain->nrain = pPRain->nrain - 1;
			pPPro->nwater = pPPro->nwater + 1;

			temp_xp = pVelPar[k].xp[ndim - 1];
			pVelPar[k].xp[ndim - 1] = 0.95 * pVelPar[k].xp[ndim - 1];
			pVelPar[k].type = 1;
			pVelPar[k].matype = 0;
			pVelPar[k].etype = 7;
			pVelPar[k].rhop = pVFluid[k].porosity * pVFluid[k].dens;
			pVelPar[k].hl = pPPro->dr;
			pVelPar[k].massini(ndim);
			pVelPar[k].vxp[0] = 0.0;
			pVelPar[k].vxp[1] = 0.0;
			pVelPar[k].vxp[2] = 0.0;
			pVelPar[k].vxp[ndim - 1] = vrain;
			pVelPar[k].vx[0] = 0.0;
			pVelPar[k].vx[1] = 0.0;
			pVelPar[k].vx[2] = 0.0;
			pVelPar[k].vxp[ndim - 1] = vrain;
			if (ndim == 2)
				pVelPar[k].prep = pVelPar[k].rhop * pPPro->dr * (-pgf.gy) * press_factor;
			else if (ndim == 3)
				pVelPar[k].prep = pVelPar[k].rhop * pPPro->dr * (-pgf.gz) * press_factor;
			pVelPar[k].satu = 0.99;

			k++;
		} while (fabs(pVelPar[k].xp[ndim - 1] - temp_xp) < 0.1 * pPPro->dr && k <= pPPro->ntotal - 1);
	}
}

//produce cement particles
void clRain_Fun::cementpro(Particle *pVelPar, Para_Rain *pPRain, Para_Pro *pPPro, const Para_Fluid *pVFluid,
						   Para_GF pgf, double vrain, double cement_base[3])
{

	int k = pPPro->ntotal - pPRain->nrain;
	int ndim = pPPro->ndim;
	int i, c_flag;
	double drp;
	double press_factor = pPRain->press_factor;
	c_flag = pPRain->cement_flag;

	if (pPRain->nrain > 0)
	{
		do
		{
			pPRain->nrain = pPRain->nrain - 1;
			pPPro->nwater = pPPro->nwater + 1;

			drp = 0.0;
			for (i = 0; i < ndim; i++)
			{
				drp = drp + (pVelPar[k].xp[i] - cement_base[i]) * (pVelPar[k].xp[i] - cement_base[i]);
			}
			drp = sqrt(drp);

			pVelPar[k].xp[ndim - 1] = pVelPar[k].xp[ndim - 1];
			pVelPar[k].type = 1;
			pVelPar[k].matype = 0;
			pVelPar[k].rhop = pVFluid[k].porosity * pVFluid[k].dens;
			pVelPar[k].hl = pPPro->dr;
			pVelPar[k].massini(ndim);
			pVelPar[k].vxp[0] = vrain * (pVelPar[k].xp[0] - cement_base[0]) / drp;
			pVelPar[k].vxp[1] = vrain * (pVelPar[k].xp[1] - cement_base[1]) / drp;
			pVelPar[k].vxp[2] = vrain * (pVelPar[k].xp[2] - cement_base[2]) / drp;
			if (ndim == 2)
				pVelPar[k].prep = pVelPar[k].rhop * pPPro->dr * (-pgf.gy) * press_factor;
			else if (ndim == 3)
				pVelPar[k].prep = pVelPar[k].rhop * pPPro->dr * (-pgf.gz) * press_factor;
			pVelPar[k].satu = 0.99;

			k++;

		} while (pVelPar[k].matype == c_flag && k < pPPro->ntotal - 1);
	}
	pPRain->cement_flag += 1;
}

//calculating center of cement
void clRain_Fun::center_cement(Particle *pVelPar, Para_Pro &pPPro, double *x)
{
	int i;
	double max_v[3], min_v[3];

	max_v[0] = -1000.0;
	max_v[1] = -1000.0;
	max_v[2] = -1000.0;
	min_v[0] = 1000.0;
	min_v[1] = 1000.0;
	min_v[2] = 1000.0;

	for (i = 0; i < pPPro.ntotal; i++)
	{
		if (pVelPar[i].type == 7)
		{
			max_v[0] = fmax(max_v[0], pVelPar[i].xp[0]);
			max_v[1] = fmax(max_v[1], pVelPar[i].xp[1]);
			max_v[2] = fmax(max_v[2], pVelPar[i].xp[2]);
			min_v[0] = fmin(min_v[0], pVelPar[i].xp[0]);
			min_v[1] = fmin(min_v[1], pVelPar[i].xp[1]);
			min_v[2] = fmin(min_v[2], pVelPar[i].xp[2]);
		}
	}

	x[0] = (max_v[0] + min_v[0]) * 0.5;
	x[1] = (max_v[1] + min_v[1]) * 0.5;
	x[2] = (max_v[2] + min_v[2]) * 0.5;
}