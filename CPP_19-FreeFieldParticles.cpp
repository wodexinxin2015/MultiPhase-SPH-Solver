/************************************************************************************
Multi-Phase and Parallelized SPH program
--parallelized and optimized for SMP (Symmetric multiprocessing ) architecture.
--Copyright (c) 2018-2023, Weijie ZHANG, GeoHohai, Hohai University.
************************************************************************************/

#include <math.h>
#include "Header_Option.h"
#include "Class_Functions.h"

inline void setvelocity(int j, double* vx, double a, double b)
{
	int i;
	for (i = 0; i < 3; i++)
	{
		if (i == j)
			vx[i] = a * vx[i];
		else
			vx[i] = b * vx[i];
	}
}

/*calculating the free-field stress*/
void clFreeFieldBndy_Fun::ff_stress_2d_solid(double *sigxx, double *sigyy, double vx_a[3], double vx_b[3],
											 double e, double poi, double dens, int type)
{
	double c_p, c_s;
	double k, g;

	g = e * 0.5  / (1 + poi);
	k = e * 0.333333333333333  / (1 - 2.0 * poi);
	c_p = sqrt((k + 4.0 * g * 0.333333333333333 ) / dens);
	c_s = sqrt(g / dens);

	if (type == 1)
	{
		*sigxx = -c_p * dens * (vx_a[0] - vx_b[0]);
		*sigyy = -c_s * dens * (vx_a[1] - vx_b[1]);
	}
	else if (type == 2)
	{
		*sigxx = -c_s * dens * (vx_a[0] - vx_b[0]);
		*sigyy = -c_p * dens * (vx_a[1] - vx_b[1]);
	}
}

void clFreeFieldBndy_Fun::ff_stress_3d_solid(double *sigxx, double *sigyy, double *sigzz, double vx_a[3], double vx_b[3],
											 double e, double poi, double dens, int type)
{
	double c_p, c_s;
	double k, g;

	g = e * 0.5 / (1 + poi);
	k = e * 0.333333333333333 / (1 - 2.0 * poi);
	c_p = sqrt((k + 4.0 * g * 0.333333333333333) / dens);
	c_s = sqrt(g / dens);

	if (type == 1)
	{
		*sigxx = -c_p * dens * (vx_a[0] - vx_b[0]);
		*sigyy = -c_s * dens * (vx_a[1] - vx_b[1]);
		*sigzz = -c_s * dens * (vx_a[2] - vx_b[2]);
	}
	else if (type == 2)
	{
		*sigxx = -c_s * dens * (vx_a[0] - vx_b[0]);
		*sigyy = -c_p * dens * (vx_a[1] - vx_b[1]);
		*sigzz = -c_s * dens * (vx_a[2] - vx_b[2]);
	}
	else if (type == 3)
	{
		*sigxx = -c_s * dens * (vx_a[0] - vx_b[0]);
		*sigyy = -c_s * dens * (vx_a[1] - vx_b[1]);
		*sigzz = -c_p * dens * (vx_a[2] - vx_b[2]);
	}
}

/*free field: setting stress and calculating acceleration*/
void clFreeFieldBndy_Fun::free_field_solid(Particle *pPar, cl_bndy_pair* pBndy_pair,
	const Para_Soil pFFBound, Par_Cell *pParCell, const Para_Pro &pPPro, int cn)
{

	int i, j, k, nc, nj;
	double tempup1[3], vxj[3];
	int ntotal = pPPro.ntotal;
	int ndim = pPPro.ndim;
	double sigj[3];
	double e = pFFBound.e;
	double poi = pFFBound.poi;

#pragma omp parallel for schedule(static) private(j, k, nc, nj, tempup1, sigj, vxj)
	for (i = 0; i < ntotal; i++)
	{
		if (pPar[i].type == 2)
		{
			for (k = 0; k < 3; k++)
			{
				tempup1[k] = 0.0;
				vxj[k] = 0.0;
			}
			//stress calculation for free boundary particles
			nc = pParCell[i].ninflu;
			for (nj = 1; nj <= nc; nj++)
			{
				j = pParCell[i].influ[nj];
				if (pPar[j].type == 0 && pPar[j].matype == 5)
				{ //Solid phase
					//initialization of variables
					sigj[0] = 0.0;
					sigj[1] = 0.0;
					sigj[2] = 0.0;

					//setting of far field boundary virtual velocity
					vxj[0] = pPar[j].interfss[0];
					vxj[1] = pPar[j].interfss[1];	
					vxj[2] = pPar[j].interfss[2];

					//calculating stress
					if (ndim == 2)
					{
						ff_stress_2d_solid(&sigj[0], &sigj[1], pPar[i].vx, vxj, e,
										   poi, pPar[j].rho, pPar[j].permtype);
					}
					else if (ndim == 3)
					{
						ff_stress_3d_solid(&sigj[0], &sigj[1], &sigj[2], pPar[i].vx, vxj, e,
										   poi, pPar[j].rho, pPar[j].permtype);
					}
					//updating acceleration
					tempup1[0] = tempup1[0] + pPar[j].mass / pPar[j].rho * sigj[0] / pPar[i].rho * pParCell[i].wij[nj][0];
					tempup1[1] = tempup1[1] + pPar[j].mass / pPar[j].rho * sigj[1] / pPar[i].rho * pParCell[i].wij[nj][1];
					tempup1[2] = tempup1[2] + pPar[j].mass / pPar[j].rho * sigj[2] / pPar[i].rho * pParCell[i].wij[nj][2];
				}
			}
			//updating acceleration
			for (k = 0; k < ndim; k++)
			{
				pPar[i].ax[k] = pPar[i].ax[k] + tempup1[k];
			}
		}
	}
}