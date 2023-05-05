/************************************************************************************
Multi-Phase and Parallelized SPH program
--parallelized and optimized for SMP (Symmetric multiprocessing ) architecture.
--Copyright (c) 2018-2023, Weijie ZHANG, GeoHohai, Hohai University.
************************************************************************************/

#include "Class_Functions.h"
#include "Header_Option.h"

/* calculation and normalization of the density */
void clDensity_Fun::density(Particle* pPar, Par_Cell* pParCell, Para_Fluid* pVFluid,
	cl_bndy_pair* pBndy_Pair, clVar_Boundary* pParti_VariBndy,
	const Para_Pro& pPPro, int bndytype)
{

	int i, j, k, nc, nj, ntotal, ndim;
	double drhot, tempup[3], tempdown[3], dt, vxj[3];

	ntotal = pPPro.ntotal;
	ndim = pPPro.ndim;
	dt = pPPro.dt;

	/*calculating the porosity for the liquid phase*/
	int l = pPPro.l;
#pragma omp parallel for schedule(static) private(j, nc, nj, tempup, tempdown)
	for (i = 0; i < ntotal; i++)
	{

		tempup[0] = 0.0;
		tempdown[0] = 0.0;
		nc = pParCell[i].ninflu;

		for (nj = 1; nj <= nc; nj++)
		{
			j = pParCell[i].influ[nj];
			if (pPar[i].type != pPar[j].type && pPar[i].type == 1 && pPar[j].type == 2)
			{
				tempup[0] = tempup[0] + pPar[j].mass * pPar[j].porosity / pPar[j].rho * pParCell[i].wij[nj][3];
				tempdown[0] = tempdown[0] + pPar[j].mass / pPar[j].rho * pParCell[i].wij[nj][3];
			}
		}

		if (pPar[i].type == 1 && fabs(tempdown[0]) > 1.0e-7)
		{
			pPar[i].poro_cop = pPar[i].porosity;
			pPar[i].porosity = tempup[0] / tempdown[0];
			if (l == 1)
				pPar[i].poro_cop = pPar[i].porosity;
		}
		else if (pPar[i].type == 1 && fabs(tempdown[0]) <= 1.0e-7)
		{
			pPar[i].poro_cop = pPar[i].porosity;
		}
	}

	/*calculation of the density according to continuity equation for stress particles*/
#pragma omp parallel for schedule(static) private(j, k, nc, nj, drhot, tempup, tempdown, vxj)
	for (i = 0; i < ntotal; i++)
	{
		//initialization of local variables
		drhot = 0.0;
		for (k = 0; k < ndim; k++) {
			tempup[k] = 0.0;
			tempdown[k] = 0.0;
		}

		nc = pParCell[i].ninflu;
		for (nj = 1; nj <= nc; nj++)
		{
			j = pParCell[i].influ[nj];
			if (pPar[i].type == pPar[j].type && pPar[i].type != 7 && pPar[i].type != 0)
			{
				//density rate
				if (pPar[i].type == 2 || pPar[i].type == 4) {
					for (k = 0; k < ndim; k++)
					{
						tempup[k] = tempup[k] + pPar[j].mass * (pPar[i].vxp[k] - pPar[j].vxp[k]) * pParCell[i].wij[nj][k];
						tempdown[k] = tempdown[k] + pPar[j].mass * (pPar[i].xp[k] - pPar[j].xp[k]) * pParCell[i].wij[nj][k];
					}
				}
				if (pPar[i].type == 1 || pPar[i].type == 3) {
					for (k = 0; k < ndim; k++)
					{
						tempup[k] = tempup[k] + pPar[j].mass * (pPar[i].vxp[k] - pPar[j].vxp[k]) * pParCell[i].wij[nj][k];
					}
				}
				
			}
			else if ((pPar[i].type == 1 || pPar[i].type == 2 || pPar[i].type == 3 || pPar[i].type == 4) && pPar[j].type == 0) {
				//settting virtual velocity
				vxj[0] = 0.0;
				vxj[1] = 0.0;
				vxj[2] = 0.0;
				if (bndytype == 1) {
					vxj[0] = pPar[i].vxp[0];
					vxj[1] = pPar[i].vxp[1];
					vxj[2] = pPar[i].vxp[2];
				}
				else if (bndytype == 2 || bndytype == 3 || bndytype == 4 || bndytype == 13 || bndytype == 14) {
					if (pPar[i].type == 1) {
						vxj[0] = pParti_VariBndy[j].vx_water[0];
						vxj[1] = pParti_VariBndy[j].vx_water[1];
						vxj[2] = pParti_VariBndy[j].vx_water[2];
					}
					else if (pPar[i].type == 2) {
						vxj[0] = pParti_VariBndy[j].vx_soil[0];
						vxj[1] = pParti_VariBndy[j].vx_soil[1];
						vxj[2] = pParti_VariBndy[j].vx_soil[2];
					}
					else if (pPar[i].type == 3) {
						vxj[0] = pParti_VariBndy[j].vx_air[0];
						vxj[1] = pParti_VariBndy[j].vx_air[1];
						vxj[2] = pParti_VariBndy[j].vx_air[2];
					}
					else if (pPar[i].type == 4) {
						vxj[0] = pParti_VariBndy[j].vx_struct[0];
						vxj[1] = pParti_VariBndy[j].vx_struct[1];
						vxj[2] = pParti_VariBndy[j].vx_struct[2];
					}
				}
				//density rate
				if (pPar[i].type == 2 || pPar[i].type == 4) {
					for (k = 0; k < ndim; k++)
					{
						tempup[k] = tempup[k] + pPar[j].mass * (pPar[i].vxp[k] - vxj[k]) * pParCell[i].wij[nj][k];
						tempdown[k] = tempdown[k] + pPar[j].mass * (pPar[i].xp[k] - pPar[j].xp[k]) * pParCell[i].wij[nj][k];
					}
				}
				if (pPar[i].type == 1 || pPar[i].type == 3) {
					for (k = 0; k < ndim; k++)
					{
						tempup[k] = tempup[k] + pPar[j].mass * (pPar[i].vxp[k] - vxj[k]) * pParCell[i].wij[nj][k];
					}
				}
			}
		}

		if (pPar[i].type != 7) {
			if (pPar[i].type == 2 || pPar[i].type == 4) {
				for (k = 0; k < ndim; k++) {
					if (fabs(tempdown[k]) > 1.0e-7)
						drhot -= tempup[k] / tempdown[k];
					else
						drhot += tempup[k];
				}
			}
			if (pPar[i].type == 1 || pPar[i].type == 3) {
				for (k = 0; k < ndim; k++) {
					drhot += tempup[k];
				}
			}

			if (pPar[i].type == 1) {
				if (fabs(pPar[i].porosity) > 0.0001) {
					pPar[i].rho = pPar[i].rhop + drhot * dt;
					pPar[i].rho = pPar[i].rho * pPar[i].poro_cop / pPar[i].porosity; //considering the porosity effect
				}
				else {
					pPar[i].rho = pPar[i].rhop + drhot * dt;
				}
			}
			else {
				pPar[i].rho = pPar[i].rhop + drhot * dt;
			}
		}
	}

	/*Improved density of normalization or density interpolation for velocity particles*/
#pragma omp parallel for schedule(static) private(j, nc, nj, tempup, tempdown)
	for (i = 0; i < ntotal; i++)
	{

		tempup[0] = 0.0;
		tempdown[0] = 0.0;
		nc = pParCell[i].ninflu;

		for (nj = 1; nj <= nc; nj++)
		{
			j = pParCell[i].influ[nj];
			if (pPar[i].type == pPar[j].type)
			{
				tempup[0] = tempup[0] + pPar[j].mass * pParCell[i].wij[nj][3];
				tempdown[0] = tempdown[0] + pPar[j].mass / pPar[j].rho * pParCell[i].wij[nj][3];
			}
		}

		if (pPar[i].type != 7 && fabs(tempdown[0]) > 1.0e-7)
			pPar[i].rho = tempup[0] / tempdown[0];
	}
}
