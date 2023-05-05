/************************************************************************************
Multi-Phase and Parallelized SPH program
--parallelized and optimized for SMP (Symmetric multiprocessing ) architecture.
--Copyright (c) 2018-2023, Weijie ZHANG, GeoHohai, Hohai University.
************************************************************************************/
#include <cmath>
#include "Header_Option.h"
#include "Header_Parameters.h"
#include "Class_Functions.h"

/* the stress derivations for water particles */
void clAcce_Fun::fluid_acceleration_water(Particle *pPar, Par_Cell *pParCell, const Para_Fluid *pVFluid,
										  const Para_Pro &pPPro, const Para_GF &pPGf, int cn)
{

	int i, j, k, nc, nj;
	int ntotal;
	double tempup[3][3], gcc[3], porosity;
	double asigi[6], asigj[6];

	ntotal = pPPro.ntotal;

	gcc[0] = pPGf.gx;
	gcc[1] = pPGf.gy;
	gcc[2] = pPGf.gz;

#pragma omp parallel for schedule(static) private(j, k, nc, nj, tempup, asigi, asigj, porosity)
	for (i = 0; i < ntotal; i++)
	{
		if (pPar[i].type == 1)
		{
			porosity = pVFluid[i].porosity;

			for (k = 0; k < 3; k++)
			{
				tempup[k][0] = 0.0;
				tempup[k][1] = 0.0;
				tempup[k][2] = 0.0;
				pPar[i].ax[k] = gcc[k];
			}

			nc = pParCell[i].ninflu;
			for (nj = 1; nj <= nc; nj++)
			{
				j = pParCell[i].influ[nj];

				if (pPar[i].type == pPar[j].type)
				{
					for (k = 0; k < 6; k++)
					{
						asigi[k] = porosity * pPar[i].sig[k];
						asigj[k] = porosity * pPar[j].sig[k];
					}

					tempup[0][0] = tempup[0][0] + pPar[j].mass * (asigj[0] / pPar[j].rho / pPar[j].rho + asigi[0] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][0];
					tempup[0][1] = tempup[0][1] + pPar[j].mass * (asigj[3] / pPar[j].rho / pPar[j].rho + asigi[3] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][1];
					tempup[0][2] = tempup[0][2] + pPar[j].mass * (asigj[5] / pPar[j].rho / pPar[j].rho + asigi[5] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][2];

					tempup[1][0] = tempup[1][0] + pPar[j].mass * (asigj[3] / pPar[j].rho / pPar[j].rho + asigi[3] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][0];
					tempup[1][1] = tempup[1][1] + pPar[j].mass * (asigj[1] / pPar[j].rho / pPar[j].rho + asigi[1] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][1];
					tempup[1][2] = tempup[1][2] + pPar[j].mass * (asigj[4] / pPar[j].rho / pPar[j].rho + asigi[4] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][2];

					tempup[2][0] = tempup[2][0] + pPar[j].mass * (asigj[5] / pPar[j].rho / pPar[j].rho + asigi[5] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][0];
					tempup[2][1] = tempup[2][1] + pPar[j].mass * (asigj[4] / pPar[j].rho / pPar[j].rho + asigi[4] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][1];
					tempup[2][2] = tempup[2][2] + pPar[j].mass * (asigj[2] / pPar[j].rho / pPar[j].rho + asigi[2] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][2];

				}
				else if (pPar[j].type == 3)
				{
					for (k = 0; k < 6; k++)
					{
						asigi[k] = porosity * pPar[i].sig[k];
						asigj[k] = porosity * pPar[j].sig[k];
					}

					tempup[0][0] = tempup[0][0] + pPar[j].mass * (asigj[0] / pPar[j].rho / pPar[j].rho + asigi[0] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][0];
					tempup[0][1] = tempup[0][1] + pPar[j].mass * (asigj[3] / pPar[j].rho / pPar[j].rho + asigi[3] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][1];
					tempup[0][2] = tempup[0][2] + pPar[j].mass * (asigj[5] / pPar[j].rho / pPar[j].rho + asigi[5] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][2];

					tempup[1][0] = tempup[1][0] + pPar[j].mass * (asigj[3] / pPar[j].rho / pPar[j].rho + asigi[3] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][0];
					tempup[1][1] = tempup[1][1] + pPar[j].mass * (asigj[1] / pPar[j].rho / pPar[j].rho + asigi[1] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][1];
					tempup[1][2] = tempup[1][2] + pPar[j].mass * (asigj[4] / pPar[j].rho / pPar[j].rho + asigi[4] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][2];

					tempup[2][0] = tempup[2][0] + pPar[j].mass * (asigj[5] / pPar[j].rho / pPar[j].rho + asigi[5] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][0];
					tempup[2][1] = tempup[2][1] + pPar[j].mass * (asigj[4] / pPar[j].rho / pPar[j].rho + asigi[4] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][1];
					tempup[2][2] = tempup[2][2] + pPar[j].mass * (asigj[2] / pPar[j].rho / pPar[j].rho + asigi[2] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][2];
				}
			}

			pPar[i].ax[0] += tempup[0][0] + tempup[0][1] + tempup[0][2];
			pPar[i].ax[1] += tempup[1][0] + tempup[1][1] + tempup[1][2];
			pPar[i].ax[2] += tempup[2][0] + tempup[2][1] + tempup[2][2];
		}
	}
}

/* the stress derivations for air particles*/
void clAcce_Fun::fluid_acceleration_air(Particle *pPar, Par_Cell *pParCell, const Para_Fluid *pVFluid,
										const Para_Pro &pPPro, const Para_GF &pPGf, int cn)
{

	int i, j, k, nc, nj;
	int ntotal;
	double tempup[3][3], gcc[3], porosity;
	double asigi[6], asigj[6];

	ntotal = pPPro.ntotal;

	gcc[0] = pPGf.gx;
	gcc[1] = pPGf.gy;
	gcc[2] = pPGf.gz;

#pragma omp parallel for schedule(static) private(j, k, nc, nj, tempup, asigi, asigj, porosity)
	for (i = 0; i < ntotal; i++)
	{
		if (pPar[i].type == 3)
		{
			porosity = pVFluid[i].porosity;

			for (k = 0; k < 3; k++)
			{
				tempup[k][0] = 0.0;
				tempup[k][1] = 0.0;
				tempup[k][2] = 0.0;
				pPar[i].ax[k] = gcc[k];
			}

			nc = pParCell[i].ninflu;
			for (nj = 1; nj <= nc; nj++)
			{
				j = pParCell[i].influ[nj];

				if (pPar[i].type == pPar[j].type)
				{
					for (k = 0; k < 6; k++)
					{
						asigi[k] = porosity * pPar[i].sig[k];
						asigj[k] = porosity * pPar[j].sig[k];
					}

					tempup[0][0] = tempup[0][0] + pPar[j].mass * (asigj[0] / pPar[j].rho / pPar[j].rho + asigi[0] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][0];
					tempup[0][1] = tempup[0][1] + pPar[j].mass * (asigj[3] / pPar[j].rho / pPar[j].rho + asigi[3] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][1];
					tempup[0][2] = tempup[0][2] + pPar[j].mass * (asigj[5] / pPar[j].rho / pPar[j].rho + asigi[5] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][2];

					tempup[1][0] = tempup[1][0] + pPar[j].mass * (asigj[3] / pPar[j].rho / pPar[j].rho + asigi[3] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][0];
					tempup[1][1] = tempup[1][1] + pPar[j].mass * (asigj[1] / pPar[j].rho / pPar[j].rho + asigi[1] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][1];
					tempup[1][2] = tempup[1][2] + pPar[j].mass * (asigj[4] / pPar[j].rho / pPar[j].rho + asigi[4] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][2];

					tempup[2][0] = tempup[2][0] + pPar[j].mass * (asigj[5] / pPar[j].rho / pPar[j].rho + asigi[5] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][0];
					tempup[2][1] = tempup[2][1] + pPar[j].mass * (asigj[4] / pPar[j].rho / pPar[j].rho + asigi[4] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][1];
					tempup[2][2] = tempup[2][2] + pPar[j].mass * (asigj[2] / pPar[j].rho / pPar[j].rho + asigi[2] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][2];
				}
				else if (pPar[j].type == 1)
				{
					for (k = 0; k < 6; k++)
					{
						asigi[k] = porosity * pPar[i].sig[k];
						asigj[k] = porosity * pPar[j].sig[k];
					}

					tempup[0][0] = tempup[0][0] + pPar[j].mass * (asigj[0] / pPar[j].rho / pPar[j].rho + asigi[0] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][0];
					tempup[0][1] = tempup[0][1] + pPar[j].mass * (asigj[3] / pPar[j].rho / pPar[j].rho + asigi[3] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][1];
					tempup[0][2] = tempup[0][2] + pPar[j].mass * (asigj[5] / pPar[j].rho / pPar[j].rho + asigi[5] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][2];

					tempup[1][0] = tempup[1][0] + pPar[j].mass * (asigj[3] / pPar[j].rho / pPar[j].rho + asigi[3] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][0];
					tempup[1][1] = tempup[1][1] + pPar[j].mass * (asigj[1] / pPar[j].rho / pPar[j].rho + asigi[1] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][1];
					tempup[1][2] = tempup[1][2] + pPar[j].mass * (asigj[4] / pPar[j].rho / pPar[j].rho + asigi[4] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][2];

					tempup[2][0] = tempup[2][0] + pPar[j].mass * (asigj[5] / pPar[j].rho / pPar[j].rho + asigi[5] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][0];
					tempup[2][1] = tempup[2][1] + pPar[j].mass * (asigj[4] / pPar[j].rho / pPar[j].rho + asigi[4] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][1];
					tempup[2][2] = tempup[2][2] + pPar[j].mass * (asigj[2] / pPar[j].rho / pPar[j].rho + asigi[2] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][2];
				}
			}

			pPar[i].ax[0] += tempup[0][0] + tempup[0][1] + tempup[0][2];
			pPar[i].ax[1] += tempup[1][0] + tempup[1][1] + tempup[1][2];
			pPar[i].ax[2] += tempup[2][0] + tempup[2][1] + tempup[2][2];
		}
	}
}

/* the minor stress derivations for water particles*/
void clAcce_Fun::fluid_acceleration_minorwater(Particle *pPar, Par_Cell *pParCell, const Para_Fluid *pVFluid,
											   const Para_Pro &pPPro, const Para_GF &pPGf, int cn)
{

	int i, j, k, nc, nj;
	int ntotal;
	double tempup[3][3], gcc[3], porosity;
	double asigi[6], asigj[6];

	ntotal = pPPro.ntotal;

	gcc[0] = pPGf.gx;
	gcc[1] = pPGf.gy;
	gcc[2] = pPGf.gz;

#pragma omp parallel for schedule(static) private(j, k, nc, nj, tempup, asigi, asigj, porosity)
	for (i = 0; i < ntotal; i++)
	{
		if (pPar[i].type == 1)
		{
			porosity = pVFluid[i].porosity;

			for (k = 0; k < 3; k++)
			{
				tempup[k][0] = 0.0;
				tempup[k][1] = 0.0;
				tempup[k][2] = 0.0;
				pPar[i].ax[k] = gcc[k];
			}

			nc = pParCell[i].ninflu;
			for (nj = 1; nj <= nc; nj++)
			{
				j = pParCell[i].influ[nj];

				if (pPar[i].type == pPar[j].type)
				{
					for (k = 0; k < 6; k++)
					{
						asigi[k] = porosity * pPar[i].sig[k];
						asigj[k] = porosity * pPar[j].sig[k];
					}

					tempup[0][0] = tempup[0][0] + pPar[j].mass * (asigj[0] / pPar[j].rho / pPar[j].rho - asigi[0] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][0];
					tempup[0][1] = tempup[0][1] + pPar[j].mass * (asigj[3] / pPar[j].rho / pPar[j].rho - asigi[3] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][1];
					tempup[0][2] = tempup[0][2] + pPar[j].mass * (asigj[5] / pPar[j].rho / pPar[j].rho - asigi[5] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][2];

					tempup[1][0] = tempup[1][0] + pPar[j].mass * (asigj[3] / pPar[j].rho / pPar[j].rho - asigi[3] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][0];
					tempup[1][1] = tempup[1][1] + pPar[j].mass * (asigj[1] / pPar[j].rho / pPar[j].rho - asigi[1] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][1];
					tempup[1][2] = tempup[1][2] + pPar[j].mass * (asigj[4] / pPar[j].rho / pPar[j].rho - asigi[4] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][2];

					tempup[2][0] = tempup[2][0] + pPar[j].mass * (asigj[5] / pPar[j].rho / pPar[j].rho - asigi[5] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][0];
					tempup[2][1] = tempup[2][1] + pPar[j].mass * (asigj[4] / pPar[j].rho / pPar[j].rho - asigi[4] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][1];
					tempup[2][2] = tempup[2][2] + pPar[j].mass * (asigj[2] / pPar[j].rho / pPar[j].rho - asigi[2] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][2];
				}
				else if (pPar[j].type == 3)
				{
					for (k = 0; k < 6; k++)
					{
						asigi[k] = porosity * pPar[i].sig[k];
						asigj[k] = porosity * pPar[j].sig[k];
					}

					tempup[0][0] = tempup[0][0] + pPar[j].mass * (asigj[0] / pPar[j].rho / pPar[j].rho - asigi[0] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][0];
					tempup[0][1] = tempup[0][1] + pPar[j].mass * (asigj[3] / pPar[j].rho / pPar[j].rho - asigi[3] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][1];
					tempup[0][2] = tempup[0][2] + pPar[j].mass * (asigj[5] / pPar[j].rho / pPar[j].rho - asigi[5] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][2];

					tempup[1][0] = tempup[1][0] + pPar[j].mass * (asigj[3] / pPar[j].rho / pPar[j].rho - asigi[3] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][0];
					tempup[1][1] = tempup[1][1] + pPar[j].mass * (asigj[1] / pPar[j].rho / pPar[j].rho - asigi[1] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][1];
					tempup[1][2] = tempup[1][2] + pPar[j].mass * (asigj[4] / pPar[j].rho / pPar[j].rho - asigi[4] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][2];

					tempup[2][0] = tempup[2][0] + pPar[j].mass * (asigj[5] / pPar[j].rho / pPar[j].rho - asigi[5] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][0];
					tempup[2][1] = tempup[2][1] + pPar[j].mass * (asigj[4] / pPar[j].rho / pPar[j].rho - asigi[4] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][1];
					tempup[2][2] = tempup[2][2] + pPar[j].mass * (asigj[2] / pPar[j].rho / pPar[j].rho - asigi[2] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][2];
				}
			}

			pPar[i].ax[0] += tempup[0][0] + tempup[0][1] + tempup[0][2];
			pPar[i].ax[1] += tempup[1][0] + tempup[1][1] + tempup[1][2];
			pPar[i].ax[2] += tempup[2][0] + tempup[2][1] + tempup[2][2];
		}
	}
}

/* the minor stress derivations for air particles*/
void clAcce_Fun::fluid_acceleration_minorair(Particle *pPar, Par_Cell *pParCell, const Para_Fluid *pVFluid,
											 const Para_Pro &pPPro, const Para_GF &pPGf, int cn)
{

	int i, j, k, nc, nj;
	int ntotal;
	double tempup[3][3], gcc[3], porosity;
	double asigi[6], asigj[6];

	ntotal = pPPro.ntotal;

	gcc[0] = pPGf.gx;
	gcc[1] = pPGf.gy;
	gcc[2] = pPGf.gz;

#pragma omp parallel for schedule(static) private(j, k, nc, nj, tempup, asigi, asigj, porosity)
	for (i = 0; i < ntotal; i++)
	{

		if (pPar[i].type == 3)
		{
			porosity = pVFluid[i].porosity;

			for (k = 0; k < 3; k++)
			{
				tempup[k][0] = 0.0;
				tempup[k][1] = 0.0;
				tempup[k][2] = 0.0;
				pPar[i].ax[k] = gcc[k];
			}

			nc = pParCell[i].ninflu;
			for (nj = 1; nj <= nc; nj++)
			{
				j = pParCell[i].influ[nj];

				if (pPar[i].type == pPar[j].type)
				{
					for (k = 0; k < 6; k++)
					{
						asigi[k] = porosity * pPar[i].sig[k];
						asigj[k] = porosity * pPar[j].sig[k];
					}

					tempup[0][0] = tempup[0][0] + pPar[j].mass * (asigj[0] / pPar[j].rho / pPar[j].rho - asigi[0] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][0];
					tempup[0][1] = tempup[0][1] + pPar[j].mass * (asigj[3] / pPar[j].rho / pPar[j].rho - asigi[3] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][1];
					tempup[0][2] = tempup[0][2] + pPar[j].mass * (asigj[5] / pPar[j].rho / pPar[j].rho - asigi[5] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][2];

					tempup[1][0] = tempup[1][0] + pPar[j].mass * (asigj[3] / pPar[j].rho / pPar[j].rho - asigi[3] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][0];
					tempup[1][1] = tempup[1][1] + pPar[j].mass * (asigj[1] / pPar[j].rho / pPar[j].rho - asigi[1] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][1];
					tempup[1][2] = tempup[1][2] + pPar[j].mass * (asigj[4] / pPar[j].rho / pPar[j].rho - asigi[4] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][2];

					tempup[2][0] = tempup[2][0] + pPar[j].mass * (asigj[5] / pPar[j].rho / pPar[j].rho - asigi[5] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][0];
					tempup[2][1] = tempup[2][1] + pPar[j].mass * (asigj[4] / pPar[j].rho / pPar[j].rho - asigi[4] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][1];
					tempup[2][2] = tempup[2][2] + pPar[j].mass * (asigj[2] / pPar[j].rho / pPar[j].rho - asigi[2] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][2];
				}
				else if (pPar[j].type == 1)
				{
					for (k = 0; k < 6; k++)
					{
						asigi[k] = porosity * pPar[i].sig[k];
						asigj[k] = porosity * pPar[j].sig[k];
					}

					tempup[0][0] = tempup[0][0] + pPar[j].mass * (asigj[0] / pPar[j].rho / pPar[j].rho - asigi[0] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][0];
					tempup[0][1] = tempup[0][1] + pPar[j].mass * (asigj[3] / pPar[j].rho / pPar[j].rho - asigi[3] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][1];
					tempup[0][2] = tempup[0][2] + pPar[j].mass * (asigj[5] / pPar[j].rho / pPar[j].rho - asigi[5] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][2];

					tempup[1][0] = tempup[1][0] + pPar[j].mass * (asigj[3] / pPar[j].rho / pPar[j].rho - asigi[3] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][0];
					tempup[1][1] = tempup[1][1] + pPar[j].mass * (asigj[1] / pPar[j].rho / pPar[j].rho - asigi[1] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][1];
					tempup[1][2] = tempup[1][2] + pPar[j].mass * (asigj[4] / pPar[j].rho / pPar[j].rho - asigi[4] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][2];

					tempup[2][0] = tempup[2][0] + pPar[j].mass * (asigj[5] / pPar[j].rho / pPar[j].rho - asigi[5] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][0];
					tempup[2][1] = tempup[2][1] + pPar[j].mass * (asigj[4] / pPar[j].rho / pPar[j].rho - asigi[4] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][1];
					tempup[2][2] = tempup[2][2] + pPar[j].mass * (asigj[2] / pPar[j].rho / pPar[j].rho - asigi[2] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][2];
				}
			}

			pPar[i].ax[0] += tempup[0][0] + tempup[0][1] + tempup[0][2];
			pPar[i].ax[1] += tempup[1][0] + tempup[1][1] + tempup[1][2];
			pPar[i].ax[2] += tempup[2][0] + tempup[2][1] + tempup[2][2];
		}
	}
}

/* the stress derivations for soil particles without regularization*/
void clAcce_Fun::soil_acceleration_noreg(Particle *pPar, Par_Cell *pParCell, const Para_Soil *pParti_ConsPara,
										 const Para_Pro &pPPro, const Para_GF &pPGf, int cn)
{

	int i, j, k, nc, nj;
	int ntotal;
	double tempup[3][3], gcc[3], prei[3], prej[3], porosity;
	double asigi[6], asigj[6], satu;

	ntotal = pPPro.ntotal;

	gcc[0] = pPGf.gx;
	gcc[1] = pPGf.gy;
	gcc[2] = pPGf.gz;

#pragma omp parallel for schedule(static) private(j, k, nc, nj, tempup, asigi, asigj, prei, prej, satu, porosity)

	for (i = 0; i < ntotal; i++)
	{

		if (pPar[i].type == 2)
		{
			porosity = pParti_ConsPara[i].porosity;

			for (k = 0; k < 3; k++)
			{
				tempup[k][0] = 0.0;
				tempup[k][1] = 0.0;
				tempup[k][2] = 0.0;
				pPar[i].ax[k] = gcc[k];
			}

			nc = pParCell[i].ninflu;
			for (nj = 1; nj <= nc; nj++)
			{
				j = pParCell[i].influ[nj];

				if (pPar[i].type == pPar[j].type)
				{
					satu = pPar[j].satu;
					for (k = 0; k < 6; k++)
					{
						asigi[k] = pPar[i].sig[k];
						asigj[k] = pPar[j].sig[k];
					}
					for (k = 0; k < 3; k++)
					{
						prei[k] = (1 - porosity) * (-satu * pPar[i].prew + (1 - satu) * pPar[i].prea);
						prej[k] = (1 - porosity) * (-satu * pPar[j].prew + (1 - satu) * pPar[j].prea);
					}

					tempup[0][0] = tempup[0][0] + pPar[j].mass * (asigj[0] / pPar[j].rho / pPar[j].rho + asigi[0] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][0];
					tempup[0][1] = tempup[0][1] + pPar[j].mass * (asigj[3] / pPar[j].rho / pPar[j].rho + asigi[3] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][1];
					tempup[0][2] = tempup[0][2] + pPar[j].mass * (asigj[5] / pPar[j].rho / pPar[j].rho + asigi[5] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][2];

					tempup[1][0] = tempup[1][0] + pPar[j].mass * (asigj[3] / pPar[j].rho / pPar[j].rho + asigi[3] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][0];
					tempup[1][1] = tempup[1][1] + pPar[j].mass * (asigj[1] / pPar[j].rho / pPar[j].rho + asigi[1] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][1];
					tempup[1][2] = tempup[1][2] + pPar[j].mass * (asigj[4] / pPar[j].rho / pPar[j].rho + asigi[4] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][2];

					tempup[2][0] = tempup[2][0] + pPar[j].mass * (asigj[5] / pPar[j].rho / pPar[j].rho + asigi[5] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][0];
					tempup[2][1] = tempup[2][1] + pPar[j].mass * (asigj[4] / pPar[j].rho / pPar[j].rho + asigi[4] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][1];
					tempup[2][2] = tempup[2][2] + pPar[j].mass * (asigj[2] / pPar[j].rho / pPar[j].rho + asigi[2] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][2];

					tempup[0][0] = tempup[0][0] + pPar[j].mass * (prej[0] / pPar[j].rho / pPar[j].rho + prei[0] / pParti_ConsPara[i].dens / pParti_ConsPara[i].dens) * pParCell[i].wij[nj][0];
					tempup[1][1] = tempup[1][1] + pPar[j].mass * (prej[1] / pPar[j].rho / pPar[j].rho + prei[1] / pParti_ConsPara[i].dens / pParti_ConsPara[i].dens) * pParCell[i].wij[nj][1];
					tempup[2][2] = tempup[2][2] + pPar[j].mass * (prej[2] / pPar[j].rho / pPar[j].rho + prei[2] / pParti_ConsPara[i].dens / pParti_ConsPara[i].dens) * pParCell[i].wij[nj][2];
				}
			}

			pPar[i].ax[0] += tempup[0][0] + tempup[0][1] + tempup[0][2];
			pPar[i].ax[1] += tempup[1][0] + tempup[1][1] + tempup[1][2];
			pPar[i].ax[2] += tempup[2][0] + tempup[2][1] + tempup[2][2];
		}
	}
}

//acceleration of strcture particles
void clAcce_Fun::structure_acceleration_noreg(Particle *pPar, Par_Cell *pParCell,
											const Para_Pro &pPPro, const Para_GF &pPGf, int cn)
{

	int i, j, k, nc, nj;
	int ntotal;
	double tempup[3][3], gcc[3], prei[3], prej[3];
	double asigi[6], asigj[6], satu;

	ntotal = pPPro.ntotal;

	gcc[0] = pPGf.gx;
	gcc[1] = pPGf.gy;
	gcc[2] = pPGf.gz;

#pragma omp parallel for schedule(static) private(j, k, nc, nj, tempup, asigi, asigj, prei, prej, satu)

	for (i = 0; i < ntotal; i++)
	{

		if (pPar[i].type == 4)
		{

			for (k = 0; k < 3; k++)
			{
				tempup[k][0] = 0.0;
				tempup[k][1] = 0.0;
				tempup[k][2] = 0.0;
				pPar[i].ax[k] = gcc[k];
			}

			nc = pParCell[i].ninflu;
			for (nj = 1; nj <= nc; nj++)
			{
				j = pParCell[i].influ[nj];

				if (pPar[i].type == pPar[j].type)
				{
					satu = pPar[j].satu;
					for (k = 0; k < 6; k++)
					{
						asigi[k] = pPar[i].sig[k];
						asigj[k] = pPar[j].sig[k];
					}
					for (k = 0; k < 3; k++)
					{
						prei[k] = -satu * pPar[i].prew + (1 - satu) * pPar[i].prea;
						prej[k] = -satu * pPar[j].prew + (1 - satu) * pPar[j].prea;
					}

					tempup[0][0] = tempup[0][0] + pPar[j].mass * (asigj[0] / pPar[j].rho / pPar[j].rho + asigi[0] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][0];
					tempup[0][1] = tempup[0][1] + pPar[j].mass * (asigj[3] / pPar[j].rho / pPar[j].rho + asigi[3] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][1];
					tempup[0][2] = tempup[0][2] + pPar[j].mass * (asigj[5] / pPar[j].rho / pPar[j].rho + asigi[5] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][2];

					tempup[1][0] = tempup[1][0] + pPar[j].mass * (asigj[3] / pPar[j].rho / pPar[j].rho + asigi[3] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][0];
					tempup[1][1] = tempup[1][1] + pPar[j].mass * (asigj[1] / pPar[j].rho / pPar[j].rho + asigi[1] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][1];
					tempup[1][2] = tempup[1][2] + pPar[j].mass * (asigj[4] / pPar[j].rho / pPar[j].rho + asigi[4] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][2];

					tempup[2][0] = tempup[2][0] + pPar[j].mass * (asigj[5] / pPar[j].rho / pPar[j].rho + asigi[5] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][0];
					tempup[2][1] = tempup[2][1] + pPar[j].mass * (asigj[4] / pPar[j].rho / pPar[j].rho + asigi[4] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][1];
					tempup[2][2] = tempup[2][2] + pPar[j].mass * (asigj[2] / pPar[j].rho / pPar[j].rho + asigi[2] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][2];

					tempup[0][0] = tempup[0][0] + pPar[j].mass * (prej[0] / pPar[j].rho / pPar[j].rho - prei[0] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][0];
					tempup[1][1] = tempup[1][1] + pPar[j].mass * (prej[1] / pPar[j].rho / pPar[j].rho - prei[1] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][1];
					tempup[2][2] = tempup[2][2] + pPar[j].mass * (prej[2] / pPar[j].rho / pPar[j].rho - prei[2] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][2];

				}
			}

			pPar[i].ax[0] += tempup[0][0] + tempup[0][1] + tempup[0][2];
			pPar[i].ax[1] += tempup[1][0] + tempup[1][1] + tempup[1][2];
			pPar[i].ax[2] += tempup[2][0] + tempup[2][1] + tempup[2][2];
		}
	}
}

void clAcce_Fun::coupled_Water_Soil_Structure(Particle *pPar, Par_Cell *pParCell, const Para_Fluid *pVFluid, const Para_Soil *pParti_ConsPara,
											  const Para_Pro &pPPro, const Para_GF &pPGf, int cn)
{

	int i, j, k, nc, nj;
	int ntotal;
	double tempup[3][3];
	double asigi[6], asigj[6];

	ntotal = pPPro.ntotal;

#pragma omp parallel for schedule(static) private(j, k, nc, nj, tempup, asigi, asigj)

	for (i = 0; i < ntotal; i++)
	{

		if (pPar[i].type == 4 || pPar[i].type == 2 || pPar[i].type == 1)
		{

			for (k = 0; k < 3; k++)
			{
				tempup[k][0] = 0.0;
				tempup[k][1] = 0.0;
				tempup[k][2] = 0.0;
			}

			nc = pParCell[i].ninflu;
			for (nj = 1; nj <= nc; nj++)
			{
				j = pParCell[i].influ[nj];

				if (pPar[i].type == 4 && pPar[j].type == 2)
				{
					for (k = 0; k < 6; k++)
					{
						asigi[k] = pPar[i].sig[k];
						asigj[k] = pPar[j].sig[k];
					}

					tempup[0][0] = tempup[0][0] + pPar[j].mass * (asigj[0] / pPar[j].rho / pPar[j].rho + asigi[0] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][0];
					tempup[0][1] = tempup[0][1] + pPar[j].mass * (asigj[3] / pPar[j].rho / pPar[j].rho + asigi[3] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][1];
					tempup[0][2] = tempup[0][2] + pPar[j].mass * (asigj[5] / pPar[j].rho / pPar[j].rho + asigi[5] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][2];

					tempup[1][0] = tempup[1][0] + pPar[j].mass * (asigj[3] / pPar[j].rho / pPar[j].rho + asigi[3] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][0];
					tempup[1][1] = tempup[1][1] + pPar[j].mass * (asigj[1] / pPar[j].rho / pPar[j].rho + asigi[1] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][1];
					tempup[1][2] = tempup[1][2] + pPar[j].mass * (asigj[4] / pPar[j].rho / pPar[j].rho + asigi[4] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][2];

					tempup[2][0] = tempup[2][0] + pPar[j].mass * (asigj[5] / pPar[j].rho / pPar[j].rho + asigi[5] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][0];
					tempup[2][1] = tempup[2][1] + pPar[j].mass * (asigj[4] / pPar[j].rho / pPar[j].rho + asigi[4] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][1];
					tempup[2][2] = tempup[2][2] + pPar[j].mass * (asigj[2] / pPar[j].rho / pPar[j].rho + asigi[2] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][2];;
				}

				if (pPar[i].type == 2 && pPar[j].type == 4)
				{
					for (k = 0; k < 6; k++)
					{
						asigi[k] = pPar[i].sig[k];
						asigj[k] = pPar[j].sig[k];
					}

					tempup[0][0] = tempup[0][0] + pPar[j].mass * (asigj[0] / pPar[j].rho / pPar[j].rho + asigi[0] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][0];
					tempup[0][1] = tempup[0][1] + pPar[j].mass * (asigj[3] / pPar[j].rho / pPar[j].rho + asigi[3] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][1];
					tempup[0][2] = tempup[0][2] + pPar[j].mass * (asigj[5] / pPar[j].rho / pPar[j].rho + asigi[5] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][2];

					tempup[1][0] = tempup[1][0] + pPar[j].mass * (asigj[3] / pPar[j].rho / pPar[j].rho + asigi[3] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][0];
					tempup[1][1] = tempup[1][1] + pPar[j].mass * (asigj[1] / pPar[j].rho / pPar[j].rho + asigi[1] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][1];
					tempup[1][2] = tempup[1][2] + pPar[j].mass * (asigj[4] / pPar[j].rho / pPar[j].rho + asigi[4] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][2];

					tempup[2][0] = tempup[2][0] + pPar[j].mass * (asigj[5] / pPar[j].rho / pPar[j].rho + asigi[5] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][0];
					tempup[2][1] = tempup[2][1] + pPar[j].mass * (asigj[4] / pPar[j].rho / pPar[j].rho + asigi[4] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][1];
					tempup[2][2] = tempup[2][2] + pPar[j].mass * (asigj[2] / pPar[j].rho / pPar[j].rho + asigi[2] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][2];
				}

				if (pPar[i].type == 4 && pPar[j].type == 1)
				{
					for (k = 0; k < 6; k++)
					{
						asigi[k] = pPar[i].sig[k];
						asigj[k] = pPar[j].sig[k];
					}

					tempup[0][0] = tempup[0][0] + pPar[j].mass * (asigj[0] / pPar[j].rho / pPar[j].rho + asigi[0] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][0];
					tempup[0][1] = tempup[0][1] + pPar[j].mass * (asigj[3] / pPar[j].rho / pPar[j].rho + asigi[3] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][1];
					tempup[0][2] = tempup[0][2] + pPar[j].mass * (asigj[5] / pPar[j].rho / pPar[j].rho + asigi[5] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][2];

					tempup[1][0] = tempup[1][0] + pPar[j].mass * (asigj[3] / pPar[j].rho / pPar[j].rho + asigi[3] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][0];
					tempup[1][1] = tempup[1][1] + pPar[j].mass * (asigj[1] / pPar[j].rho / pPar[j].rho + asigi[1] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][1];
					tempup[1][2] = tempup[1][2] + pPar[j].mass * (asigj[4] / pPar[j].rho / pPar[j].rho + asigi[4] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][2];

					tempup[2][0] = tempup[2][0] + pPar[j].mass * (asigj[5] / pPar[j].rho / pPar[j].rho + asigi[5] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][0];
					tempup[2][1] = tempup[2][1] + pPar[j].mass * (asigj[4] / pPar[j].rho / pPar[j].rho + asigi[4] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][1];
					tempup[2][2] = tempup[2][2] + pPar[j].mass * (asigj[2] / pPar[j].rho / pPar[j].rho + asigi[2] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][2];
				}

				if (pPar[i].type == 1 && pPar[j].type == 4)
				{
					for (k = 0; k < 6; k++)
					{
						asigi[k] = pPar[i].sig[k];
						asigj[k] = pPar[j].sig[k];
					}

					tempup[0][0] = tempup[0][0] + pPar[j].mass * (asigj[0] / pPar[j].rho / pPar[j].rho - asigi[0] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][0];
					tempup[0][1] = tempup[0][1] + pPar[j].mass * (asigj[3] / pPar[j].rho / pPar[j].rho - asigi[3] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][1];
					tempup[0][2] = tempup[0][2] + pPar[j].mass * (asigj[5] / pPar[j].rho / pPar[j].rho - asigi[5] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][2];

					tempup[1][0] = tempup[1][0] + pPar[j].mass * (asigj[3] / pPar[j].rho / pPar[j].rho - asigi[3] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][0];
					tempup[1][1] = tempup[1][1] + pPar[j].mass * (asigj[1] / pPar[j].rho / pPar[j].rho - asigi[1] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][1];
					tempup[1][2] = tempup[1][2] + pPar[j].mass * (asigj[4] / pPar[j].rho / pPar[j].rho - asigi[4] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][2];

					tempup[2][0] = tempup[2][0] + pPar[j].mass * (asigj[5] / pPar[j].rho / pPar[j].rho - asigi[5] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][0];
					tempup[2][1] = tempup[2][1] + pPar[j].mass * (asigj[4] / pPar[j].rho / pPar[j].rho - asigi[4] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][1];
					tempup[2][2] = tempup[2][2] + pPar[j].mass * (asigj[2] / pPar[j].rho / pPar[j].rho - asigi[2] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][2];
				}
			}

			pPar[i].ax[0] += tempup[0][0] + tempup[0][1] + tempup[0][2];
			pPar[i].ax[1] += tempup[1][0] + tempup[1][1] + tempup[1][2];
			pPar[i].ax[2] += tempup[2][0] + tempup[2][1] + tempup[2][2];
		}
	}
}

//calculating stress effect of boudary particles on moving particles
void clAcce_Fun::boundary_moving_stress_effect(Particle *pPar, Par_Cell *pParCell, const Para_Pro &pPPro, int cn)
{

	int i, j, k, nc, nj;
	int ntotal;
	double tempup[3][3];
	double asigi[6], asigj[6], satu;

	ntotal = pPPro.ntotal;

#pragma omp parallel for schedule(static) private(j, k, nc, nj, tempup, asigi, asigj, satu)
	for (i = 0; i < ntotal; i++)
	{

		if (pPar[i].type != 0 && pPar[i].type != 7)
		{

			for (k = 0; k < 3; k++)
			{
				tempup[k][0] = 0.0;
				tempup[k][1] = 0.0;
				tempup[k][2] = 0.0;
			}

			nc = pParCell[i].ninflu;
			for (nj = 1; nj <= nc; nj++)
			{
				j = pParCell[i].influ[nj];
				if (pPar[j].type == 0)
				{
					for (k = 0; k < 6; k++)
					{
						asigi[k] = 0.0;
						asigj[k] = 0.0;
					}

					//values of asigi and asigj
					if (pPar[i].type == 1 || pPar[i].type == 3)
					{
						asigi[0] = pPar[i].sig[0];
						asigj[0] = asigi[0];
						asigi[1] = pPar[i].sig[1];
						asigj[1] = asigi[1];
						asigi[2] = pPar[i].sig[2];
						asigj[2] = asigi[2];

						asigi[3] = pPar[i].sig[3];
						asigj[3] = asigi[3];
						asigi[4] = pPar[i].sig[4];
						asigj[4] = asigi[4];
						asigi[5] = pPar[i].sig[5];
						asigj[5] = asigi[5];
					}
					else if (pPar[i].type == 2 || pPar[i].type == 4)
					{
						satu = pPar[i].satu;
						asigi[0] = pPar[i].sig[0] - satu * pPar[i].prew + (1 - satu) * pPar[i].prea;
						asigj[0] = asigi[0];
						asigi[1] = pPar[i].sig[1] - satu * pPar[i].prew + (1 - satu) * pPar[i].prea;
						asigj[1] = asigi[1];
						asigi[2] = pPar[i].sig[2] - satu * pPar[i].prew + (1 - satu) * pPar[i].prea;
						asigj[2] = asigi[2];

						asigi[3] = pPar[i].sig[3];
						asigj[3] = asigi[3];
						asigi[4] = pPar[i].sig[4];
						asigj[4] = asigi[4];
						asigi[5] = pPar[i].sig[5];
						asigj[5] = asigi[5];
					}
					//acceleration calculation
					tempup[0][0] = tempup[0][0] + pPar[j].mass * (asigj[0] / pPar[j].rho / pPar[j].rho + asigi[0] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][0];
					tempup[0][1] = tempup[0][1] + pPar[j].mass * (asigj[3] / pPar[j].rho / pPar[j].rho + asigi[3] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][1];
					tempup[0][2] = tempup[0][2] + pPar[j].mass * (asigj[5] / pPar[j].rho / pPar[j].rho + asigi[5] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][2];

					tempup[1][0] = tempup[1][0] + pPar[j].mass * (asigj[3] / pPar[j].rho / pPar[j].rho + asigi[3] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][0];
					tempup[1][1] = tempup[1][1] + pPar[j].mass * (asigj[1] / pPar[j].rho / pPar[j].rho + asigi[1] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][1];
					tempup[1][2] = tempup[1][2] + pPar[j].mass * (asigj[4] / pPar[j].rho / pPar[j].rho + asigi[4] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][2];

					tempup[2][0] = tempup[2][0] + pPar[j].mass * (asigj[5] / pPar[j].rho / pPar[j].rho + asigi[5] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][0];
					tempup[2][1] = tempup[2][1] + pPar[j].mass * (asigj[4] / pPar[j].rho / pPar[j].rho + asigi[4] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][1];
					tempup[2][2] = tempup[2][2] + pPar[j].mass * (asigj[2] / pPar[j].rho / pPar[j].rho + asigi[2] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][2];
				}
			}

			pPar[i].ax[0] += tempup[0][0] + tempup[0][1] + tempup[0][2];
			pPar[i].ax[1] += tempup[1][0] + tempup[1][1] + tempup[1][2];
			pPar[i].ax[2] += tempup[2][0] + tempup[2][1] + tempup[2][2];
		}
	}
}


//calculating stress effect of boudary particles on moving particles
void clAcce_Fun::boundary_stress_effect_tran(Particle* pPar, Par_Cell* pParCell, clVar_Boundary* pParti_VariBndy, 
	const Para_Pro& pPPro, int cn)
{

	int i, j, k, nc, nj;
	int ntotal;
	double tempup[3][3];
	double asigi[6], asigj[6], satu;

	ntotal = pPPro.ntotal;

#pragma omp parallel for schedule(static) private(j, k, nc, nj, tempup, asigi, asigj, satu)
	for (i = 0; i < ntotal; i++)
	{

		if (pPar[i].type != 0 && pPar[i].type != 7)
		{

			for (k = 0; k < 3; k++)
			{
				tempup[k][0] = 0.0;
				tempup[k][1] = 0.0;
				tempup[k][2] = 0.0;
			}

			nc = pParCell[i].ninflu;
			for (nj = 1; nj <= nc; nj++)
			{
				j = pParCell[i].influ[nj];
				if (pPar[j].type == 0)
				{
					for (k = 0; k < 6; k++)
					{
						asigi[k] = 0.0;
						asigj[k] = 0.0;
					}

					if (pPar[i].type == 1) {
						asigi[0] = pPar[i].sig[0];
						asigi[1] = pPar[i].sig[1];
						asigi[2] = pPar[i].sig[2];
						asigi[3] = pPar[i].sig[3];
						asigi[4] = pPar[i].sig[4];
						asigi[5] = pPar[i].sig[5];

						asigj[0] = pParti_VariBndy[j].sig_water[0];
						asigj[1] = pParti_VariBndy[j].sig_water[1];
						asigj[2] = pParti_VariBndy[j].sig_water[2];
						asigj[3] = pParti_VariBndy[j].sig_water[3];
						asigj[4] = pParti_VariBndy[j].sig_water[4];
						asigj[5] = pParti_VariBndy[j].sig_water[5];
					}
					else if (pPar[i].type == 2) {
						satu = pPar[i].satu;
						asigi[0] = pPar[i].sig[0] - satu * pPar[i].prew + (1 - satu) * pPar[i].prea;
						asigi[1] = pPar[i].sig[1] - satu * pPar[i].prew + (1 - satu) * pPar[i].prea;
						asigi[2] = pPar[i].sig[2] - satu * pPar[i].prew + (1 - satu) * pPar[i].prea;
						asigi[3] = pPar[i].sig[3];
						asigi[4] = pPar[i].sig[4];
						asigi[5] = pPar[i].sig[5];

						asigj[0] = pParti_VariBndy[j].sig_soil[0];
						asigj[1] = pParti_VariBndy[j].sig_soil[1];
						asigj[2] = pParti_VariBndy[j].sig_soil[2];
						asigj[3] = pParti_VariBndy[j].sig_soil[3];
						asigj[4] = pParti_VariBndy[j].sig_soil[4];
						asigj[5] = pParti_VariBndy[j].sig_soil[5];
					}
					else if (pPar[i].type == 3) {
						asigi[0] = pPar[i].sig[0];
						asigi[1] = pPar[i].sig[1];
						asigi[2] = pPar[i].sig[2];
						asigi[3] = pPar[i].sig[3];
						asigi[4] = pPar[i].sig[4];
						asigi[5] = pPar[i].sig[5];

						asigj[0] = pParti_VariBndy[j].sig_air[0];
						asigj[1] = pParti_VariBndy[j].sig_air[1];
						asigj[2] = pParti_VariBndy[j].sig_air[2];
						asigj[3] = pParti_VariBndy[j].sig_air[3];
						asigj[4] = pParti_VariBndy[j].sig_air[4];
						asigj[5] = pParti_VariBndy[j].sig_air[5];
					}
					else if (pPar[i].type == 4) {
						satu = 1.0;

						asigi[0] = pPar[i].sig[0];
						asigi[1] = pPar[i].sig[1];
						asigi[2] = pPar[i].sig[2];
						asigi[3] = pPar[i].sig[3];
						asigi[4] = pPar[i].sig[4];
						asigi[5] = pPar[i].sig[5];

						asigj[0] = pParti_VariBndy[j].sig_struct[0];
						asigj[1] = pParti_VariBndy[j].sig_struct[1];
						asigj[2] = pParti_VariBndy[j].sig_struct[2];
						asigj[3] = pParti_VariBndy[j].sig_struct[3];
						asigj[4] = pParti_VariBndy[j].sig_struct[4];
						asigj[5] = pParti_VariBndy[j].sig_struct[5];
					}

					
					//acceleration calculation
					tempup[0][0] = tempup[0][0] + pPar[j].mass * (asigj[0] / pPar[j].rho / pPar[j].rho + asigi[0] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][0];
					tempup[0][1] = tempup[0][1] + pPar[j].mass * (asigj[3] / pPar[j].rho / pPar[j].rho + asigi[3] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][1];
					tempup[0][2] = tempup[0][2] + pPar[j].mass * (asigj[5] / pPar[j].rho / pPar[j].rho + asigi[5] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][2];

					tempup[1][0] = tempup[1][0] + pPar[j].mass * (asigj[3] / pPar[j].rho / pPar[j].rho + asigi[3] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][0];
					tempup[1][1] = tempup[1][1] + pPar[j].mass * (asigj[1] / pPar[j].rho / pPar[j].rho + asigi[1] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][1];
					tempup[1][2] = tempup[1][2] + pPar[j].mass * (asigj[4] / pPar[j].rho / pPar[j].rho + asigi[4] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][2];

					tempup[2][0] = tempup[2][0] + pPar[j].mass * (asigj[5] / pPar[j].rho / pPar[j].rho + asigi[5] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][0];
					tempup[2][1] = tempup[2][1] + pPar[j].mass * (asigj[4] / pPar[j].rho / pPar[j].rho + asigi[4] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][1];
					tempup[2][2] = tempup[2][2] + pPar[j].mass * (asigj[2] / pPar[j].rho / pPar[j].rho + asigi[2] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][2];
				}
			}

			pPar[i].ax[0] += tempup[0][0] + tempup[0][1] + tempup[0][2];
			pPar[i].ax[1] += tempup[1][0] + tempup[1][1] + tempup[1][2];
			pPar[i].ax[2] += tempup[2][0] + tempup[2][1] + tempup[2][2];
		}
	}
}

//damping effect for the initializing of geostress
void clAcce_Fun::soil_damping_geostress(Particle* pPar, Par_Cell* pParCell, const Para_Pro& pPPro, int cn)
{
	int i, j, k, nc, nj;
	int ntotal, ndim;
	double cd, vx1[3], vx2[3], dt, dr, damp[3];

	ntotal = pPPro.ntotal;
	ndim = pPPro.ndim;
	dr= pPPro.dr;
	dt= pPPro.dt;

#pragma omp parallel for schedule(static) private(j, k, nc, nj, cd, vx1, vx2, damp)
	for (i = 0; i < ntotal; i++)
	{
		if (pPar[i].type == 2 || pPar[i].type == 4) {
			//initialization
			cd = pPar[i].cd;
			for (k = 0; k < ndim; k++) damp[k] = 0.0;

			//damping effect
			nc = pParCell[i].ninflu;
			for (nj = 1; nj <= nc; nj++)
			{
				j = pParCell[i].influ[nj];

				if (pPar[i].type == pPar[j].type)
				{
					for (k = 0; k < ndim; k++)
						damp[k] = damp[k]
						+ pPar[j].mass / pPar[j].rho * cd * (pPar[i].vxp[k] - pPar[j].vxp[k]) * pParCell[i].wij[nj][3];
				}
			}

			//refecting to acceleration
			for (k = 0; k < ndim; k++)
			{
				vx1[k] = pPar[i].vxp[k] + 0.5 * (pPar[i].axp[k] + pPar[i].ax[k]) * dt;
				vx2[k] = pPar[i].vxp[k] + 0.5 * (pPar[i].axp[k] + pPar[i].ax[k] - damp[k]) * dt;
				if (vx1[k] * vx2[k] <= 0.0)
				{
					pPar[i].ax[k] = -2.0 * pPar[i].vxp[k] / dt - pPar[i].axp[k];
				}
				else
					pPar[i].ax[k] -= damp[k];
			}
		}
	}
}

//calculating stress effect of boudary particles on moving particles
void clAcce_Fun::boundary_effect_free_struct(Particle* pPar, Par_Cell* pParCell, const Para_Pro& pPPro, int cn)
{

	int i, j, k, nc, nj;
	int ntotal;
	double tempup[3][3];
	double asigi[6], asigj[6], satu;

	ntotal = pPPro.ntotal;

#pragma omp parallel for schedule(static) private(j, k, nc, nj, tempup, asigi, asigj, satu)
	for (i = 0; i < ntotal; i++)
	{

		if (pPar[i].type != 0 && pPar[i].type != 7)
		{

			for (k = 0; k < 3; k++)
			{
				tempup[k][0] = 0.0;
				tempup[k][1] = 0.0;
				tempup[k][2] = 0.0;
			}

			nc = pParCell[i].ninflu;
			for (nj = 1; nj <= nc; nj++)
			{
				j = pParCell[i].influ[nj];
				if (pPar[j].type == 0)
				{
					for (k = 0; k < 6; k++)
					{
						asigi[k] = 0.0;
						asigj[k] = 0.0;
					}

					//values of asigi and asigj
					if (pPar[i].type == 1 || pPar[i].type == 3)
					{
						asigi[0] = pPar[i].sig[0];
						asigj[0] = asigi[0];
						asigi[1] = pPar[i].sig[1];
						asigj[1] = asigi[1];
						asigi[2] = pPar[i].sig[2];
						asigj[2] = asigi[2];

						asigi[3] = pPar[i].sig[3];
						asigj[3] = -asigi[3];
						asigi[4] = pPar[i].sig[4];
						asigj[4] = -asigi[4];
						asigi[5] = pPar[i].sig[5];
						asigj[5] = -asigi[5];
					}
					else if (pPar[i].type == 2 || pPar[i].type == 4)
					{
						satu = pPar[i].satu;
						asigi[0] = pPar[i].sig[0] - satu * pPar[i].prew + (1 - satu) * pPar[i].prea;
						asigj[0] = asigi[0];
						asigi[1] = pPar[i].sig[1] - satu * pPar[i].prew + (1 - satu) * pPar[i].prea;
						asigj[1] = asigi[1];
						asigi[2] = pPar[i].sig[2] - satu * pPar[i].prew + (1 - satu) * pPar[i].prea;
						asigj[2] = asigi[2];

						asigi[3] = pPar[i].sig[3];
						asigj[3] = asigi[3];
						asigi[4] = pPar[i].sig[4];
						asigj[4] = asigi[4];
						asigi[5] = pPar[i].sig[5];
						asigj[5] = asigi[5];
					}
					//acceleration calculation
					tempup[0][0] = tempup[0][0] + pPar[j].mass * (asigj[0] / pPar[j].rho / pPar[j].rho + asigi[0] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][0];
					tempup[0][1] = tempup[0][1] + pPar[j].mass * (asigj[3] / pPar[j].rho / pPar[j].rho + asigi[3] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][1];
					tempup[0][2] = tempup[0][2] + pPar[j].mass * (asigj[5] / pPar[j].rho / pPar[j].rho + asigi[5] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][2];

					tempup[1][0] = tempup[1][0] + pPar[j].mass * (asigj[3] / pPar[j].rho / pPar[j].rho + asigi[3] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][0];
					tempup[1][1] = tempup[1][1] + pPar[j].mass * (asigj[1] / pPar[j].rho / pPar[j].rho + asigi[1] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][1];
					tempup[1][2] = tempup[1][2] + pPar[j].mass * (asigj[4] / pPar[j].rho / pPar[j].rho + asigi[4] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][2];

					tempup[2][0] = tempup[2][0] + pPar[j].mass * (asigj[5] / pPar[j].rho / pPar[j].rho + asigi[5] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][0];
					tempup[2][1] = tempup[2][1] + pPar[j].mass * (asigj[4] / pPar[j].rho / pPar[j].rho + asigi[4] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][1];
					tempup[2][2] = tempup[2][2] + pPar[j].mass * (asigj[2] / pPar[j].rho / pPar[j].rho + asigi[2] / pPar[i].rho / pPar[i].rho) * pParCell[i].wij[nj][2];
				}
			}

			pPar[i].ax[0] += tempup[0][0] + tempup[0][1] + tempup[0][2];
			pPar[i].ax[1] += tempup[1][0] + tempup[1][1] + tempup[1][2];
			pPar[i].ax[2] += tempup[2][0] + tempup[2][1] + tempup[2][2];
		}
	}
}

