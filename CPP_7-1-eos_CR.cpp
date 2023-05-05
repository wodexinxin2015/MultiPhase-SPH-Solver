/************************************************************************************
Multi-Phase and Parallelized SPH program
--parallelized and optimized for SMP (Symmetric multiprocessing ) architecture.
--Copyright (c) 2018-2023, Weijie ZHANG, GeoHohai, Hohai University.
************************************************************************************/

#include <cmath>
#include "Class_Functions.h"
#include "Header_Option.h"

/* calculation of the pressure according the equation of state*/
void clStraStre_Fun::eos(Particle *pPar, Par_Cell *pParCell, Para_Fluid *pVFluid, const Para_Pro &pPPro, int cn)
{
	int i, j, nc, nj;
	int ntotal;
	double tempup1, tempup2, tempdown1, tempdown2, tempup3, tempdown3, tempup4, tempdown4;
	double tempupstre[6], tempdownstre;

	ntotal = pPPro.ntotal;

	/* calculation of the pressure for water and air*/
#pragma omp parallel for schedule(static)
	for (i = 0; i < ntotal; i++)
	{
		/*water*/
		if (pPar[i].type == 1)
		{
			pPar[i].pre = pPar[i].prep * (pow(pPar[i].rho / pVFluid[i].dens / pVFluid[i].porosity, 7.0) - 1.0);
			pPar[i].pre = pPar[i].pre + pPar[i].prep;
		}
		/*air*/
		if (pPar[i].type == 3)
		{
			pPar[i].pre = 1.403 * pPar[i].prep * (pPar[i].rho / pVFluid[i].dens / pVFluid[i].porosity, -1.0);
			pPar[i].pre = pPar[i].pre + pPar[i].prep;
		}
		/*solid boundary*/
		if (pPar[i].type == 0 && pPar[i].matype == 1)
		{
			pPar[i].pre = pPar[i].prep * (pow(pPar[i].rho / pPar[i].rhop, 7.0) - 1.0);
			pPar[i].pre = pPar[i].pre + pPar[i].prep;
		}
	}

	/* calculation of the water and air pressure for soil, water pressure for boundary*/
#pragma omp parallel for schedule(static) private(j, nc, nj, tempup1, \
																tempup2, tempdown1, tempdown2, tempup3, tempdown3, tempup4, tempdown4, tempupstre, tempdownstre)
	for (i = 0; i < ntotal; i++)
	{

		tempup1 = 0.0;
		tempup2 = 0.0;
		tempup3 = 0.0;
		tempup4 = 0.0;
		tempdown1 = 0.0;
		tempdown2 = 0.0;
		tempdown3 = 0.0;
		tempdown4 = 0.0;
		tempdownstre = 0.0;

		for (j = 0; j < 6; j++)
		{
			tempupstre[j] = 0.0;
		}

		nc = pParCell[i].ninflu;
		for (nj = 1; nj <= nc; nj++)
		{
			j = pParCell[i].influ[nj];
			if (pPar[i].type != pPar[j].type)
			{
				if (pPar[i].type == 2 && pPar[j].type == 1)
				{
					tempup1 = tempup1 + pPar[j].mass * pPar[j].pre * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempdown1 = tempdown1 + pPar[j].mass * pParCell[i].wij[nj][3] / pPar[j].rho;
				}

				if (pPar[i].type == 2 && pPar[j].type == 3)
				{
					tempup2 = tempup2 + pPar[j].mass * pPar[j].pre * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempdown2 = tempdown2 + pPar[j].mass * pParCell[i].wij[nj][3] / pPar[j].rho;
				}

				if (pPar[i].type == 4 && pPar[j].type == 1)
				{
					tempup3 = tempup3 + pPar[j].mass * pPar[j].pre * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempdown3 = tempdown3 + pPar[j].mass * pParCell[i].wij[nj][3] / pPar[j].rho;
				}

				if (pPar[i].type == 4 && pPar[j].type == 3)
				{
					tempup4 = tempup4 + pPar[j].mass * pPar[j].pre * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempdown4 = tempdown4 + pPar[j].mass * pParCell[i].wij[nj][3] / pPar[j].rho;
				}

				if (pPar[i].type == 0 && pPar[j].type == 1)
				{
					tempup1 = tempup1 + pPar[j].mass * pPar[j].pre * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempdown1 = tempdown1 + pPar[j].mass * pParCell[i].wij[nj][3] / pPar[j].rho;
				}

				if (pPar[i].type == 0 && pPar[j].type == 3)
				{
					tempup2 = tempup2 + pPar[j].mass * pPar[j].pre * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempdown2 = tempdown2 + pPar[j].mass * pParCell[i].wij[nj][3] / pPar[j].rho;
				}

				if (pPar[i].type == 0 && pPar[j].type == 2)
				{
					tempupstre[0] = tempupstre[0] + pPar[j].mass * pPar[j].sig[0] * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempupstre[1] = tempupstre[1] + pPar[j].mass * pPar[j].sig[1] * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempupstre[2] = tempupstre[2] + pPar[j].mass * pPar[j].sig[2] * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempupstre[3] = tempupstre[3] + pPar[j].mass * pPar[j].sig[3] * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempupstre[4] = tempupstre[4] + pPar[j].mass * pPar[j].sig[4] * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempupstre[5] = tempupstre[5] + pPar[j].mass * pPar[j].sig[5] * pParCell[i].wij[nj][3] / pPar[j].rho;

					tempdownstre = tempdownstre + pPar[j].mass * pParCell[i].wij[nj][3] / pPar[j].rho;
				}
			}
		}

		if (pPar[i].type == 2)
		{
			if (tempdown1 != 0.0)
				pPar[i].prew = tempup1 / tempdown1;
			else
				pPar[i].prew = tempup1;

			if (tempdown2 != 0.0)
				pPar[i].prea = tempup2 / tempdown2;
			else
				pPar[i].prea = tempup2;
		}
		else if (pPar[i].type == 0)
		{
			if (tempdown1 != 0.0)
				pPar[i].prew = tempup1 / tempdown1;
			else
				pPar[i].prew = tempup1;

			if (tempdown2 != 0.0)
				pPar[i].prea = tempup2 / tempdown2;
			else
				pPar[i].prea = tempup2;

			if (tempdownstre != 0.0)
			{
				pPar[i].sig[0] = tempupstre[0] / tempdownstre;
				pPar[i].sig[1] = tempupstre[1] / tempdownstre;
				pPar[i].sig[2] = tempupstre[2] / tempdownstre;
				pPar[i].sig[3] = tempupstre[3] / tempdownstre;
				pPar[i].sig[4] = tempupstre[4] / tempdownstre;
				pPar[i].sig[5] = tempupstre[5] / tempdownstre;
			}
			else
			{
				pPar[i].sig[0] = tempupstre[0];
				pPar[i].sig[1] = tempupstre[1];
				pPar[i].sig[2] = tempupstre[2];
				pPar[i].sig[3] = tempupstre[3];
				pPar[i].sig[4] = tempupstre[4];
				pPar[i].sig[5] = tempupstre[5];
			}
		}
		else if (pPar[i].type == 4)
		{
			if (tempdown3 != 0.0)
				pPar[i].prew = tempup3 / tempdown3;
			else
				pPar[i].prew = tempup3;

			if (tempdown4 != 0.0)
				pPar[i].prea = tempup4 / tempdown4;
			else
				pPar[i].prea = tempup4;
		}
	}
}
