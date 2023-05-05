/************************************************************************************
Multi-Phase and Parallelized SPH program
--parallelized and optimized for SMP (Symmetric multiprocessing ) architecture.
--Copyright (c) 2018-2023, Weijie ZHANG, GeoHohai, Hohai University.
************************************************************************************/

#include "Header_Option.h"
#include "Class_Functions.h"

/* XSPH Correction */
void clVel_Fun::veocity_update(Particle *pPar, const Para_Pro &pPPro, const Para_Boundary &pPBou, int flv)
{
	int i, k;
	double ath;

	double dt = pPPro.dt;
	int ntotal = pPPro.ntotal;
	int ndim = pPPro.ndim;
	int inip, l;

	if (pPBou.if_move == 0)
	{ //static boundary and seismic virtual boundary
#pragma omp parallel for schedule(static) private(k, ath)
		for (i = 0; i < ntotal; i++)
		{
			if (pPar[i].type == 0 && pPar[i].matype < 2)
			{
				for (k = 0; k < ndim; k++)
				{
					pPar[i].interf[k] = 0.0;
					ath = (pPar[i].interf[k] + pPar[i].interf[k]) * 0.5;
					pPar[i].interfss[k] = pPar[i].interfss[k] + ath * dt;
				}
			}
			else if (pPar[i].type == 0 && pPar[i].matype == 3)
			{
				for (k = 0; k < ndim; k++)
				{
					ath = (pPar[i].interf[k] + pPar[i].interf[k]) * 0.5;
					pPar[i].interfss[k] = pPar[i].interfss[k] + ath * dt;
				}
			}
			else if (pPar[i].type == 0 && pPar[i].matype == 5 && flv == 1)
			{
				for (k = 0; k < ndim; k++)
				{
					ath = (pPar[i].interf[k] + pPar[i].interf[k]) * 0.5;
					pPar[i].interfss[k] = pPar[i].interfss[k] + ath * dt;
				}
			}
			else if (pPar[i].type != 0)
			{
				for (k = 0; k < ndim; k++)
				{
					ath = (pPar[i].axp[k] + pPar[i].ax[k]) * 0.5;
					pPar[i].vx[k] = pPar[i].vxp[k] + ath * dt;
				}
			}
		}
	}
	else if (pPBou.if_move == 1)
	{ //moving boundary by velocity
		l = pPPro.l;
		inip = pPPro.inip;
#pragma omp parallel for schedule(static) private(k, ath)
		for (i = 0; i < ntotal; i++)
		{
			if ((pPar[i].type == 0 && pPar[i].matype == 2) && l > inip)
			{
				for (k = 0; k < ndim; k++)
				{
					pPar[i].interfss[k] = pPBou.move_vx[k];
				}
			}
			else if (pPar[i].type == 0)
			{
				for (k = 0; k < ndim; k++)
				{
					pPar[i].interfss[k] = 0.0;
				}
			}
			else if (pPar[i].type != 0)
			{
				for (k = 0; k < ndim; k++)
				{
					ath = (pPar[i].axp[k] + pPar[i].ax[k]) * 0.5;
					pPar[i].vx[k] = pPar[i].vxp[k] + ath * dt;
				}
			}
		}
	}
	else if (pPBou.if_move == 2)
	{ //moving soil by velocity
		l = pPPro.l;
		inip = pPPro.inip;
#pragma omp parallel for schedule(static) private(k, ath)
		for (i = 0; i < ntotal; i++)
		{

			if ((pPar[i].type == 2 && pPar[i].etype == 5) && l > inip)
			{
				for (k = 0; k < ndim; k++)
				{
					pPar[i].axp[k] = 0.0;
					pPar[i].vxp[k] = pPBou.move_vx[k];
					pPar[i].vx[k] = pPBou.move_vx[k];
					pPar[i].ax[k] = 0.0;
				}
			}
			for (k = 0; k < ndim; k++)
			{
				ath = (pPar[i].axp[k] + pPar[i].ax[k])  * 0.5;
				pPar[i].vx[k] = pPar[i].vxp[k] + ath * dt;
			}
		}
	}
	else if (pPBou.if_move == 3)
	{ //moving structure by velocity
		l = pPPro.l;
		inip = pPPro.inip;
#pragma omp parallel for schedule(static) private(k, ath)
		for (i = 0; i < ntotal; i++)
		{

			if ((pPar[i].type == 4 && pPar[i].etype == 5) && l > inip)
			{
				for (k = 0; k < ndim; k++)
				{
					pPar[i].axp[k] = 0.0;
					pPar[i].vxp[k] = pPBou.move_vx[k];
					pPar[i].vx[k] = pPBou.move_vx[k];
					pPar[i].ax[k] = 0.0;
				}
			}
			for (k = 0; k < ndim; k++)
			{
				ath = (pPar[i].axp[k] + pPar[i].ax[k])  * 0.5;
				pPar[i].vx[k] = pPar[i].vxp[k] + ath * dt;
			}
		}
	}
}

//xsph correction
void clVel_Fun::xsph(Particle* pPar, Par_Cell* pParCell, cl_bndy_pair* pBndy_pair, clVar_Boundary* pParti_VariBndy,
	const Para_Pro& pPPro, double(*vx_temp)[3], int bndytype)
{
	int i, j, k, nc, nj;
	double xishu, tempup1[3], tempdown;
	double vxj[3];

	int ntotal = pPPro.ntotal;
	int ndim = pPPro.ndim;

#pragma omp parallel for schedule(static) private(j, k, nc, nj, xishu, tempup1, tempdown, vxj)
	for (i = 0; i < ntotal; i++)
	{
		tempdown = 0.0;
		xishu = 0.0;
		for (k = 0; k < ndim; k++)
		{
			tempup1[k] = 0.0;
			vx_temp[i][k] = 0.0;
		}

		nc = pParCell[i].ninflu;
		for (nj = 1; nj <= nc; nj++)
		{
			j = pParCell[i].influ[nj];
			if (pPar[i].type == pPar[j].type && (pPar[i].type != 0 && pPar[i].type != 7))
			{

				if (pPar[i].type == 1)
					xishu = pPPro.xsph[0];
				else if (pPar[i].type == 2)
					xishu = pPPro.xsph[1];
				else if (pPar[i].type == 3)
					xishu = pPPro.xsph[2];
				else if (pPar[i].type == 4)
					xishu = pPPro.xsph[1];

				for (k = 0; k < ndim; k++)
					tempup1[k] = tempup1[k] + xishu * pPar[j].mass * (pPar[j].vx[k] - pPar[i].vx[k]) * pParCell[i].wij[nj][3] / pPar[j].rho;

				tempdown = tempdown + pPar[j].mass * pParCell[i].wij[nj][3] / pPar[j].rho;
			}
			else if ((pPar[i].type == 1 || pPar[i].type == 2 || pPar[i].type == 3 || pPar[i].type == 4) && pPar[j].type == 0) {
				//settting virtual velocity
				vxj[0] = 0.0;
				vxj[1] = 0.0;
				vxj[2] = 0.0;
				if (bndytype == 1) {
					vxj[0] = pPar[i].vx[0];
					vxj[1] = pPar[i].vx[1];
					vxj[2] = pPar[i].vx[2];
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

				if (pPar[i].type == 1)
					xishu = 0.006;
				else if (pPar[i].type == 2)
					xishu = pPPro.xsph[1] / 5.0;
				else if (pPar[i].type == 3)
					xishu = 0.006;
				else if (pPar[i].type == 4)
					xishu = pPPro.xsph[1] / 5.0;

				for (k = 0; k < ndim; k++)
					tempup1[k] = tempup1[k] + xishu * pPar[j].mass * (pPar[j].vx[k] - pPar[i].vx[k]) * pParCell[i].wij[nj][3] / pPar[j].rho;

				tempdown = tempdown + pPar[j].mass * pParCell[i].wij[nj][3] / pPar[j].rho;
			}
		}

		if (pPar[i].type != 0 && pPar[i].type != 7)
		{
			if (fabs(tempdown) > 1.0e-7 && (pPar[i].type == 2 || pPar[i].type == 4)) {
				for (k = 0; k < ndim; k++) vx_temp[i][k] = pPar[i].vx[k] + tempup1[k] / tempdown;
			}
			else {
				for (k = 0; k < ndim; k++) vx_temp[i][k] = pPar[i].vx[k] + tempup1[k];
			}
			
		}
	}

#pragma omp parallel for schedule(static) private(k)
	for (i = 0; i < ntotal; i++)
	{
		if (pPar[i].type != 0 && pPar[i].type != 7) {
			for (k = 0; k < ndim; k++)
				pPar[i].vx[k] = vx_temp[i][k];
		}
	}
}
