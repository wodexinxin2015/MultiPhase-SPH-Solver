/************************************************************************************
Multi-Phase and Parallelized SPH program
--parallelized and optimized for SMP (Symmetric multiprocessing ) architecture.
--Copyright (c) 2018-2023, Weijie ZHANG, GeoHohai, Hohai University.
************************************************************************************/
#include <cmath>
#include "Header_Option.h"
#include "Class_Functions.h"

//vibration acceleration
void clAcce_Fun::vibration_acceleration(Particle *pPar, const Para_Vibra &pPVib, const Para_Pro &pPPro, int cn)
{

	int i, k;
	int ntotal, ndim, nstep;
	double time1, time2, accvib[3];
	double max_acc;

	ntotal = pPPro.ntotal;
	ndim = pPPro.ndim;

	time1 = pPPro.t - pPVib.sttime;
	time2 = pPPro.t - pPVib.edtime;

	if (time1 >= 0.0 && time2 <= 0.0)
	{
		nstep = (int)(time1 / pPVib.ttime * (double)pPVib.nsteps);
		if (nstep > pPVib.nsteps - 1) 
			nstep = pPVib.nsteps - 1;
		accvib[0] = pPVib.loads[nstep][0];
		accvib[1] = pPVib.loads[nstep][1];
		accvib[2] = pPVib.loads[nstep][2];

		max_acc = fmax(0.0, fabs(accvib[0]));
		max_acc = fmax(max_acc, fabs(accvib[1]));
		max_acc = fmax(max_acc, fabs(accvib[2]));

		//applying vibration loads
		if (max_acc < 1000.0) {
#pragma omp parallel for schedule(static) private(k)
			for (i = 0; i < ntotal; i++)
			{
				if (pPar[i].type == 0 && (pPar[i].matype == 3 || pPar[i].matype == 5))
				{
					for (k = 0; k < ndim; k++)
					{
						pPar[i].interf[k] = accvib[k];
					}
				}
			}
		}
		else {
#pragma omp parallel for schedule(static) private(k)
			for (i = 0; i < ntotal; i++)
			{

				if (pPar[i].type == 0 && (pPar[i].matype == 3 || pPar[i].matype == 5))
				{
					for (k = 0; k < ndim; k++)
					{
						pPar[i].interf[k] = 0.0;
						pPar[i].interfss[k] = 0.0;
					}
				}
			}
		}
	}
	else
	{
#pragma omp parallel for schedule(static) private(k)
		for (i = 0; i < ntotal; i++)
		{

			if (pPar[i].type == 0 && (pPar[i].matype == 3 || pPar[i].matype == 5))
			{
				for (k = 0; k < ndim; k++)
				{
					pPar[i].interf[k] = 0.0;
					pPar[i].interfss[k] = 0.0;
				}
			}
		}
	}
}

//vibration loads for all particles
void clAcce_Fun::vibration_velocity(Particle *pPar, const Para_Vibra &pPVib, const Para_Pro &pPPro, int cn)
{

	int i, k;
	int ntotal, ndim, nstep;
	double time1, time2, velvib[3];
	double max_vel;

	ntotal = pPPro.ntotal;
	ndim = pPPro.ndim;

	time1 = pPPro.t - pPVib.sttime;
	time2 = pPPro.t - pPVib.edtime;

	if (time1 >= 0.0 && time2 <= 0.0)
	{
		nstep = (int)(time1 / pPVib.ttime * (double)pPVib.nsteps);
		if (nstep > (pPVib.nsteps - 1))
			nstep = pPVib.nsteps - 1;
		velvib[0] = pPVib.loads[nstep][0];
		velvib[1] = pPVib.loads[nstep][1];
		velvib[2] = pPVib.loads[nstep][2];

		max_vel = fmax(0.0, fabs(velvib[0]));
		max_vel = fmax(max_vel, fabs(velvib[1]));
		max_vel = fmax(max_vel, fabs(velvib[2]));

		//applying vibration loads
		if (max_vel < 1000.0) {
#pragma omp parallel for schedule(static) private(k)
			for (i = 0; i < ntotal; i++)
			{

				if (pPar[i].type == 0 && (pPar[i].matype == 4 || pPar[i].matype == 5))
				{
					for (k = 0; k < ndim; k++)
					{
						pPar[i].interfss[k] = velvib[k];
					}
				}
			}
		}
		else {
#pragma omp parallel for schedule(static) private(k)
			for (i = 0; i < ntotal; i++)
			{

				if (pPar[i].type == 0 && (pPar[i].matype == 4 || pPar[i].matype == 5))
				{
					for (k = 0; k < ndim; k++)
					{
						pPar[i].interf[k] = 0.0;
						pPar[i].interfss[k] = 0.0;
					}
				}
			}
		}
	}
	else
	{
#pragma omp parallel for schedule(static) private(k)
		for (i = 0; i < ntotal; i++)
		{

			if (pPar[i].type == 0 && (pPar[i].matype == 4 || pPar[i].matype == 5))
			{
				for (k = 0; k < ndim; k++)
				{
					pPar[i].interf[k] = 0.0;
					pPar[i].interfss[k] = 0.0;
				}
			}
		}
	}
}

//damping effect for the dynamic problems
void clAcce_Fun::soil_damping_vibration(Particle* pPar, Par_Cell* pParCell, const Para_Soil* pParti_ConsPara,
	StiffMat* pParStiff, const Para_Pro& pPPro)
{
	int i, j, k, nc, nj;
	int ntotal, ndim;
	double cd, alpha, beta, stiff, dr, damp[3], vx1[3], vx2[3], dt, tempdown;

	ntotal = pPPro.ntotal;
	ndim = pPPro.ndim;
	dr= pPPro.dr;
	dt= pPPro.dt;

#pragma omp parallel for schedule(static) private(j, k, nc, nj, cd, alpha, beta, stiff, damp, vx1, vx2, tempdown)
	for (i = 0; i < ntotal; i++)
	{
		if (pPar[i].type == 2 || pPar[i].type == 4) {
			//initialization
			cd = 0.0;
			for (k = 0; k < ndim; k++) damp[k] = 0.0;
			tempdown = 0.0;

			//calculate cd 
			stiff = fmax(0.0, pParStiff[i].depp0[0][0]);
			stiff = fmax(stiff, pParStiff[i].depp0[1][1]);
			stiff = fmax(stiff, pParStiff[i].depp0[2][2]);
			alpha = pParti_ConsPara[i].damp1;
			beta = pParti_ConsPara[i].damp2;
			cd = alpha + beta * sqrt(stiff / pPar[i].rhop / pPar[i].hl / pPar[i].hl);

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
					tempdown = tempdown + pPar[j].mass / pPar[j].rho * pParCell[i].wij[nj][3];
				}
			}

			//refecting to acceleration
			for (k = 0; k < ndim; k++)
			{
				if (fabs(tempdown) < 1.0e-7) tempdown = 1.0e-7;

				vx1[k] = pPar[i].vxp[k] + 0.5 * (pPar[i].axp[k] + pPar[i].ax[k]) * dt;
				vx2[k] = pPar[i].vxp[k] + 0.5 * (pPar[i].axp[k] + pPar[i].ax[k] - damp[k] / tempdown) * dt;
				if (vx1[k] * vx2[k] <= 0.0)
				{
					pPar[i].ax[k] = -2.0 * pPar[i].vxp[k] / dt - pPar[i].axp[k];
				}
				else
					pPar[i].ax[k] -= damp[k] / tempdown;
			}
		}
	}
}
