/************************************************************************************
Multi-Phase and Parallelized SPH program
--parallelized and optimized for SMP (Symmetric multiprocessing ) architecture.
--Copyright (c) 2018-2023, Weijie ZHANG, GeoHohai, Hohai University.
************************************************************************************/

#include <math.h>
#include "Header_Parameters.h"
#include "Header_Option.h"
#include "Class_Functions.h"

/* the calculation for the interaction force between fluid phase and solid phase*/
//only for the Darcy's law and laminar flow
/*Biot, 1941 and Prevost, 1979*/
void clInterFor_Fun::interact1(Particle *pPar, Par_Cell *pParCell, const Para_Inter &pInter,
							   const Para_Fluid *pVFluid, const Para_Pro &pPPro, Para_GF pgf, int cn)
{

	int i, j, k, nc, nj, ntotal, ndim;
	double tempup[3], inprod, dx[3], dvx[3], vx1[3], vx2[3];

	double xishu = pInter.interc;
	double dt = pPPro.dt;
	double g_acc;

	ntotal = pPPro.ntotal;
	ndim = pPPro.ndim;
	g_acc = sqrt(pgf.gx * pgf.gx + pgf.gy * pgf.gy + pgf.gz * pgf.gz);

#pragma omp parallel for schedule(static) private(j, k, nc, nj, dx, dvx, inprod, vx1, vx2, tempup)
	for (i = 0; i < ntotal; i++)
	{

		for (j = 0; j < 3; j++)
		{
			tempup[j] = 0.0;
		}

		nc = pParCell[i].ninflu;
		for (nj = 1; nj <= nc; nj++)
		{
			j = pParCell[i].influ[nj];

			if (pPar[i].type != pPar[j].type)
			{
				inprod = 0.0;
				for (k = 0; k < ndim; k++)
				{
					dx[k] = pPar[i].xp[k] - pPar[j].xp[k];
					dvx[k] = pPar[i].vxp[k] - pPar[j].vxp[k];
					inprod = inprod + dx[k] * dvx[k];
				}

				if (inprod < 0)
				{
					//saturation effect
					if (pPar[i].type == 2 && pPar[j].type == 1)
					{
						for (k = 0; k < ndim; k++)
						{
							tempup[k] = tempup[k] + pPar[i].satu * pPar[i].satu * pPar[i].porosity * pPar[i].porosity * pVFluid[j].dens * g_acc * dvx[k] / pPar[i].cop * pPar[j].mass / pPar[j].rho * pParCell[i].wij[nj][3];
						}
						//Mori and Uzuoka
					}
					else if (pPar[i].type == 1 && pPar[j].type == 2)
					{
						for (k = 0; k < ndim; k++)
						{
							tempup[k] = tempup[k] + pPar[j].satu * pPar[j].satu * pPar[j].porosity * pPar[j].porosity * pVFluid[i].dens * g_acc * dvx[k] / pPar[j].cop * pPar[j].mass / pPar[j].rho * pParCell[i].wij[nj][3];
						}
						//Mori and Uzuoka
					}
					else if (pPar[i].type == 2 && pPar[j].type == 3)
					{
						for (k = 0; k < ndim; k++)
						{
							tempup[k] = tempup[k] + (1 - pPar[i].satu) * (1 - pPar[i].satu) * pPar[i].porosity * pPar[i].porosity * pVFluid[j].dens * g_acc * dvx[k] * 2000 * pPar[j].mass / pPar[j].rho * pParCell[i].wij[nj][3];
						}
						//Mori and Uzuoka
					}
					else if (pPar[i].type == 3 && pPar[j].type == 2)
					{
						for (k = 0; k < ndim; k++)
						{
							tempup[k] = tempup[k] + (1 - pPar[j].satu) * (1 - pPar[j].satu) * pPar[j].porosity * pPar[j].porosity * pVFluid[i].dens * g_acc * dvx[k] * 2000 * pPar[j].mass / pPar[j].rho * pParCell[i].wij[nj][3];
						}
						//Mori and Uzuoka
					}
				}
			}
		}

		if (pPar[i].type != 0 && pPar[i].type != 7)
		{
			for (k = 0; k < ndim; k++)
			{
				pPar[i].interf[k] = tempup[k];
				vx1[k] = pPar[i].vxp[k] + 0.5 * (pPar[i].axp[k] + pPar[i].ax[k]) * dt;
				vx2[k] = pPar[i].vxp[k] + 0.5 * (pPar[i].axp[k] + pPar[i].ax[k] - pPar[i].interf[k] * xishu / pPar[i].rho) * dt;
				if (vx1[k] * vx2[k] <= 0.0)
				{
					pPar[i].interf[k] = 0.0;
					pPar[i].ax[k] = -2.0 * pPar[i].vxp[k] / dt - pPar[i].axp[k];
				}
				else
					pPar[i].ax[k] -= pPar[i].interf[k] * xishu / pPar[i].rho;
			}
		}
	}
}

/* the calculation for the interaction force between fluid phase and solid phase*/
//inlcuding the interaction forece of laminar flow and turbulent flow
/*from Tsuji et al. (Powder Technology) 1993 and Shi ZM et al. (Computers and Geotechnics) 2018*/
void clInterFor_Fun::interact2(Particle *pPar, Par_Cell *pParCell, const Para_Inter &pInter,
							   const Para_Soil *pParti_ConsPara, const Para_Fluid *pVFluid, const Para_Pro &pPPro, Para_GF pgf, int cn)
{

	int i, j, k, nc, nj, ntotal, ndim;
	double tempup[3], inprod, dx[3], dvx[3], vx1[3], vx2[3];
	double re, cd, mag_v, beta;

	double xishu = pInter.interc;
	double dt = pPPro.dt;
	double g_acc;

	ntotal = pPPro.ntotal;
	ndim = pPPro.ndim;
	g_acc = sqrt(pgf.gx * pgf.gx + pgf.gy * pgf.gy + pgf.gz * pgf.gz);

#pragma omp parallel for schedule(static) private(j, k, nc, nj, dx, dvx, inprod, vx1, vx2, tempup, re, cd, mag_v, beta)
	for (i = 0; i < ntotal; i++)
	{

		for (j = 0; j < 3; j++)
		{
			tempup[j] = 0.0;
		}

		nc = pParCell[i].ninflu;
		for (nj = 1; nj <= nc; nj++)
		{
			j = pParCell[i].influ[nj];

			if (pPar[i].type != pPar[j].type)
			{
				inprod = 0.0;
				mag_v = 0.0;
				for (k = 0; k < ndim; k++)
				{
					dx[k] = pPar[i].xp[k] - pPar[j].xp[k];
					dvx[k] = pPar[i].vxp[k] - pPar[j].vxp[k];
					mag_v = mag_v + dvx[k] * dvx[k];
					inprod = inprod + dx[k] * dvx[k];
				}
				mag_v = sqrt(mag_v);

				//saturation effect
				if (pPar[i].type == 2 && pPar[j].type == 1)
				{
					if (inprod < 0)
					{
						//calculating Renold number
						re = mag_v * pVFluid[j].dens * pPar[i].porosity * pParti_ConsPara[i].ds / pVFluid[j].vis;

						//calculating cd
						if (re > 1000.0)
							cd = 0.43;
						else if (fabs(re) > 0.0001)
							cd = 24.0 * (1.0 + 0.15 * pow(re, 0.687)) / re;
						else
							cd = 0.0;

						//calculating beta
						if (pPar[i].porosity <= 0.8)
						{
							beta = 150.0 * (1.0 - pPar[i].porosity) * (1.0 - pPar[i].porosity) / pPar[i].porosity / pPar[i].porosity / pParti_ConsPara[i].ds / pParti_ConsPara[i].ds * pVFluid[j].vis;
							beta = beta + 1.75 * (1.0 - pPar[i].porosity) / pParti_ConsPara[i].ds / pPar[i].porosity * pVFluid[j].dens * mag_v;
						}
						else
						{
							beta = 0.75 * cd * mag_v * pVFluid[j].dens * (1.0 - pPar[i].porosity) / pParti_ConsPara[i].ds * pow(pPar[i].porosity, -2.7);
						}

						//calculating interaction force
						for (k = 0; k < ndim; k++)
						{
							tempup[k] = tempup[k] + pPar[j].mass / pPar[j].rho * pParCell[i].wij[nj][3] * beta * (pPar[i].vxp[k] - pPar[j].vxp[k]);
						}
					}
				}
				else if (pPar[i].type == 1 && pPar[j].type == 2)
				{
					if (inprod < 0)
					{
						//calculating Renold number
						re = mag_v * pVFluid[i].dens * pPar[j].porosity * pParti_ConsPara[j].ds / pVFluid[i].vis;

						//calculating cd
						if (re > 1000.0)
							cd = 0.43;
						else if (fabs(re) > 0.0001)
							cd = 24.0 * (1.0 + 0.15 * pow(re, 0.687)) / re;
						else
							cd = 0.0;

						//calculating beta
						if (pPar[j].porosity <= 0.8)
						{
							beta = 150.0 * (1.0 - pPar[j].porosity) * (1.0 - pPar[j].porosity) / pPar[j].porosity / pPar[j].porosity / pParti_ConsPara[j].ds / pParti_ConsPara[j].ds * pVFluid[i].vis;
							beta = beta + 1.75 * (1.0 - pPar[j].porosity) / pParti_ConsPara[j].ds / pPar[j].porosity * pVFluid[i].dens * mag_v;
						}
						else
						{
							beta = 0.75 * cd * mag_v * pVFluid[i].dens * (1.0 - pPar[j].porosity) / pParti_ConsPara[j].ds * pow(pPar[j].porosity, -2.7);
						}

						//calculating interaction force
						for (k = 0; k < ndim; k++)
						{
							tempup[k] = tempup[k] + pPar[j].mass / pPar[j].rho * pParCell[i].wij[nj][3] * beta * (pPar[i].vxp[k] - pPar[j].vxp[k]);
						}
					}
				}
				else if (pPar[i].type == 2 && pPar[j].type == 3)
				{
					if (inprod < 0)
					{
						for (k = 0; k < ndim; k++)
						{
							tempup[k] = tempup[k] + (1 - pPar[i].satu) * (1 - pPar[i].satu) * pPar[i].porosity * pPar[i].porosity * pVFluid[j].dens * g_acc * dvx[k] * 2000 * pPar[j].mass / pPar[j].rho * pParCell[i].wij[nj][3];
						}
					}
				}
				else if (pPar[i].type == 3 && pPar[j].type == 2)
				{
					if (inprod < 0)
					{
						for (k = 0; k < ndim; k++)
						{
							tempup[k] = tempup[k] + (1 - pPar[j].satu) * (1 - pPar[j].satu) * pPar[j].porosity * pPar[j].porosity * pVFluid[i].dens * g_acc * dvx[k] * 2000 * pPar[j].mass / pPar[j].rho * pParCell[i].wij[nj][3];
						}
					}
				}
			}
		}

		if (pPar[i].type != 0 && pPar[i].type != 7)
		{
			for (k = 0; k < ndim; k++)
			{
				pPar[i].interf[k] = tempup[k];
				vx1[k] = pPar[i].vxp[k] + 0.5 * (pPar[i].axp[k] + pPar[i].ax[k]) * dt;
				vx2[k] = pPar[i].vxp[k] + 0.5 * (pPar[i].axp[k] + pPar[i].ax[k] - pPar[i].interf[k] * xishu / pPar[i].rho) * dt;
				if (vx1[k] * vx2[k] <= 0.0)
				{
					pPar[i].interf[k] = 0.0;
					pPar[i].ax[k] = -2.0 * pPar[i].vxp[k] / dt - pPar[i].axp[k];
				}
				else
					pPar[i].ax[k] -= pPar[i].interf[k] * xishu / pPar[i].rho;
			}
		}
	}
}

/* the calculation for the interaction force between fluid phase and solid phase*/
//inlcuding the interaction forece of laminar flow and turbulent flow
/*from Peng C. et al. (Computers and Geotechnics) 2017*/
void clInterFor_Fun::interact3(Particle *pPar, Par_Cell *pParCell, const Para_Inter &pInter,
							   const Para_Soil *pParti_ConsPara, const Para_Fluid *pVFluid, const Para_Pro &pPPro, Para_GF pgf, int cn)
{

	int i, j, k, nc, nj, ntotal, ndim;
	double tempup[3], inprod, dx[3], dvx[3], vx1[3], vx2[3];
	double re, cd, mag_v, beta;

	double xishu = pInter.interc;
	double dt = pPPro.dt;
	double g_acc;

	ntotal = pPPro.ntotal;
	ndim = pPPro.ndim;
	g_acc = sqrt(pgf.gx * pgf.gx + pgf.gy * pgf.gy + pgf.gz * pgf.gz);

#pragma omp parallel for schedule(static) private(j, k, nc, nj, dx, dvx, inprod, vx1, vx2, tempup, re, cd, mag_v, beta)
	for (i = 0; i < ntotal; i++)
	{

		for (j = 0; j < 3; j++)
		{
			tempup[j] = 0.0;
		}

		nc = pParCell[i].ninflu;
		for (nj = 1; nj <= nc; nj++)
		{
			j = pParCell[i].influ[nj];

			if (pPar[i].type != pPar[j].type)
			{
				inprod = 0.0;
				mag_v = 0.0;
				for (k = 0; k < ndim; k++)
				{
					dx[k] = pPar[i].xp[k] - pPar[j].xp[k];
					dvx[k] = pPar[i].vxp[k] - pPar[j].vxp[k];
					mag_v = mag_v + dvx[k] * dvx[k];
					inprod = inprod + dx[k] * dvx[k];
				}
				mag_v = sqrt(mag_v);

				//saturation effect
				if (pPar[i].type == 2 && pPar[j].type == 1)
				{
					if (inprod < 0)
					{
						//calculating Renold number
						re = mag_v * pVFluid[j].dens * pParti_ConsPara[i].ds / pVFluid[j].vis;

						//calculating beta
						if (re > 60.0)
						{
							cd = pPar[i].porosity * pPar[i].porosity * pPar[i].porosity * pParti_ConsPara[i].ds * pParti_ConsPara[i].ds * 0.0066666667 / (1.0 - pPar[i].porosity) / (1.0 - pPar[i].porosity);
							beta = pVFluid[j].vis / cd + 1.75 / sqrt(150.0) * pPar[j].rho * mag_v / sqrt(cd) / pow(pPar[i].porosity, 1.5);
						}
						else
						{
							cd = pPar[i].porosity * pPar[i].porosity * pPar[i].porosity * pParti_ConsPara[i].ds * pParti_ConsPara[i].ds * 0.0066666667 / (1.0 - pPar[i].porosity) / (1.0 - pPar[i].porosity);
							beta = pVFluid[j].vis / cd;
						}

						//calculating interaction force
						for (k = 0; k < ndim; k++)
						{
							tempup[k] = tempup[k] + pPar[j].mass / pPar[j].rho * pParCell[i].wij[nj][3] * beta * (pPar[i].vxp[k] - pPar[j].vxp[k]);
						}
					}
				}
				else if (pPar[i].type == 1 && pPar[j].type == 2)
				{
					if (inprod < 0)
					{
						//calculating Renold number
						re = mag_v * pVFluid[i].dens * pParti_ConsPara[j].ds / pVFluid[i].vis;

						//calculating beta
						if (re > 60.0)
						{
							cd = pPar[j].porosity * pPar[j].porosity * pPar[j].porosity * pParti_ConsPara[j].ds * pParti_ConsPara[j].ds * 0.0066666667 / (1.0 - pPar[j].porosity) / (1.0 - pPar[j].porosity);
							beta = pVFluid[i].vis / cd + 1.75 / sqrt(150.0) * pPar[i].rho * mag_v / sqrt(cd) / pow(pPar[j].porosity, 1.5);
						}
						else
						{
							cd = pPar[j].porosity * pPar[j].porosity * pPar[j].porosity * pParti_ConsPara[j].ds * pParti_ConsPara[j].ds * 0.0066666667 / (1.0 - pPar[j].porosity) / (1.0 - pPar[j].porosity);
							beta = pVFluid[i].vis / cd;
						}

						//calculating interaction force
						for (k = 0; k < ndim; k++)
						{
							tempup[k] = tempup[k] + pPar[j].mass / pPar[j].rho * pParCell[i].wij[nj][3] * beta * (pPar[i].vxp[k] - pPar[j].vxp[k]);
						}
					}
				}
				else if (pPar[i].type == 2 && pPar[j].type == 3)
				{
					if (inprod < 0)
					{
						for (k = 0; k < ndim; k++)
						{
							tempup[k] = tempup[k] + (1 - pPar[i].satu) * (1 - pPar[i].satu) * pPar[i].porosity * pPar[i].porosity * pVFluid[j].dens * g_acc * dvx[k] * 2000 * pPar[j].mass / pPar[j].rho * pParCell[i].wij[nj][3];
						}
					}
				}
				else if (pPar[i].type == 3 && pPar[j].type == 2)
				{
					if (inprod < 0)
					{
						for (k = 0; k < ndim; k++)
						{
							tempup[k] = tempup[k] + (1 - pPar[j].satu) * (1 - pPar[j].satu) * pPar[j].porosity * pPar[j].porosity * pVFluid[i].dens * g_acc * dvx[k] * 2000 * pPar[j].mass / pPar[j].rho * pParCell[i].wij[nj][3];
						}
					}
				}
			}
		}

		if (pPar[i].type != 0 && pPar[i].type != 7)
		{
			for (k = 0; k < ndim; k++)
			{
				pPar[i].interf[k] = tempup[k];
				vx1[k] = pPar[i].vxp[k] + 0.5 * (pPar[i].axp[k] + pPar[i].ax[k]) * dt;
				vx2[k] = pPar[i].vxp[k] + 0.5 * (pPar[i].axp[k] + pPar[i].ax[k] - pPar[i].interf[k] * xishu / pPar[i].rho) * dt;
				if (vx1[k] * vx2[k] <= 0.0)
				{
					pPar[i].interf[k] = 0.0;
					pPar[i].ax[k] = -2.0 * pPar[i].vxp[k] / dt - pPar[i].axp[k];
				}
				else
					pPar[i].ax[k] -= pPar[i].interf[k] * xishu / pPar[i].rho;
			}
		}
	}
}