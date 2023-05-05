/************************************************************************************
Multi-Phase and Parallelized SPH program
--parallelized and optimized for SMP (Symmetric multiprocessing ) architecture.
--Copyright (c) 2018-2023, Weijie ZHANG, GeoHohai, Hohai University.
************************************************************************************/

#include <math.h>
#include "Header_Parameters.h"
#include "Header_Option.h"
#include "Class_Functions.h"

/* strain rate --->stress rate: using constitutive model*/
void clStraStre_Fun::strain_to_stress(Particle *pPar, const Para_Soil *pParti_ConsPara, StiffMat *pParStiff, const Para_Fluid *pVFluid,
									  const Para_Structure &pStruct, const Para_Pro &pPPro, int cn, int flts, double alpha)
{
	int i, ki, kj;
	double dept[6][6], cijkl[6][6], erfs[6], depsmax, tao, fai;
	int ntotal;
	double viscos1, adps[6], deps[6], strain_min;
	double dt = pPPro.dt;

	ntotal = pPPro.ntotal;

#pragma omp parallel for schedule(static) private(ki, kj, dept, deps, erfs, viscos1, depsmax, tao, adps, cijkl, fai, strain_min)
	for (i = 0; i < ntotal; i++)
	{
		/*water and air*/
		if (pPar[i].type == 1)
		{
			strain_min = pVFluid[i].strain_min;
			depsmax = (pPar[i].deps[0] + pPar[i].deps[1] + pPar[i].deps[2]) / 3.0;
			pPar[i].eps[0] = pPar[i].deps[0] - depsmax;
			pPar[i].eps[1] = pPar[i].deps[1] - depsmax;
			pPar[i].eps[2] = pPar[i].deps[2] - depsmax;
			depsmax = 0.666666666666667 * (pPar[i].eps[0] * pPar[i].eps[0] + pPar[i].eps[1] * pPar[i].eps[1] 
				+ pPar[i].eps[2] * pPar[i].eps[2] + pPar[i].deps[3] * pPar[i].deps[3] * 0.25 + 
				pPar[i].deps[4] * pPar[i].deps[4] * 0.25 + pPar[i].deps[5] * pPar[i].deps[5] * 0.25);
			depsmax = sqrt(depsmax);
			//Bingham
			if (pVFluid[i].bhtype == 1)
			{
				fai = pVFluid[i].fai * 0.0174532925199;
				tao = pPar[i].pre * tan(fai) + pVFluid[i].c;
				if (depsmax < strain_min)
					viscos1 = pVFluid[i].vis + tao / strain_min;
				else
					viscos1 = pVFluid[i].vis + tao / depsmax;
			}
			else
				viscos1 = pVFluid[i].vis;

			deps[0] = 2.0 * pPar[i].deps[0] - 2.0 * (pPar[i].deps[0] + pPar[i].deps[1] + pPar[i].deps[2]) / 3.0;
			deps[1] = 2.0 * pPar[i].deps[1] - 2.0 * (pPar[i].deps[0] + pPar[i].deps[1] + pPar[i].deps[2]) / 3.0;
			deps[2] = 2.0 * pPar[i].deps[2] - 2.0 * (pPar[i].deps[0] + pPar[i].deps[1] + pPar[i].deps[2]) / 3.0;
			deps[3] = pPar[i].deps[3];
			deps[4] = pPar[i].deps[4];
			deps[5] = pPar[i].deps[5];

			for (ki = 0; ki < 6; ki++)
			{
				if (ki < 3)
					pPar[i].sig[ki] = -pPar[i].pre - viscos1 * deps[ki];
				else
					pPar[i].sig[ki] = -viscos1 * deps[ki];
			}
		}
		else if (pPar[i].type == 3)
		{

			viscos1 = pVFluid[i].vis;

			deps[0] = 2.0 * pPar[i].deps[0] - 2.0 * (pPar[i].deps[0] + pPar[i].deps[1] + pPar[i].deps[2]) / 3.0;
			deps[1] = 2.0 * pPar[i].deps[1] - 2.0 * (pPar[i].deps[0] + pPar[i].deps[1] + pPar[i].deps[2]) / 3.0;
			deps[2] = 2.0 * pPar[i].deps[2] - 2.0 * (pPar[i].deps[0] + pPar[i].deps[1] + pPar[i].deps[2]) / 3.0;
			deps[3] = pPar[i].deps[3];
			deps[4] = pPar[i].deps[4];
			deps[5] = pPar[i].deps[5];

			for (ki = 0; ki < 6; ki++)
			{
				if (ki < 3)
					pPar[i].sig[ki] = -pPar[i].pre - viscos1 * deps[ki];
				else
					pPar[i].sig[ki] = -viscos1 * deps[ki];
			}
		}
		else if (pPar[i].type == 2)
		{ /*soil*/
			for (ki = 0; ki < 6; ki++)
			{
				erfs[ki] = 0.0;
				adps[ki] = 0.0;
				for (kj = 0; kj < 6; kj++)
				{
					dept[ki][kj] = 0.0;
				}
			}

			/*transformed stress for the 3D Space Step 1*/
			if (flts >= 1)
				TSSMP_before(&pPar[i], alpha);

			/* Constitutive model */
			if (pPPro.l < pPPro.inip)
			{
				elastic(dept, pPar[i], pParti_ConsPara[i]);
				for (ki = 0; ki < 6; ki++)
				{
					for (kj = 0; kj < 6; kj++)
					{
						pPar[i].dsig[ki] += dept[ki][kj] * pPar[i].deps[kj];
						pParStiff[i].depp0[ki][kj] = dept[ki][kj];
					}
				}
			}
			else
			{

				if (pPar[i].matype == 0 || pPar[i].matype == 10 || pPar[i].matype == 20 || pPar[i].matype == 50)
				{ //EM
					elastic(dept, pPar[i], pParti_ConsPara[i]);
					for (ki = 0; ki < 6; ki++)
					{
						for (kj = 0; kj < 6; kj++)
						{
							pPar[i].dsig[ki] += dept[ki][kj] * pPar[i].deps[kj];
							pParStiff[i].depp0[ki][kj] = dept[ki][kj];
						}
					}
				}
				else if (pPar[i].matype == 1 || pPar[i].matype == 11 || pPar[i].matype == 21 || pPar[i].matype == 51)
				{ //SB
					sub_camclay(dept, pPar[i], pParti_ConsPara[i], adps);
					for (ki = 0; ki < 6; ki++)
					{
						for (kj = 0; kj < 6; kj++)
						{
							pPar[i].dsig[ki] += dept[ki][kj] * pPar[i].deps[kj];
							pParStiff[i].depp0[ki][kj] = dept[ki][kj];
						}
					}
				}
				else if (pPar[i].matype == 2 || pPar[i].matype == 12 || pPar[i].matype == 22 || pPar[i].matype == 52)
				{ //Drucker-Prager model of DB Leaves
					dpmodel1(dept, pPar[i], pParti_ConsPara[i], adps);
					for (ki = 0; ki < 6; ki++)
					{
						for (kj = 0; kj < 6; kj++)
						{
							pPar[i].dsig[ki] += dept[ki][kj] * pPar[i].deps[kj];
							pParStiff[i].depp0[ki][kj] = dept[ki][kj];
						}
					}
					DPModel_check(&pPar[i], pParti_ConsPara[i], dt);
				}
				else if (pPar[i].matype == 3 || pPar[i].matype == 13 || pPar[i].matype == 23 || pPar[i].matype == 53)
				{ //UMCC
					sub_camunsatu(dept, erfs, pPar[i], pParti_ConsPara[i], adps);
					for (ki = 0; ki < 6; ki++)
					{
						for (kj = 0; kj < 6; kj++)
						{
							pPar[i].dsig[ki] += dept[ki][kj] * pPar[i].deps[kj];
							pParStiff[i].depp0[ki][kj] = dept[ki][kj];
						}
					}
				}
				else if (pPar[i].matype == 4 || pPar[i].matype == 14 || pPar[i].matype == 24 || pPar[i].matype == 54)
				{ //AMCC or TSAMCC
					sub_camanisotropy(dept, erfs, pPar[i], pParti_ConsPara[i], adps);
					for (ki = 0; ki < 6; ki++)
					{
						for (kj = 0; kj < 6; kj++)
						{
							pPar[i].dsig[ki] += dept[ki][kj] * pPar[i].deps[kj];
							pParStiff[i].depp0[ki][kj] = dept[ki][kj];
						}
					}
				}
				else if (pPar[i].matype == 5 || pPar[i].matype == 15 || pPar[i].matype == 25 || pPar[i].matype == 55)
				{ //Drucker-Prager model from Bui
					dpmodel_bui(dept, pPar[i], pParti_ConsPara[i], adps);
					for (ki = 0; ki < 6; ki++)
					{
						for (kj = 0; kj < 6; kj++)
						{
							pPar[i].dsig[ki] += dept[ki][kj] * pPar[i].deps[kj];
							pParStiff[i].depp0[ki][kj] = dept[ki][kj];
						}
					}
					DPModel_check(&pPar[i], pParti_ConsPara[i], dt);
				}
				else if (pPar[i].matype == 6 || pPar[i].matype == 16 || pPar[i].matype == 26 || pPar[i].matype == 56)
				{
					elastic(dept, pPar[i], pParti_ConsPara[i]);
					for (ki = 0; ki < 6; ki++)
					{
						for (kj = 0; kj < 6; kj++)
						{
							pPar[i].dsig[ki] += dept[ki][kj] * pPar[i].deps[kj];
							pParStiff[i].depp0[ki][kj] = dept[ki][kj];
						}
					}
				}
				else if (pPar[i].matype == 7 || pPar[i].matype == 17 || pPar[i].matype == 27 || pPar[i].matype == 57)
				{ //VE
					viscous_elastic(&pPar[i], pParti_ConsPara[i]);
				}
				else if (pPar[i].matype == 8 || pPar[i].matype == 18 || pPar[i].matype == 28 || pPar[i].matype == 58)
				{ //focc
					fractional_mcc(dept, cijkl, pPar[i], pParti_ConsPara[i]);
					for (ki = 0; ki < 6; ki++)
					{
						for (kj = 0; kj < 6; kj++)
						{
							pPar[i].dsig[ki] += dept[ki][kj] * pPar[i].deps[kj];
							pParStiff[i].depp0[ki][kj] = dept[ki][kj];
						}
					}
				}
				else if (pPar[i].matype == 9 || pPar[i].matype == 19 || pPar[i].matype == 29 || pPar[i].matype == 59)
				{ //M-C failure criteria
					elastic_plastic(dept, pPar[i], pParti_ConsPara[i], adps);
					for (ki = 0; ki < 6; ki++)
					{
						for (kj = 0; kj < 6; kj++)
						{
							pPar[i].dsig[ki] += dept[ki][kj] * pPar[i].deps[kj];
							pParStiff[i].depp0[ki][kj] = dept[ki][kj];
						}
					}
				}
			}

			/*Correction of stress rate*/
			pPar[i].dsig[0] = pPar[i].dsig[0] 
				+ (pPar[i].weps[0][0] * pPar[i].sig[0] + pPar[i].weps[0][1] * pPar[i].sig[3] + pPar[i].weps[0][2] * pPar[i].sig[5]) 
				- (pPar[i].sig[0] * pPar[i].weps[0][0] + pPar[i].sig[3] * pPar[i].weps[1][0] + pPar[i].sig[5] * pPar[i].weps[2][0]);
			pPar[i].dsig[1] = pPar[i].dsig[1] 
				+ (pPar[i].weps[1][0] * pPar[i].sig[3] + pPar[i].weps[1][1] * pPar[i].sig[1] + pPar[i].weps[1][2] * pPar[i].sig[4]) 
				- (pPar[i].sig[3] * pPar[i].weps[0][1] + pPar[i].sig[1] * pPar[i].weps[1][1] + pPar[i].sig[4] * pPar[i].weps[2][1]);
			pPar[i].dsig[2] = pPar[i].dsig[2] 
				+ (pPar[i].weps[2][0] * pPar[i].sig[5] + pPar[i].weps[2][1] * pPar[i].sig[4] + pPar[i].weps[2][2] * pPar[i].sig[2]) 
				- (pPar[i].sig[5] * pPar[i].weps[0][2] + pPar[i].sig[4] * pPar[i].weps[1][2] + pPar[i].sig[2] * pPar[i].weps[2][2]);
			pPar[i].dsig[3] = pPar[i].dsig[3] 
				+ (pPar[i].weps[0][0] * pPar[i].sig[3] + pPar[i].weps[0][1] * pPar[i].sig[1] + pPar[i].weps[0][2] * pPar[i].sig[4]) 
				- (pPar[i].sig[0] * pPar[i].weps[0][1] + pPar[i].sig[3] * pPar[i].weps[1][1] + pPar[i].sig[5] * pPar[i].weps[2][1]);
			pPar[i].dsig[4] = pPar[i].dsig[4] 
				+ (pPar[i].weps[1][0] * pPar[i].sig[5] + pPar[i].weps[1][1] * pPar[i].sig[4] + pPar[i].weps[1][2] * pPar[i].sig[2]) 
				- (pPar[i].sig[3] * pPar[i].weps[0][2] + pPar[i].sig[1] * pPar[i].weps[1][2] + pPar[i].sig[4] * pPar[i].weps[2][2]);
			pPar[i].dsig[5] = pPar[i].dsig[5] 
				+ (pPar[i].weps[0][0] * pPar[i].sig[5] + pPar[i].weps[0][1] * pPar[i].sig[4] + pPar[i].weps[0][2] * pPar[i].sig[2]) 
				- (pPar[i].sig[0] * pPar[i].weps[0][2] + pPar[i].sig[3] * pPar[i].weps[1][2] + pPar[i].sig[5] * pPar[i].weps[2][2]);

			if (pPar[i].matype == 8 || pPar[i].matype == 18 || pPar[i].matype == 28)
			{ //focc
				fractional_plastic(pPar[i].dsig, pPar[i].deps, adps, cijkl);
			}

			/*transformed stress for the 3D Space Step 2*/
			if (flts >= 1)
				TSSMP_after(&pPar[i]);

			/*Update of stress*/
			for (ki = 0; ki < 6; ki++)
			{
				pPar[i].adps[ki] = pPar[i].adps[ki] + adps[ki] * dt;
				pPar[i].eps[ki] = pPar[i].eps[ki] + pPar[i].deps[ki] * dt;
				pPar[i].sig[ki] = pPar[i].sig[ki] + pPar[i].dsig[ki] * dt + erfs[ki];
			}

			pPar[i].deq();
			pPar[i].shearstrain();
			pPar[i].poro_cal(dt);
		}
		else if (pPar[i].type == 4)
		{ /*structure*/
			for (ki = 0; ki < 6; ki++)
			{
				for (kj = 0; kj < 6; kj++)
				{
					dept[ki][kj] = 0.0;
				}
			}

			elastic_struct(dept, pStruct);
			for (ki = 0; ki < 6; ki++)
			{
				for (kj = 0; kj < 6; kj++)
				{
					pPar[i].dsig[ki] += dept[ki][kj] * pPar[i].deps[kj];
					pParStiff[i].depp0[ki][kj] = dept[ki][kj];
				}
			}

			/*Correction of stress rate*/
			pPar[i].dsig[0] = pPar[i].dsig[0] + (pPar[i].weps[0][0] * pPar[i].sig[0] + pPar[i].weps[0][1] * pPar[i].sig[3] + pPar[i].weps[0][2] * pPar[i].sig[5]) - (pPar[i].sig[0] * pPar[i].weps[0][0] + pPar[i].sig[3] * pPar[i].weps[1][0] + pPar[i].sig[5] * pPar[i].weps[2][0]);

			pPar[i].dsig[1] = pPar[i].dsig[1] + (pPar[i].weps[1][0] * pPar[i].sig[3] + pPar[i].weps[1][1] * pPar[i].sig[1] + pPar[i].weps[1][2] * pPar[i].sig[4]) - (pPar[i].sig[3] * pPar[i].weps[0][1] + pPar[i].sig[1] * pPar[i].weps[1][1] + pPar[i].sig[4] * pPar[i].weps[2][1]);

			pPar[i].dsig[2] = pPar[i].dsig[2] + (pPar[i].weps[2][0] * pPar[i].sig[5] + pPar[i].weps[2][1] * pPar[i].sig[4] + pPar[i].weps[2][2] * pPar[i].sig[2]) - (pPar[i].sig[5] * pPar[i].weps[0][2] + pPar[i].sig[4] * pPar[i].weps[1][2] + pPar[i].sig[2] * pPar[i].weps[2][2]);

			pPar[i].dsig[3] = pPar[i].dsig[3] + (pPar[i].weps[0][0] * pPar[i].sig[3] + pPar[i].weps[0][1] * pPar[i].sig[1] + pPar[i].weps[0][2] * pPar[i].sig[4]) - (pPar[i].sig[0] * pPar[i].weps[0][1] + pPar[i].sig[3] * pPar[i].weps[1][1] + pPar[i].sig[5] * pPar[i].weps[2][1]);

			pPar[i].dsig[4] = pPar[i].dsig[4] + (pPar[i].weps[1][0] * pPar[i].sig[5] + pPar[i].weps[1][1] * pPar[i].sig[4] + pPar[i].weps[1][2] * pPar[i].sig[2]) - (pPar[i].sig[3] * pPar[i].weps[0][2] + pPar[i].sig[1] * pPar[i].weps[1][2] + pPar[i].sig[4] * pPar[i].weps[2][2]);

			pPar[i].dsig[5] = pPar[i].dsig[5] + (pPar[i].weps[0][0] * pPar[i].sig[5] + pPar[i].weps[0][1] * pPar[i].sig[4] + pPar[i].weps[0][2] * pPar[i].sig[2]) - (pPar[i].sig[0] * pPar[i].weps[0][2] + pPar[i].sig[3] * pPar[i].weps[1][2] + pPar[i].sig[5] * pPar[i].weps[2][2]);

			/*Update of stress*/
			for (ki = 0; ki < 6; ki++)
			{
				pPar[i].eps[ki] = pPar[i].eps[ki] + pPar[i].deps[ki] * dt;
				pPar[i].sig[ki] = pPar[i].sig[ki] + pPar[i].dsig[ki] * dt;
			}

			pPar[i].deq();
			pPar[i].shearstrain();
			pPar[i].poro_cal(dt);
		}
	}
}

/*moisture characteristic curve for unsaturated soil: sr->suct*/
double moisturesuct(double satu, double dsatu, const Para_Soil &Soil)
{
	double temp, suctn;

	if ((satu - Soil.sitar) < 0.0001)
		temp = tan(pi * (Soil.sitas - satu - 0.0001) / (Soil.sitas - Soil.sitar) / 2.0);
	temp = tan(pi * (Soil.sitas - satu) / (Soil.sitas - Soil.sitar) / 2.0);

	if (dsatu > 0.0) //wetting curve
		suctn = log(1 + exp(Soil.c2 * Soil.sw) * temp) / Soil.c2;
	else
		suctn = log(1 + exp(Soil.c1 * Soil.sd) * temp) / Soil.c1;

	return suctn;
}

/*moisture characteristic curve for unsaturated soil: suct->sr*/
double moisturesr(double suct, double dsuct, const Para_Soil &Soil)
{
	double temp, sr;

	temp = 2 * (Soil.sitas - Soil.sitar) / pi;

	if (dsuct < 0.0) //wetting curve
		sr = Soil.sitas - temp * atan((exp(Soil.c2 * suct) - 1) / exp(Soil.c2 * Soil.sw));
	else
		sr = Soil.sitas - temp * atan((exp(Soil.c1 * suct) - 1) / exp(Soil.c1 * Soil.sd));

	return sr;
}

void clInterFor_Fun::cal_permeability(Particle *pPar, const Para_Inter &pPInter, const Para_Soil *pParti_ConsPara, const Para_Pro &pPPro, int cn)
{
	/* Calculation of  the permeability */

	int i, ntotal;
	ntotal = pPPro.ntotal;

	/* calculation of the pressure for water and air*/
#pragma omp parallel for schedule(static)
	for (i = 0; i < ntotal; i++)
	{
		if (pPar[i].type == 2)
		{
			if (pPInter.permtype == 0)
				pPar[i].cop = pParti_ConsPara[i].cop;
			//coefficient of permeability determined by the constitutive model
			else if (pPInter.permtype == 1)
				pPar[i].cop = pPInter.pers;
			//homogeneous dike
			else if (pPInter.permtype == 2)
				pPar[i].permeablel(pParti_ConsPara[i].sitas, pParti_ConsPara[i].sitar, pPInter.pers, pPInter.perl);
			//homogeneous dike with variable permeability changing with saturation
			else if (pPInter.permtype == 3)
				pPar[i].permeabletl(pPInter.pers, pPInter.perl, pPInter.perm);
			//heterogeneous dike with three layers
		}
	}
}

//stress regularizaiton by Peng et al. 2019
void clStraStre_Fun::stress_regu(Particle* pPar, Par_Cell* pParCell, const Para_Pro& pPPro)
{
	int i, j, nc, nj, ntotal;
	double temp_sig_up[6], temp_down;

	ntotal = pPPro.ntotal;

#pragma omp parallel for schedule(static) private(j, nc, nj, temp_sig_up, temp_down)
	for (i = 0; i < ntotal; i++) {
		for (j = 0; j < 6; j++) temp_sig_up[j] = 0.0;
		temp_down = 0.0;

		nc = pParCell[i].ninflu;
		for (nj = 1; nj <= nc; nj++)
		{
			j = pParCell[i].influ[nj];
			if (pPar[i].type == pPar[j].type) {
				temp_sig_up[0] = temp_sig_up[0] + pPar[j].mass / pPar[j].rho * pParCell[i].wij[nj][3] * pPar[j].sig[0];
				temp_sig_up[1] = temp_sig_up[1] + pPar[j].mass / pPar[j].rho * pParCell[i].wij[nj][3] * pPar[j].sig[1];
				temp_sig_up[2] = temp_sig_up[2] + pPar[j].mass / pPar[j].rho * pParCell[i].wij[nj][3] * pPar[j].sig[2];
				temp_sig_up[3] = temp_sig_up[3] + pPar[j].mass / pPar[j].rho * pParCell[i].wij[nj][3] * pPar[j].sig[3];
				temp_sig_up[4] = temp_sig_up[4] + pPar[j].mass / pPar[j].rho * pParCell[i].wij[nj][3] * pPar[j].sig[4];
				temp_sig_up[5] = temp_sig_up[5] + pPar[j].mass / pPar[j].rho * pParCell[i].wij[nj][3] * pPar[j].sig[5];

				temp_down = temp_down + pPar[j].mass / pPar[j].rho * pParCell[i].wij[nj][3];
			}
		}
		if (fabs(temp_down) > 1.0e-6) {
			for (j = 0; j < 6; j++) pPar[i].sig[j] = temp_sig_up[j] / temp_down;
		}
	}
}

void clStraStre_Fun::strain_regu(Particle* pPar, Par_Cell* pParCell, const Para_Pro& pPPro)
{
	int i, j, nc, nj, ntotal;
	double temp_eps_up[6], temp_down;

	ntotal = pPPro.ntotal;

#pragma omp parallel for schedule(static) private(j, nc, nj, temp_eps_up, temp_down)
	for (i = 0; i < ntotal; i++) {
		for (j = 0; j < 6; j++) temp_eps_up[j] = 0.0;
		temp_down = 0.0;

		nc = pParCell[i].ninflu;
		for (nj = 1; nj <= nc; nj++)
		{
			j = pParCell[i].influ[nj];
			if (pPar[i].type == pPar[j].type) {
				temp_eps_up[0] = temp_eps_up[0] + pPar[j].mass / pPar[j].rho * pParCell[i].wij[nj][3] * pPar[j].eps[0];
				temp_eps_up[1] = temp_eps_up[1] + pPar[j].mass / pPar[j].rho * pParCell[i].wij[nj][3] * pPar[j].eps[1];
				temp_eps_up[2] = temp_eps_up[2] + pPar[j].mass / pPar[j].rho * pParCell[i].wij[nj][3] * pPar[j].eps[2];
				temp_eps_up[3] = temp_eps_up[3] + pPar[j].mass / pPar[j].rho * pParCell[i].wij[nj][3] * pPar[j].eps[3];
				temp_eps_up[4] = temp_eps_up[4] + pPar[j].mass / pPar[j].rho * pParCell[i].wij[nj][3] * pPar[j].eps[4];
				temp_eps_up[5] = temp_eps_up[5] + pPar[j].mass / pPar[j].rho * pParCell[i].wij[nj][3] * pPar[j].eps[5];

				temp_down = temp_down + pPar[j].mass / pPar[j].rho * pParCell[i].wij[nj][3];
			}
		}

		if (fabs(temp_down) > 1.0e-6) {
			for (j = 0; j < 6; j++) pPar[i].eps[j] = temp_eps_up[j] / temp_down;
		}
	}
}

//free surface check and stress correction
void clStraStre_Fun::free_surface_check(Particle* pPar, Par_Cell* pParCell, const Para_Pro& pPPro, double max, double min) {
	int i, j, k, nc, nj, ntotal;
	double color, threshold_min, threshold_max, incline, xishu;

	ntotal = pPPro.ntotal;
	threshold_max = max;
	threshold_min = min;
	incline = threshold_max - threshold_min;

#pragma omp parallel for schedule(static) private(j, k, nc, nj, color, xishu)
	for (i = 0; i < ntotal; i++) {
		//initialization of color value
		color = 0.0;

		//color summation of neighboring particles
		nc = pParCell[i].ninflu;
		for (nj = 1; nj <= nc; nj++)
		{
			j = pParCell[i].influ[nj];

			if (pPar[i].type == pPar[j].type && (pPar[i].type != 0 && pPar[i].type != 7) ) {
				color = color + pPar[j].mass / pPar[j].rho * 1.0 * pParCell[i].wij[nj][3];
			}
			else if(pPar[i].type != pPar[j].type && (pPar[i].type > 0 && pPar[i].type < 5)) {
				if(pPar[j].type == 0)
					color = color + pPar[j].mass / pPar[j].rho * 1.0 * pParCell[i].wij[nj][3];
			}
		}

		//comparison between color and threshold, stress correction
		if (pPar[i].type == 1) {
			if (color < threshold_min / 2.0) {
				for (k = 0; k < 3; k++) pPar[i].sig[k] = -1.0;
				for (k = 3; k < 6; k++) pPar[i].sig[k] = 0.0;
			}
			else if (color >= threshold_min / 2.0 && color < threshold_max / 2.0) {
				xishu = (2.0 * color - threshold_min) / incline;
				for (k = 0; k < 6; k++) pPar[i].sig[k] = xishu * pPar[i].sig[k];
			}
		}
		else if (pPar[i].type == 2) {
			if (color < threshold_min) {
				for (k = 0; k < 3; k++) pPar[i].sig[k] = -1.0;
				for (k = 3; k < 6; k++) pPar[i].sig[k] = 0.0;
			}
			else if (color >= threshold_min && color < threshold_max) {
				xishu = (color - threshold_min) / incline * 1.0;
				for (k = 0; k < 6; k++) pPar[i].sig[k] = xishu * pPar[i].sig[k];
			}
			pPar[i].stress_adjust();
			pPar[i].deq();
		}
	}
}