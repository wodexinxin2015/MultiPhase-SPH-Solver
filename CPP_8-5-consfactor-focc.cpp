/************************************************************************************
Multi-Phase and Parallelized SPH program
--parallelized and optimized for SMP(Symmetric multiprocessing) architecture.
--Copyright(c) 2016 - 2018, Weijie ZHANG, GeoHohai, Hohai University.
************************************************************************************/
#include <stdio.h>
#include <math.h>
#include "Class_Functions.h"
#include "Header_Parameters.h"
#include "Header_Option.h"

//calculating the stiffniss matrix [D]: dsig=[D]*deps
void clStraStre_Fun::fractional_mcc(double (*dept)[6], double (*cijkl)[6], Particle &pPar, const Para_Soil &pSoil)
{
	int ki, kj, kk;
	double sigct[6], sij[6], sm, q, s_val, skl, smn[3][3], sita;
	double poi, h1, h2, ramda, mc, alpha, beta, c, zeta, etao, g0;
	double p0, pa, m, ita, df, dg, pc, qc, e, mf, order;
	double mij[6], nij[6], hp, pijkl[6][6], total[6][6], adjoint[6][6];
	double g, k, e0, det_a, dlam, dmu, dmu2;
	double ramda_cc, ymf, tt1, tt2, f11, f22, f33, f12, f23, f31, fkk;

	//setting e, p0 and pa
	e = pPar.e;
	p0 = 98000;
	pa = 98000;

	//setting parameters
	poi = pSoil.poi;
	h1 = pSoil.damp1;
	h2 = pSoil.damp2;
	ramda = pSoil.ramda;
	mc = pSoil.zmf;
	alpha = pSoil.ann;
	beta = pSoil.bnn;
	c = pSoil.c;
	zeta = pSoil.zet;
	etao = pSoil.etao;
	g0 = pSoil.g0;

	//setting stress
	for (ki = 0; ki < 6; ki++)
	{
		sigct[ki] = -pPar.sig[ki];

		for (kj = 0; kj < 6; kj++)
		{
			cijkl[ki][kj] = 0.0;
			pijkl[ki][kj] = 0.0;
		}
	}

	//calculating p, sij and q
	sm = (sigct[0] + sigct[1] + sigct[2]) / 0.33333333;
	for (ki = 0; ki < 6; ki++)
	{
		if (ki < 3)
			sij[ki] = sigct[ki] - sm;
		else
			sij[ki] = sigct[ki];
	}
	q = 0.5 * (sij[0] * sij[0] + sij[1] * sij[1] + sij[2] * sij[2]) + sij[3] * sij[3] + sij[4] * sij[4] + sij[5] * sij[5];
	q = sqrt(3.0 * q);

	ita = q / sm;
	if (ita == 0.0)
	{
		ita = 0.000001;
		q = ita * sm;
	}
	//elastic parameter
	g = g0 * (2.97 - e) * (2.97 - e) / (1.0 + e) * sqrt(sm * pa);
	k = 2.0 * (1 + poi) * 0.33333333 / (1.0 - 2.0 * poi) * g;
	e0 = 2.0 * g * (1.0 + poi);

	//calculating ramda_cc from traditional cma-clay for the loading criteron
	ramda_cc = 0.0;
	ymf = mc;

	dlam = e0 * poi / (1.0 + poi) / (1.0 - 2.0 * poi);
	dmu = e0 * 0.5 / (1.0 + poi);
	dmu2 = dmu * 2.0;

	tt1 = 1.0 / (ymf * ymf * sm * sm + q * q);
	tt2 = 2.0 * ymf * ymf * sm * tt1 - 1.0 / sm;

	f11 = tt2 * 0.33333333 + sij[0] * tt1 * 3.0;
	f22 = tt2 * 0.33333333 + sij[1] * tt1 * 3.0;
	f33 = tt2 * 0.33333333 + sij[2] * tt1 * 3.0;
	f12 = 3.0 * sij[3] * tt1;
	f23 = 3.0 * sij[4] * tt1;
	f31 = 3.0 * sij[5] * tt1;

	fkk = f11 + f22 + f33;

	ramda_cc = (dmu2 * f11 + dlam * fkk) * pPar.deps[0] + (dmu2 * f22 + dlam * fkk) * pPar.deps[1] + (dmu2 * f33 + dlam * fkk) * pPar.deps[2] + dmu * (f12 * pPar.deps[3] + f23 * pPar.deps[4] + f31 * pPar.deps[5]);

	//calculating Lode's angle
	smn[0][0] = sij[0];
	smn[0][1] = sij[3];
	smn[0][2] = sij[5];

	smn[1][0] = sij[3];
	smn[1][1] = sij[1];
	smn[1][2] = sij[4];

	smn[2][0] = sij[5];
	smn[2][1] = sij[4];
	smn[2][2] = sij[2];

	s_val = 0.0;
	for (ki = 0; ki < 3; ki++)
	{
		for (kj = 0; kj < 3; kj++)
		{
			for (kk = 0; kk < 3; kk++)
			{
				s_val = s_val + smn[ki][kj] * smn[kj][kk] * smn[kk][ki];
			}
		}
	}
	skl = 0.0;
	for (ki = 0; ki < 3; ki++)
	{
		for (kj = 0; kj < 3; kj++)
		{
			skl = skl + smn[ki][kj] * smn[kj][ki];
		}
	}
	if (skl == 0.0)
		skl = 0.000001;

	sita = 9.0 * s_val * 0.5 / pow(1.5 * skl, 1.5);
	if (sita < -1.0)
		sita = -1.0;
	else if (sita > 1.0)
		sita = 1.0;
	sita = asin(sita) * 0.33333333;

	//calculating M
	c = pow(c, 4);
	m = 2 * c / (1 + c + (1 - c) * sin(3.0 * sita));
	m = pow(m, 0.25);
	m = mc * m;

	//calculating pc and qc
	pc = pa * pow((etao - e) / ramda, 1 / zeta);
	qc = q + m * (sm - pc);

	//calculating df and dg
	df = (sm - p0 / 2) * mc * mc / q;
	dg = ((sm - pc) + (2 - alpha) * (pc - sm * 0.5 - q * q / mc / mc / sm / 2)) / ((q - qc) + (2 - alpha) * qc);
	dg = pow(m, alpha + 1) * dg;

	//calculating mij and nij
	mij[0] = 1.0 / sqrt(1.0 + df * df) * (3.0 * sij[0] * 0.5 / q + df * 0.33333333);
	mij[1] = 1.0 / sqrt(1.0 + df * df) * (3.0 * sij[1] * 0.5 / q + df * 0.33333333);
	mij[2] = 1.0 / sqrt(1.0 + df * df) * (3.0 * sij[2] * 0.5 / q + df * 0.33333333);
	mij[3] = 1.0 / sqrt(1.0 + df * df) * (3.0 * sij[3] * 0.5 / q);
	mij[4] = 1.0 / sqrt(1.0 + df * df) * (3.0 * sij[4] * 0.5 / q);
	mij[5] = 1.0 / sqrt(1.0 + df * df) * (3.0 * sij[5] * 0.5 / q);

	nij[0] = 1.0 / sqrt(1.0 + dg * dg) * (3.0 * sij[0] * 0.5 / q + dg * 0.33333333);
	nij[1] = 1.0 / sqrt(1.0 + dg * dg) * (3.0 * sij[1] * 0.5 / q + dg * 0.33333333);
	nij[2] = 1.0 / sqrt(1.0 + dg * dg) * (3.0 * sij[2] * 0.5 / q + dg * 0.33333333);
	nij[3] = 1.0 / sqrt(1.0 + dg * dg) * (3.0 * sij[3] * 0.5 / q);
	nij[4] = 1.0 / sqrt(1.0 + dg * dg) * (3.0 * sij[4] * 0.5 / q);
	nij[5] = 1.0 / sqrt(1.0 + dg * dg) * (3.0 * sij[5] * 0.5 / q);

	//calculating hardening parameter
	order = beta * (e - etao + ramda * pow(sm / pa, zeta));
	mf = mc * pow(2.71828, -order);
	hp = (h1 - h2 * e) * (1 - e) * g * (mf - ita) / ita * pow(2.71828, order);

	//calculating elastic matrix

	cijkl[0][0] = 1.0 / e0;
	cijkl[1][1] = 1.0 / e0;
	cijkl[2][2] = 1.0 / e0;

	cijkl[0][1] = -poi / e0;
	cijkl[0][2] = -poi / e0;
	cijkl[1][0] = -poi / e0;
	cijkl[1][2] = -poi / e0;
	cijkl[2][0] = -poi / e0;
	cijkl[2][1] = -poi / e0;

	cijkl[3][3] = 1.0 / g;
	cijkl[4][4] = 1.0 / g;
	cijkl[5][5] = 1.0 / g;

	//calculating plastic matrix
	for (ki = 0; ki < 6; ki++)
	{
		for (kj = 0; kj < 6; kj++)
		{
			pijkl[ki][kj] = nij[ki] * mij[kj] / hp;
		}
	}

	//forming total matrix
	if (ramda_cc < 0.0)
	{
		for (ki = 0; ki < 6; ki++)
		{
			for (kj = 0; kj < 6; kj++)
			{
				total[ki][kj] = cijkl[ki][kj] + pijkl[ki][kj];
			}
		}
	}
	else
	{
		for (ki = 0; ki < 6; ki++)
		{
			for (kj = 0; kj < 6; kj++)
			{
				total[ki][kj] = cijkl[ki][kj];
			}
		}
	}

	det_a = detA(total, 6);
	getAStart(total, 6, adjoint);

	for (ki = 0; ki < 6; ki++)
	{
		for (kj = 0; kj < 6; kj++)
		{
			if (det_a != 0.0)
				dept[ki][kj] = adjoint[ki][kj] / det_a;
			else
			{
				det_a = detA(cijkl, 6);
				getAStart(cijkl, 6, adjoint);
				dept[ki][kj] = adjoint[ki][kj] / det_a;
			}
		}
	}
}

//calculating plastic strain
void clStraStre_Fun::fractional_plastic(double *dsig, double *deps, double *adps, double cijkl[6][6])
{
	int ki, kj;
	double eeps[6];

	for (ki = 0; ki < 6; ki++)
	{
		eeps[ki] = 0.0;
		for (kj = 0; kj < 6; kj++)
		{
			eeps[ki] = eeps[ki] + cijkl[ki][kj] * dsig[kj];
		}
	}

	for (ki = 0; ki < 6; ki++)
		adps[ki] = deps[ki] - eeps[ki];
}

void clStraStre_Fun::Display(double arr[][6], int n)
{
	int i, j;
	printf(" ");
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			printf("%7.5lf  ", arr[i][j]);
		}
		printf("\n");
	}
}

double clStraStre_Fun::detA(double arr[][6], int n)
{
	double bak[6][6];
	int i, j, k, c;
	double sum = 0;

	if (n == 1)
	{
		return arr[0][0];
	}

	for (i = 0; i < n; i++)
	{
		for (j = 1; j < n; j++)
		{
			for (c = 0, k = 0; k < n; k++)
			{
				if (k == i)
				{
					continue;
				}
				bak[j - 1][c++] = arr[j][k];
			}
		}
		sum += (i % 2 == 0 ? 1 : -1) * arr[0][i] * detA(bak, n - 1);
	}
	return sum;
}

//adjoint matrix
void clStraStre_Fun::getAStart(double arcs[6][6], int n, double ans[6][6])
{
	if (n == 1)
	{
		ans[0][0] = 1.0;
		return;
	}

	int i, j, k, t;
	double temp[6][6];

	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			for (k = 0; k < n - 1; k++)
			{
				for (t = 0; t < n - 1; t++)
				{
					temp[k][t] = arcs[k >= i ? k + 1 : k][t >= j ? t + 1 : t];
				}
			}

			ans[j][i] = detA(temp, n - 1);

			if ((i + j) % 2 == 1)
			{
				ans[j][i] = -ans[j][i];
			}
		}
	}
}
