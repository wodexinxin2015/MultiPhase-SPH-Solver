/************************************************************************************
Multi-Phase and Parallelized SPH program
--parallelized and optimized for SMP (Symmetric multiprocessing ) architecture.
--Copyright (c) 2018-2023, Weijie ZHANG, GeoHohai, Hohai University.
************************************************************************************/

#include <math.h>
#include "Class_NNPS.h"
#include "Class_particle.h"
#include "Class_ParticleProperty.h"
#include "Class_Problem.h"
#include "Header_Parameters.h"
#include "Class_Functions.h"

/* constitutive model */
/*elastic model*/
void clStraStre_Fun::elastic(double (*dept)[6], const Particle &pStrePar, const Para_Soil &pSoil)
{
	double dlam, dmu, dmu2;
	double e0, sm, sigct[6];
	double zkapa = pSoil.kapa;
	double poi = pSoil.poi;
	double eps0 = pSoil.eps0;
	int ki;

	for (ki = 0; ki < 6; ki++)
	{
		sigct[ki] = -pStrePar.sig[ki];
	}

	sm = (sigct[0] + sigct[1] + sigct[2]) * 0.3333333333;

	if (sm < 1.0)
		sm = 1.0;

	if (pSoil.e > 1.0)
		e0 = pSoil.e;
	else
	{
		e0 = 3 * (1 - 2 * poi) * (1 + eps0) * fabs(sm) / zkapa;
	}

	dlam = e0 * poi / (1.0 + poi) / (1.0 - 2.0 * poi);
	dmu = e0 * 0.5 / (1.0 + poi);
	dmu2 = dmu * 2.0;

	dept[0][0] = dlam + dmu2;
	dept[1][1] = dlam + dmu2;
	dept[2][2] = dlam + dmu2;
	dept[0][1] = dlam;
	dept[0][2] = dlam;
	dept[1][0] = dlam;
	dept[1][2] = dlam;
	dept[2][0] = dlam;
	dept[2][1] = dlam;
	dept[3][3] = dmu;
	dept[4][4] = dmu;
	dept[5][5] = dmu;
}

//elastic model for structure
void clStraStre_Fun::elastic_struct(double (*dept)[6], const Para_Structure &pSoil)
{
	double dlam, dmu, dmu2;
	double e0 = pSoil.e;
	double poi = pSoil.poi;

	e0 = pSoil.e;
	dlam = e0 * poi / (1.0 + poi) / (1.0 - 2.0 * poi);
	dmu = e0 * 0.5 / (1.0 + poi);
	dmu2 = dmu * 2.0;

	dept[0][0] = dlam + dmu2;
	dept[1][1] = dlam + dmu2;
	dept[2][2] = dlam + dmu2;
	dept[0][1] = dlam;
	dept[0][2] = dlam;
	dept[1][0] = dlam;
	dept[1][2] = dlam;
	dept[2][0] = dlam;
	dept[2][1] = dlam;
	dept[3][3] = dmu;
	dept[4][4] = dmu;
	dept[5][5] = dmu;
}

/*viscous elastic model*/
void clStraStre_Fun::viscous_elastic(Particle *pPar, const Para_Soil &pSoil)
{
	double k, epkk, g;
	double e0, sm, sigct[6];
	double zkapa = pSoil.kapa;
	double poi = pSoil.poi;
	double eps0 = pSoil.eps0;
	int ki;

	for (ki = 0; ki < 6; ki++)
	{
		sigct[ki] = -pPar->sig[ki];
	}

	sm = (sigct[0] + sigct[1] + sigct[2]) * 0.3333333333;

	if (sm < 1.0)
		sm = 1.0;

	if (pSoil.e > 1.0)
		e0 = pSoil.e;
	else
	{
		e0 = 3 * (1 - 2 * poi) * (1 + eps0) * fabs(sm) / zkapa;
	}

	k = e0 * 0.3333333333 / (1.0 - 2.0 * poi);
	g = pSoil.c;

	epkk = pPar->deps[0] + pPar->deps[1] + pPar->deps[2];

	pPar->dsig[0] = 2.0 * g * (pPar->deps[0] - epkk * 0.3333333333) + k * epkk;
	pPar->dsig[1] = 2.0 * g * (pPar->deps[1] - epkk * 0.3333333333) + k * epkk;
	pPar->dsig[2] = 2.0 * g * (pPar->deps[2] - epkk * 0.3333333333) + k * epkk;

	pPar->dsig[3] = 2.0 * g * pPar->deps[3];
	pPar->dsig[4] = 2.0 * g * pPar->deps[4];
	pPar->dsig[5] = 2.0 * g * pPar->deps[5];
}

/* unsaturated sub loading cam-clay model */
void clStraStre_Fun::sub_camclay(double (*dept)[6], const Particle &pStrePar, const Para_Soil &pSoil, double *adps)
{
	int ki;
	double dlam, dmu, dmu2, sjj;
	double e0, poi, ann, zmf, zramda, zkapa, eps0, roue;
	double ymf, sigct[6], sm, hee, sx, sy, sz, dx, dy, dz, dxy, dyz, dzx;
	double tt1, tt2, ttt, f11, f22, f33, f12, f23, f31, hdd, fkk, fij, cc7, d, df, dg, ramda;
	double ww1, ww2, ww3, ww4, ww5, ww0, zz1, zz2, zz3, zz4, zz5, zz0;

	ann = pSoil.ann;
	poi = pSoil.poi;
	zmf = pSoil.zmf;
	zramda = pSoil.ramda;
	zkapa = pSoil.kapa;
	eps0 = pSoil.eps0;

	for (ki = 0; ki < 6; ki++)
	{
		sigct[ki] = -pStrePar.sig[ki];
	}
	ymf = 3.0 * (zmf - 1.0) / (zmf + 2.0) / pSoil.fs;
	sm = (sigct[0] + sigct[1] + sigct[2]) * 0.3333333333;

	if (sm < 1.0)
		sm = 1.0;

	if (pSoil.e > 1.0)
		e0 = pSoil.e;
	else
	{
		e0 = 3 * (1 - 2 * poi) * (1 + eps0) * fabs(sm) / zkapa;
	}

	hee = (zramda - zkapa) / (1.0 + eps0);

	dlam = e0 * poi / (1.0 + poi) / (1.0 - 2.0 * poi);
	dmu = e0 * 0.5 / (1.0 + poi);
	dmu2 = dmu * 2.0;

	dept[0][0] = dlam + dmu2;
	dept[1][1] = dlam + dmu2;
	dept[2][2] = dlam + dmu2;
	dept[0][1] = dlam;
	dept[0][2] = dlam;
	dept[1][0] = dlam;
	dept[1][2] = dlam;
	dept[2][0] = dlam;
	dept[2][1] = dlam;
	dept[3][3] = dmu;
	dept[4][4] = dmu;
	dept[5][5] = dmu;

	sx = sigct[0];
	sy = sigct[1];
	sz = sigct[2];
	dx = sx - sm;
	dy = sy - sm;
	dz = sz - sm;
	dxy = sigct[3];
	dyz = sigct[4];
	dzx = sigct[5];
	sjj = 0.5 * (dx * dx + dy * dy + dz * dz) + dxy * dxy + dyz * dyz + dzx * dzx;

	tt1 = 1.0 / (ymf * ymf * sm * sm + 3.0 * sjj);
	tt2 = 2.0 * ymf * ymf * sm * tt1 - 1.0 / sm;
	f11 = tt2 / 3.0 + dx * tt1 * 3.0;
	f22 = tt2 / 3.0 + dy * tt1 * 3.0;
	f33 = tt2 / 3.0 + dz * tt1 * 3.0;
	f12 = 3.0 * dxy * tt1;
	f23 = 3.0 * dyz * tt1;
	f31 = 3.0 * dzx * tt1;

	roue = (1 + eps0) * log(pSoil.sme0 / sm);
	if (roue < 0.0)
		roue = 0.0;
	ttt = ann * roue * roue;
	hdd = f11 + f22 + f33 + ttt / sm;

	fkk = f11 + f22 + f33;
	fij = 2.0 * (f11 * f11 + f22 * f22 + f33 * f33) + f12 * f12 + f23 * f23 + f31 * f31;
	cc7 = dlam * (fkk * fkk) + dmu * fij + hdd / hee;
	d = 1.0 / cc7;

	df = dlam * fkk;
	dg = dlam * fkk;

	ramda = (dmu2 * f11 + dlam * fkk) * pStrePar.deps[0] + (dmu2 * f22 + dlam * fkk) * pStrePar.deps[1] + (dmu2 * f33 + dlam * fkk) * pStrePar.deps[2] + dmu * (f12 * pStrePar.deps[3] + f23 * pStrePar.deps[4] + f31 * pStrePar.deps[5]);

	if (ramda < 0.0)
	{
		ww0 = dg + dmu2 * f11;
		ww1 = dg + dmu2 * f22;
		ww2 = dg + dmu2 * f33;
		ww3 = dmu * f12;
		ww4 = dmu * f23;
		ww5 = dmu * f31;
		zz0 = df + dmu2 * f11;
		zz1 = df + dmu2 * f22;
		zz2 = df + dmu2 * f33;
		zz3 = dmu * f12;
		zz4 = dmu * f23;
		zz5 = dmu * f31;

		dept[0][0] = dept[0][0] - d * ww0 * zz0;
		dept[0][1] = dept[0][1] - d * ww0 * zz1;
		dept[0][2] = dept[0][2] - d * ww0 * zz2;
		dept[0][3] = dept[0][3] - d * ww0 * zz3;
		dept[0][4] = dept[0][4] - d * ww0 * zz4;
		dept[0][5] = dept[0][5] - d * ww0 * zz5;

		dept[1][0] = dept[1][0] - d * ww1 * zz0;
		dept[1][1] = dept[1][1] - d * ww1 * zz1;
		dept[1][2] = dept[1][2] - d * ww1 * zz2;
		dept[1][3] = dept[1][3] - d * ww1 * zz3;
		dept[1][4] = dept[1][4] - d * ww1 * zz4;
		dept[1][5] = dept[1][5] - d * ww1 * zz5;

		dept[2][0] = dept[2][0] - d * ww2 * zz0;
		dept[2][1] = dept[2][1] - d * ww2 * zz1;
		dept[2][2] = dept[2][2] - d * ww2 * zz2;
		dept[2][3] = dept[2][3] - d * ww2 * zz3;
		dept[2][4] = dept[2][4] - d * ww2 * zz4;
		dept[2][5] = dept[2][5] - d * ww2 * zz5;

		dept[3][0] = dept[3][0] - d * ww3 * zz0;
		dept[3][1] = dept[3][1] - d * ww3 * zz1;
		dept[3][2] = dept[3][2] - d * ww3 * zz2;
		dept[3][3] = dept[3][3] - d * ww3 * zz3;
		dept[3][4] = dept[3][4] - d * ww3 * zz4;
		dept[3][5] = dept[3][5] - d * ww3 * zz5;

		dept[4][0] = dept[4][0] - d * ww4 * zz0;
		dept[4][1] = dept[4][1] - d * ww4 * zz1;
		dept[4][2] = dept[4][2] - d * ww4 * zz2;
		dept[4][3] = dept[4][3] - d * ww4 * zz3;
		dept[4][4] = dept[4][4] - d * ww4 * zz4;
		dept[4][5] = dept[4][5] - d * ww4 * zz5;

		dept[5][0] = dept[5][0] - d * ww5 * zz0;
		dept[5][1] = dept[5][1] - d * ww5 * zz1;
		dept[5][2] = dept[5][2] - d * ww5 * zz2;
		dept[5][3] = dept[5][3] - d * ww5 * zz3;
		dept[5][4] = dept[5][4] - d * ww5 * zz4;
		dept[5][5] = dept[5][5] - d * ww5 * zz5;

		adps[0] = ramda * d * f11;
		adps[1] = ramda * d * f22;
		adps[2] = ramda * d * f33;
		adps[3] = ramda * d * f12;
		adps[4] = ramda * d * f23;
		adps[5] = ramda * d * f31;
	}
}

/*Drucker-Prager model*/ //drucker-prager model of DB Leaves
void clStraStre_Fun::dpmodel1(double (*dept)[6], const Particle &pStrePar, const Para_Soil &pSoil, double *adps)
{
	int ki, kj;
	double dlam, dmu, dmu2, sjj, sj, ramda, ttt, e0, eps0, poi;
	double sigct[6], sm, sx, sy, sz, dx, dy, dz, dxy, dyz, dzx;
	double f11, f22, f33, f12, f23, f31, hdd, fkk, fij, cc7, d, df, dg;
	double g11, g22, g33, g12, g23, g31;
	double ww1, ww2, ww3, ww4, ww5, ww0, zz1, zz2, zz3, zz4, zz5, zz0;
	double fai, ci, yalf, ykap;

	fai = atan(tan(pSoil.fai * 0.0174532925166667) / pSoil.fs);
	ci = pSoil.c / pSoil.fs;
	poi = pSoil.poi;
	eps0 = pSoil.eps0;

	for (ki = 0; ki < 6; ki++)
	{
		sigct[ki] = -pStrePar.sig[ki];
		for (kj = 0; kj < 6; kj++)
			dept[ki][kj] = 0.0;
	}
	sm = (sigct[0] + sigct[1] + sigct[2]) * 0.3333333333;

	if (sm < 1.0)
		sm = 1.0;
	e0 = 3 * (1 - 2 * poi) * (1 + eps0) * fabs(sm) * 50.0;

	dlam = e0 * poi / (1.0 + poi) / (1.0 - 2.0 * poi);
	dmu = e0 * 0.5 / (1.0 + poi);
	dmu2 = dmu * 2.0;

	dept[0][0] = dlam + dmu2;
	dept[1][1] = dlam + dmu2;
	dept[2][2] = dlam + dmu2;
	dept[0][1] = dlam;
	dept[0][2] = dlam;
	dept[1][0] = dlam;
	dept[1][2] = dlam;
	dept[2][0] = dlam;
	dept[2][1] = dlam;
	dept[3][3] = dmu;
	dept[4][4] = dmu;
	dept[5][5] = dmu;

	sx = sigct[0];
	sy = sigct[1];
	sz = sigct[2];
	dx = sx - sm;
	dy = sy - sm;
	dz = sz - sm;
	dxy = sigct[3];
	dyz = sigct[4];
	dzx = sigct[5];
	sjj = 0.5 * (dx * dx + dy * dy + dz * dz) + dxy * dxy + dyz * dyz + dzx * dzx;
	sj = sqrt(sjj);

	yalf = sin(fai) / sqrt(9.0 + 3.0 * sin(fai) * sin(fai));
	ykap = 3.0 * ci * cos(fai) / sqrt(9.0 + 3.0 * sin(fai) * sin(fai));

	if ((sj - 3.0 * yalf * sm - ykap) >= 0.0)
	{

		f11 = -yalf + dx / sj * 0.5;
		f22 = -yalf + dy / sj * 0.5;
		f33 = -yalf + dz / sj * 0.5;
		f12 = dxy / sj * 0.5;
		f23 = dyz / sj * 0.5;
		f31 = dzx / sj * 0.5;
		fkk = f11 + f22 + f33;

		g11 = dx / sj * 0.5;
		g22 = dy / sj * 0.5;
		g33 = dz / sj * 0.5;
		g12 = dxy / sj * 0.5;
		g23 = dyz / sj * 0.5;
		g31 = dzx / sj * 0.5;
		hdd = g11 + g22 + g33;

		fij = 2.0 * (f11 * g11 + f22 * g22 + f33 * g33) + f12 * g12 + f23 * g23 + f31 * g31;
		cc7 = dlam * (fkk * hdd) + dmu * fij;
		d = 1.0 / cc7;
		df = dlam * fkk;
		dg = dlam * hdd;

		ramda = (dmu2 * f11 - 3.0 * dlam * yalf) * pStrePar.deps[0] + (dmu2 * f22 - 3.0 * dlam * yalf) * pStrePar.deps[1] + (dmu2 * f33 - 3.0 * dlam * yalf) * pStrePar.deps[2] + dmu * (f12 * pStrePar.deps[3] + f23 * pStrePar.deps[4] + f31 * pStrePar.deps[5]);

		ttt = ramda * d;

		if (ttt < 0.0)
		{
			ww0 = dg + dmu2 * g11;
			ww1 = dg + dmu2 * g22;
			ww2 = dg + dmu2 * g33;
			ww3 = dmu * g12;
			ww4 = dmu * g23;
			ww5 = dmu * g31;
			zz0 = df + dmu2 * f11;
			zz1 = df + dmu2 * f22;
			zz2 = df + dmu2 * f33;
			zz3 = dmu * f12;
			zz4 = dmu * f23;
			zz5 = dmu * f31;

			dept[0][0] = dept[0][0] - d * ww0 * zz0;
			dept[0][1] = dept[0][1] - d * ww0 * zz1;
			dept[0][2] = dept[0][2] - d * ww0 * zz2;
			dept[0][3] = dept[0][3] - d * ww0 * zz3;
			dept[0][4] = dept[0][4] - d * ww0 * zz4;
			dept[0][5] = dept[0][5] - d * ww0 * zz5;

			dept[1][0] = dept[1][0] - d * ww1 * zz0;
			dept[1][1] = dept[1][1] - d * ww1 * zz1;
			dept[1][2] = dept[1][2] - d * ww1 * zz2;
			dept[1][3] = dept[1][3] - d * ww1 * zz3;
			dept[1][4] = dept[1][4] - d * ww1 * zz4;
			dept[1][5] = dept[1][5] - d * ww1 * zz5;

			dept[2][0] = dept[2][0] - d * ww2 * zz0;
			dept[2][1] = dept[2][1] - d * ww2 * zz1;
			dept[2][2] = dept[2][2] - d * ww2 * zz2;
			dept[2][3] = dept[2][3] - d * ww2 * zz3;
			dept[2][4] = dept[2][4] - d * ww2 * zz4;
			dept[2][5] = dept[2][5] - d * ww2 * zz5;

			dept[3][0] = dept[3][0] - d * ww3 * zz0;
			dept[3][1] = dept[3][1] - d * ww3 * zz1;
			dept[3][2] = dept[3][2] - d * ww3 * zz2;
			dept[3][3] = dept[3][3] - d * ww3 * zz3;
			dept[3][4] = dept[3][4] - d * ww3 * zz4;
			dept[3][5] = dept[3][5] - d * ww3 * zz5;
			dept[4][0] = dept[4][0] - d * ww4 * zz0;
			dept[4][1] = dept[4][1] - d * ww4 * zz1;
			dept[4][2] = dept[4][2] - d * ww4 * zz2;
			dept[4][3] = dept[4][3] - d * ww4 * zz3;
			dept[4][4] = dept[4][4] - d * ww4 * zz4;
			dept[4][5] = dept[4][5] - d * ww4 * zz5;

			dept[5][0] = dept[5][0] - d * ww5 * zz0;
			dept[5][1] = dept[5][1] - d * ww5 * zz1;
			dept[5][2] = dept[5][2] - d * ww5 * zz2;
			dept[5][3] = dept[5][3] - d * ww5 * zz3;
			dept[5][4] = dept[5][4] - d * ww5 * zz4;
			dept[5][5] = dept[5][5] - d * ww5 * zz5;

			adps[0] = ramda * d * f11;
			adps[1] = ramda * d * f22;
			adps[2] = ramda * d * f33;
			adps[3] = ramda * d * f12;
			adps[4] = ramda * d * f23;
			adps[5] = ramda * d * f31;
		}
	}
}

/*Drucker-Prager model*/ //drucker-prager model from Bui
void clStraStre_Fun::dpmodel_bui(double(*dept)[6], const Particle& pStrePar, const Para_Soil& pSoil, double* adps)
{
	int ki, kj;
	double dlam, dmu, dmu2, sjj, sj, ramda, ttt, e0, eps0, poi;
	double sigct[6], sm, sx, sy, sz, dx, dy, dz, dxy, dyz, dzx;
	double f11, f22, f33, f12, f23, f31, hdd, fkk, fij, cc7, d, df, dg;
	double g11, g22, g33, g12, g23, g31;
	double ww1, ww2, ww3, ww4, ww5, ww0, zz1, zz2, zz3, zz4, zz5, zz0;
	double fai, ci, yalf, ykap, fai_dil, yalf_dil;

	fai = atan(tan(pSoil.fai * 0.0174532925166667) / pSoil.fs);
	fai_dil = atan(tan(pSoil.beta * 0.0174532925166667) / pSoil.fs);
	ci = pSoil.c / pSoil.fs;
	poi = pSoil.poi;
	eps0 = pSoil.eps0;

	for (ki = 0; ki < 6; ki++)
	{
		sigct[ki] = -pStrePar.sig[ki];
		for (kj = 0; kj < 6; kj++)
			dept[ki][kj] = 0.0;
	}
	sm = (sigct[0] + sigct[1] + sigct[2]) * 0.3333333333;

	if (sm < 1.0)
		sm = 1.0;
	e0 = 3 * (1 - 2 * poi) * (1 + eps0) * fabs(sm) * 50.0;

	dlam = e0 * poi / (1.0 + poi) / (1.0 - 2.0 * poi);
	dmu = e0 * 0.5 / (1.0 + poi);
	dmu2 = dmu * 2.0;

	dept[0][0] = dlam + dmu2;
	dept[1][1] = dlam + dmu2;
	dept[2][2] = dlam + dmu2;
	dept[0][1] = dlam;
	dept[0][2] = dlam;
	dept[1][0] = dlam;
	dept[1][2] = dlam;
	dept[2][0] = dlam;
	dept[2][1] = dlam;
	dept[3][3] = dmu;
	dept[4][4] = dmu;
	dept[5][5] = dmu;

	sx = sigct[0];
	sy = sigct[1];
	sz = sigct[2];
	dx = sx - sm;
	dy = sy - sm;
	dz = sz - sm;
	dxy = sigct[3];
	dyz = sigct[4];
	dzx = sigct[5];
	sjj = 0.5 * (dx * dx + dy * dy + dz * dz) + dxy * dxy + dyz * dyz + dzx * dzx;
	sj = sqrt(sjj);

	yalf = sin(fai) / sqrt(9.0 + 3.0 * sin(fai) * sin(fai));
	ykap = 3.0 * ci * cos(fai) / sqrt(9.0 + 3.0 * sin(fai) * sin(fai));

	if ((sj - 3.0 * yalf * sm - ykap) >= 0.0)
	{

		f11 = -yalf + dx / sj * 0.5;
		f22 = -yalf + dy / sj * 0.5;
		f33 = -yalf + dz / sj * 0.5;
		f12 = dxy / sj * 0.5;
		f23 = dyz / sj * 0.5;
		f31 = dzx / sj * 0.5;
		fkk = f11 + f22 + f33;

		yalf_dil = sin(fai_dil) / sqrt(9.0 + 3.0 * sin(fai_dil) * sin(fai_dil));
		g11 = -yalf_dil + dx / sj * 0.5;
		g22 = -yalf_dil + dy / sj * 0.5;
		g33 = -yalf_dil + dz / sj * 0.5;
		g12 = dxy / sj * 0.5;
		g23 = dyz / sj * 0.5;
		g31 = dzx / sj * 0.5;
		hdd = g11 + g22 + g33;

		fij = 2.0 * (f11 * g11 + f22 * g22 + f33 * g33) + f12 * g12 + f23 * g23 + f31 * g31;
		cc7 = dlam * (fkk * hdd) + dmu * fij + hdd;
		d = 1.0 / cc7;
		df = dlam * fkk;
		dg = dlam * hdd;

		ramda = (dmu2 * f11 - 3.0 * dlam * yalf) * pStrePar.deps[0] + (dmu2 * f22 - 3.0 * dlam * yalf) * pStrePar.deps[1] + (dmu2 * f33 - 3.0 * dlam * yalf) * pStrePar.deps[2] + dmu * (f12 * pStrePar.deps[3] + f23 * pStrePar.deps[4] + f31 * pStrePar.deps[5]);

		ttt = ramda * d;

		if (ttt < 0.0)
		{
			ww0 = dg + dmu2 * g11;
			ww1 = dg + dmu2 * g22;
			ww2 = dg + dmu2 * g33;
			ww3 = dmu * g12;
			ww4 = dmu * g23;
			ww5 = dmu * g31;
			zz0 = df + dmu2 * f11;
			zz1 = df + dmu2 * f22;
			zz2 = df + dmu2 * f33;
			zz3 = dmu * f12;
			zz4 = dmu * f23;
			zz5 = dmu * f31;

			dept[0][0] = dept[0][0] - d * ww0 * zz0;
			dept[0][1] = dept[0][1] - d * ww0 * zz1;
			dept[0][2] = dept[0][2] - d * ww0 * zz2;
			dept[0][3] = dept[0][3] - d * ww0 * zz3;
			dept[0][4] = dept[0][4] - d * ww0 * zz4;
			dept[0][5] = dept[0][5] - d * ww0 * zz5;

			dept[1][0] = dept[1][0] - d * ww1 * zz0;
			dept[1][1] = dept[1][1] - d * ww1 * zz1;
			dept[1][2] = dept[1][2] - d * ww1 * zz2;
			dept[1][3] = dept[1][3] - d * ww1 * zz3;
			dept[1][4] = dept[1][4] - d * ww1 * zz4;
			dept[1][5] = dept[1][5] - d * ww1 * zz5;

			dept[2][0] = dept[2][0] - d * ww2 * zz0;
			dept[2][1] = dept[2][1] - d * ww2 * zz1;
			dept[2][2] = dept[2][2] - d * ww2 * zz2;
			dept[2][3] = dept[2][3] - d * ww2 * zz3;
			dept[2][4] = dept[2][4] - d * ww2 * zz4;
			dept[2][5] = dept[2][5] - d * ww2 * zz5;

			dept[3][0] = dept[3][0] - d * ww3 * zz0;
			dept[3][1] = dept[3][1] - d * ww3 * zz1;
			dept[3][2] = dept[3][2] - d * ww3 * zz2;
			dept[3][3] = dept[3][3] - d * ww3 * zz3;
			dept[3][4] = dept[3][4] - d * ww3 * zz4;
			dept[3][5] = dept[3][5] - d * ww3 * zz5;
			dept[4][0] = dept[4][0] - d * ww4 * zz0;
			dept[4][1] = dept[4][1] - d * ww4 * zz1;
			dept[4][2] = dept[4][2] - d * ww4 * zz2;
			dept[4][3] = dept[4][3] - d * ww4 * zz3;
			dept[4][4] = dept[4][4] - d * ww4 * zz4;
			dept[4][5] = dept[4][5] - d * ww4 * zz5;

			dept[5][0] = dept[5][0] - d * ww5 * zz0;
			dept[5][1] = dept[5][1] - d * ww5 * zz1;
			dept[5][2] = dept[5][2] - d * ww5 * zz2;
			dept[5][3] = dept[5][3] - d * ww5 * zz3;
			dept[5][4] = dept[5][4] - d * ww5 * zz4;
			dept[5][5] = dept[5][5] - d * ww5 * zz5;

			adps[0] = ramda * d * f11;
			adps[1] = ramda * d * f22;
			adps[2] = ramda * d * f33;
			adps[3] = ramda * d * f12;
			adps[4] = ramda * d * f23;
			adps[5] = ramda * d * f31;
		}
	}
}

/*elastic plastic M-C model*/
void clStraStre_Fun::elastic_plastic(double (*dept)[6], const Particle &pStrePar, const Para_Soil &pSoil, double *adps)
{
	int ki, kj;
	double dlam, dmu, dmu2, sjj, sj, ramda, ttt, e0, eps0, poi;
	double sigct[6], sm, sx, sy, sz, dx, dy, dz, dxy, dyz, dzx;
	double f11, f22, f33, f12, f23, f31, hdd, fkk, fij, cc7, d, df, dg;
	double g11, g22, g33, g12, g23, g31;
	double ww1, ww2, ww3, ww4, ww5, ww0, zz1, zz2, zz3, zz4, zz5, zz0;
	double fai, ci, yalf, ykap;

	fai = atan(tan(pSoil.fai * 0.0174532925166667) / pSoil.fs);
	ci = pSoil.c / pSoil.fs;
	poi = pSoil.poi;
	eps0 = pSoil.eps0;

	for (ki = 0; ki < 6; ki++)
	{
		sigct[ki] = -pStrePar.sig[ki];
		for (kj = 0; kj < 6; kj++)
			dept[ki][kj] = 0.0;
	}
	sm = (sigct[0] + sigct[1] + sigct[2]) * 0.333333333;

	if (sm < 1.0)
		sm = 1.0;
	e0 = 3 * (1 - 2 * poi) * (1 + eps0) * fabs(sm) * 50.0;

	dlam = e0 * poi / (1.0 + poi) / (1.0 - 2.0 * poi);
	dmu = e0 * 0.5 / (1.0 + poi);
	dmu2 = dmu * 2.0;

	dept[0][0] = dlam + dmu2;
	dept[1][1] = dlam + dmu2;
	dept[2][2] = dlam + dmu2;
	dept[0][1] = dlam;
	dept[0][2] = dlam;
	dept[1][0] = dlam;
	dept[1][2] = dlam;
	dept[2][0] = dlam;
	dept[2][1] = dlam;
	dept[3][3] = dmu;
	dept[4][4] = dmu;
	dept[5][5] = dmu;

	sx = sigct[0];
	sy = sigct[1];
	sz = sigct[2];
	dx = sx - sm;
	dy = sy - sm;
	dz = sz - sm;
	dxy = sigct[3];
	dyz = sigct[4];
	dzx = sigct[5];
	sjj = 0.5 * (dx * dx + dy * dy + dz * dz) + dxy * dxy + dyz * dyz + dzx * dzx;
	sj = sqrt(sjj);

	yalf = sin(fai) / sqrt(9.0 + 3.0 * sin(fai) * sin(fai));
	ykap = 3.0 * ci * cos(fai) / sqrt(9.0 + 3.0 * sin(fai) * sin(fai));

	if ((sj - 3.0 * yalf * sm - ykap) >= 0.0)
	{

		f11 = -yalf + dx / sj * 0.5;
		f22 = -yalf + dy / sj * 0.5;
		f33 = -yalf + dz / sj * 0.5;
		f12 = dxy / sj * 0.5;
		f23 = dyz / sj * 0.5;
		f31 = dzx / sj * 0.5;
		fkk = f11 + f22 + f33;

		g11 = dx / sj * 0.5;
		g22 = dy / sj * 0.5;
		g33 = dz / sj * 0.5;
		g12 = dxy / sj * 0.5;
		g23 = dyz / sj * 0.5;
		g31 = dzx / sj * 0.5;
		hdd = g11 + g22 + g33;

		fij = 2.0 * (f11 * g11 + f22 * g22 + f33 * g33) + f12 * g12 + f23 * g23 + f31 * g31;
		cc7 = dlam * (fkk * hdd) + dmu * fij;
		d = 1.0 / cc7;
		df = dlam * fkk;
		dg = dlam * hdd;

		ramda = (dmu2 * f11 - 3.0 * dlam * yalf) * pStrePar.deps[0] + (dmu2 * f22 - 3.0 * dlam * yalf) * pStrePar.deps[1] + (dmu2 * f33 - 3.0 * dlam * yalf) * pStrePar.deps[2] + dmu * (f12 * pStrePar.deps[3] + f23 * pStrePar.deps[4] + f31 * pStrePar.deps[5]);

		ttt = ramda * d;

		if (ttt < 0.0)
		{
			ww0 = dg + dmu2 * g11;
			ww1 = dg + dmu2 * g22;
			ww2 = dg + dmu2 * g33;
			ww3 = dmu * g12;
			ww4 = dmu * g23;
			ww5 = dmu * g31;
			zz0 = df + dmu2 * f11;
			zz1 = df + dmu2 * f22;
			zz2 = df + dmu2 * f33;
			zz3 = dmu * f12;
			zz4 = dmu * f23;
			zz5 = dmu * f31;

			dept[0][0] = dept[0][0] - d * ww0 * zz0;
			dept[0][1] = dept[0][1] - d * ww0 * zz1;
			dept[0][2] = dept[0][2] - d * ww0 * zz2;
			dept[0][3] = dept[0][3] - d * ww0 * zz3;
			dept[0][4] = dept[0][4] - d * ww0 * zz4;
			dept[0][5] = dept[0][5] - d * ww0 * zz5;

			dept[1][0] = dept[1][0] - d * ww1 * zz0;
			dept[1][1] = dept[1][1] - d * ww1 * zz1;
			dept[1][2] = dept[1][2] - d * ww1 * zz2;
			dept[1][3] = dept[1][3] - d * ww1 * zz3;
			dept[1][4] = dept[1][4] - d * ww1 * zz4;
			dept[1][5] = dept[1][5] - d * ww1 * zz5;

			dept[2][0] = dept[2][0] - d * ww2 * zz0;
			dept[2][1] = dept[2][1] - d * ww2 * zz1;
			dept[2][2] = dept[2][2] - d * ww2 * zz2;
			dept[2][3] = dept[2][3] - d * ww2 * zz3;
			dept[2][4] = dept[2][4] - d * ww2 * zz4;
			dept[2][5] = dept[2][5] - d * ww2 * zz5;

			dept[3][0] = dept[3][0] - d * ww3 * zz0;
			dept[3][1] = dept[3][1] - d * ww3 * zz1;
			dept[3][2] = dept[3][2] - d * ww3 * zz2;
			dept[3][3] = dept[3][3] - d * ww3 * zz3;
			dept[3][4] = dept[3][4] - d * ww3 * zz4;
			dept[3][5] = dept[3][5] - d * ww3 * zz5;
			dept[4][0] = dept[4][0] - d * ww4 * zz0;
			dept[4][1] = dept[4][1] - d * ww4 * zz1;
			dept[4][2] = dept[4][2] - d * ww4 * zz2;
			dept[4][3] = dept[4][3] - d * ww4 * zz3;
			dept[4][4] = dept[4][4] - d * ww4 * zz4;
			dept[4][5] = dept[4][5] - d * ww4 * zz5;

			dept[5][0] = dept[5][0] - d * ww5 * zz0;
			dept[5][1] = dept[5][1] - d * ww5 * zz1;
			dept[5][2] = dept[5][2] - d * ww5 * zz2;
			dept[5][3] = dept[5][3] - d * ww5 * zz3;
			dept[5][4] = dept[5][4] - d * ww5 * zz4;
			dept[5][5] = dept[5][5] - d * ww5 * zz5;

			adps[0] = ramda * d * f11;
			adps[1] = ramda * d * f22;
			adps[2] = ramda * d * f33;
			adps[3] = ramda * d * f12;
			adps[4] = ramda * d * f23;
			adps[5] = ramda * d * f31;
		}
	}
}

void clStraStre_Fun::DPModel_check(Particle* pStrePar, const Para_Soil& pSoil, double dt) {
	int ki;
	double temp_sig[6], i_1, alpha_fai, k_c, rn;
	double dx, dy, dz, dxy, dyz, dzx, sj_2, sj;
	double fai = atan(tan(pSoil.fai * 0.0174532925166667) / pSoil.fs);
	double ci = pSoil.c / pSoil.fs;
	double sm = (pStrePar->sig[0] + pStrePar->sig[1] + pStrePar->sig[2]) / 3.0;

	dx = pStrePar->sig[0] - sm;
	dy = pStrePar->sig[1] - sm;
	dz = pStrePar->sig[2] - sm;
	dxy = pStrePar->sig[3];
	dyz = pStrePar->sig[4];
	dzx = pStrePar->sig[5];
	sj_2 = 0.5 * (dx * dx + dy * dy + dz * dz) + dxy * dxy + dyz * dyz + dzx * dzx;
	sj = sqrt(sj_2);

	for (ki = 0; ki < 6; ki++)
	{
		temp_sig[ki] = pStrePar->sig[ki] + pStrePar->dsig[ki] * dt;
	}
	i_1 = temp_sig[0] + temp_sig[1] + temp_sig[2];

	alpha_fai = sin(fai) / sqrt(9.0 + 3.0 * sin(fai) * sin(fai));
	k_c = 3.0 * ci * cos(fai) / sqrt(9.0 + 3.0 * sin(fai) * sin(fai));
	rn = (-alpha_fai * i_1 + k_c) / sj;

	if (rn < 1.0) { //Stress-scaling back procedure
		temp_sig[0] = rn * (temp_sig[0] - i_1 * 0.3333333) + i_1 * 0.3333333;
		temp_sig[1] = rn * (temp_sig[1] - i_1 * 0.3333333) + i_1 * 0.3333333;
		temp_sig[2] = rn * (temp_sig[2] - i_1 * 0.3333333) + i_1 * 0.3333333;
		temp_sig[3] = rn * temp_sig[3];
		temp_sig[4] = rn * temp_sig[4];
		temp_sig[5] = rn * temp_sig[5];
	}

	for (ki = 0; ki < 6; ki++)
	{
		pStrePar->dsig[ki] = (temp_sig[ki] - pStrePar->sig[ki]) / dt;
	}
}

void clStraStre_Fun::MCModel_check(Particle* pStrePar, const Para_Soil& pSoil, double dt) {
	int ki;
	double temp_sig[6], i_1, alpha_fai, k_c, rn;
	double dx, dy, dz, dxy, dyz, dzx, sj_2, sj;
	double fai = atan(tan(pSoil.fai * 0.0174532925166667) / pSoil.fs);
	double ci = pSoil.c / pSoil.fs;
	double sm = (pStrePar->sig[0] + pStrePar->sig[1] + pStrePar->sig[2]) / 3.0;

	dx = pStrePar->sig[0] - sm;
	dy = pStrePar->sig[1] - sm;
	dz = pStrePar->sig[2] - sm;
	dxy = pStrePar->sig[3];
	dyz = pStrePar->sig[4];
	dzx = pStrePar->sig[5];
	sj_2 = 0.5 * (dx * dx + dy * dy + dz * dz) + dxy * dxy + dyz * dyz + dzx * dzx;
	sj = sqrt(sj_2);

	for (ki = 0; ki < 6; ki++)
	{
		temp_sig[ki] = pStrePar->sig[ki] + pStrePar->dsig[ki] * dt;
	}
	i_1 = temp_sig[0] + temp_sig[1] + temp_sig[2];

	alpha_fai = sin(fai) / sqrt(9.0 + 3.0 * sin(fai) * sin(fai));
	k_c = 3.0 * ci * cos(fai) / sqrt(9.0 + 3.0 * sin(fai) * sin(fai));
	rn = (-alpha_fai * i_1 + k_c) / sj;

	if (rn < 1.0) { //Stress-scaling back procedure
		temp_sig[0] = rn * (temp_sig[0] - i_1 * 0.3333333) + i_1 * 0.3333333;
		temp_sig[1] = rn * (temp_sig[1] - i_1 * 0.3333333) + i_1 * 0.3333333;
		temp_sig[2] = rn * (temp_sig[2] - i_1 * 0.3333333) + i_1 * 0.3333333;
		temp_sig[3] = rn * temp_sig[3];
		temp_sig[4] = rn * temp_sig[4];
		temp_sig[5] = rn * temp_sig[5];
	}

	for (ki = 0; ki < 6; ki++)
	{
		pStrePar->dsig[ki] = (temp_sig[ki] - pStrePar->sig[ki]) / dt;
	}
}