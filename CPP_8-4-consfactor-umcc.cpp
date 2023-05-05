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
/* unsaturated sub loading cam-clay model */
void clStraStre_Fun::sub_camunsatu(double (*dept)[6], double *erfs, const Particle &pStrePar, const Para_Soil &pSoil, double *adps)
{
	int ki;
	double dlam, dmu, dmu2, sjj;
	double e0, roue, rous;
	double ymf, sigct[6], sm, hee, sx, sy, sz, dx, dy, dz, dxy, dyz, dzx;
	double tt1, tt2, ttt, f11, f22, f33, f12, f23, f31, hdd, fkk, fij, cc7, d, df, dg, ramda, omake;
	double ww1, ww2, ww3, ww4, ww5, ww0, zz1, zz2, zz3, zz4, zz5, zz0, qq;

	double zkapa = pSoil.kapa;
	double poi = pSoil.poi;
	double eps0 = pSoil.eps0;

	qq = (pSoil.epsr - pSoil.epss) / (pSoil.sitas - pSoil.sitar);

	for (ki = 0; ki < 6; ki++)
	{
		sigct[ki] = -pStrePar.sig[ki];
	}
	ymf = 3.0 * (pSoil.zmf - 1.0) / (pSoil.zmf + 2.0) / pSoil.fs;
	sm = (sigct[0] + sigct[1] + sigct[2]) * 0.333333333;

	if (sm < 1.0)
		sm = 1.0;

	if (pSoil.e > 1.0)
		e0 = pSoil.e;
	else
	{
		e0 = 3 * (1 - 2 * poi) * (1 + eps0) * fabs(sm) / zkapa;
	}

	hee = (pSoil.ramda - pSoil.kapa) / (1.0 + eps0);

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
	f11 = tt2 * 0.333333333 + dx * tt1 * 3.0;
	f22 = tt2 * 0.333333333 + dy * tt1 * 3.0;
	f33 = tt2 * 0.333333333 + dz * tt1 * 3.0;
	f12 = 3.0 * dxy * tt1;
	f23 = 3.0 * dyz * tt1;
	f31 = 3.0 * dzx * tt1;

	rous = qq * (pSoil.sitas - pStrePar.satu);
	roue = hee * (1.0 + pSoil.eps0) * log(pSoil.sme0 / sm);
	ttt = pSoil.ann * roue + pSoil.bnn * rous;
	ttt = pow(ttt, pSoil.beta);

	hdd = f11 + f22 + f33 + ttt / sm;

	fkk = f11 + f22 + f33;
	fij = 2.0 * (f11 * f11 + f22 * f22 + f33 * f33) + f12 * f12 + f23 * f23 + f31 * f31;
	cc7 = dlam * (fkk * fkk) + dmu * fij + hdd / hee;
	d = 1.0 / cc7;

	df = dlam * fkk;
	dg = dlam * fkk;

	ramda = (dmu2 * f11 + dlam * fkk) * pStrePar.deps[0] + (dmu2 * f22 + dlam * fkk) * pStrePar.deps[1] + (dmu2 * f33 + dlam * fkk) * pStrePar.deps[2] + dmu * (f12 * pStrePar.deps[3] + f23 * pStrePar.deps[4] + f31 * pStrePar.deps[5]);

	omake = qq * pStrePar.dsatu / (pSoil.ramda - pSoil.kapa);
	ramda = ramda + omake;

	if (ramda < 0.0)
	{

		erfs[0] = omake * d * (dept[0][0] * f11 + dept[0][1] * f22 + dept[0][2] * f33);
		erfs[1] = omake * d * (dept[1][0] * f11 + dept[1][1] * f22 + dept[1][2] * f33);
		erfs[2] = omake * d * (dept[2][0] * f11 + dept[2][1] * f22 + dept[2][2] * f33);
		erfs[3] = omake * d * (dept[3][3] * f12);
		erfs[4] = omake * d * (dept[4][4] * f23);
		erfs[5] = omake * d * (dept[5][5] * f31);

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
