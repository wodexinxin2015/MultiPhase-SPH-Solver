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

/* consititutive model */
/* anisotrpic cam-clay model*/
void clStraStre_Fun::sub_camanisotropy(double (*dept)[6], double *erfs, const Particle &pStrePar, const Para_Soil &pSoil, double *adps)
{
	/*int ki;
	double e0, dlam, dmu, dmu2, sjj, rt2, rt3, zkt0, sj;
	double eta[6], sigct[6], etaval;
	double mm, sm, hee, sx, sy, sz, dx, dy, dz, dxy, dyz, dzx;
	double tt1, tt2, ttt, f11, f22, f33, f12, f23, f31, hdd, fkk, fij, cc7, d, df, dg, ramda, omake;
	double ww1, ww2, ww3, ww4, ww5, ww0, zz1, zz2, zz3, zz4, zz5, zz0, qq;

	double zkapa = pSoil.kapa;
	double poi = pSoil.poi;
	double eps0 = pSoil.eps0;

	for (ki = 0; ki < 6; ki++) {
		sigct[ki] = -pStrePar.sig[ki];
	}

	//mean stress and critical state transformation
	//mm=ita=q/p' at critical state
	mm = 3.0*(pSoil.zmf - 1.0) / (pSoil.zmf + 2.0) / pSoil.fs;
	sm = (sigct[0] + sigct[1] + sigct[2]) / 3.0;

	if (sm < 1.0) sm = 1.0;

	//elastic module
	if (pSoil.e > 1.0) e0 = pSoil.e;
	else {
		e0 = 3 * (1 - 2 * poi)*(1 + eps0) * fabs(sm) / zkapa;
	}

	//elastic stiffness matrix
	rt3 = 1.7320508100; //sqrt(2)
	rt2 = 1.4142135624; //sqrt(3)
	hee = (pSoil.ramda - pSoil.kapa) / (1.0 + eps0);


	dlam = e0 * poi / (1.0 + poi) / (1.0 - 2.0 * poi);
	dmu = e0 / 2.0 / (1.0 + poi);
	dmu2 = dmu*2.0;

	zkt0 = e0 / 3.0 / (1.0 - 2.0 * poi);

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

	//deviatoric stress tensor
	sx = sigct[0];
	sy = sigct[1];
	sz = sigct[2];
	dx = sx - sm;
	dy = sy - sm;
	dz = sz - sm;
	dxy = sigct[3];
	dyz = sigct[4];
	dzx = sigct[5];
	sjj = 0.5*(dx*dx + dy*dy + dz*dz) + dxy*dxy + dyz*dyz + dzx*dzx;
	sj = sqrt(sjj);

	//eta=sij/sm
	eta[0] = dx / sm;
	eta[1] = dy / sm;
	eta[2] = dz / sm;
	eta[3] = dxy / sm;
	eta[4] = dyz / sm;
	eta[5] = dzx / sm;

	//etaval=sqrt(3/2*eta*eta)
	etaval = 0.0;
	for (ki = 0; ki < 3; ki++) {
		etaval = etaval + eta[ki] * eta[ki] + 2.0*eta[ki + 3] * eta[ki + 3];
	}
	etaval = sqrt(etaval);*/
}
