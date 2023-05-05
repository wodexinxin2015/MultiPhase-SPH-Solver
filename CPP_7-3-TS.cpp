/************************************************************************************
Multi-Phase and Parallelized SPH program
--parallelized and optimized for SMP (Symmetric multiprocessing ) architecture.
--Copyright (c) 2018-2023, Weijie ZHANG, GeoHohai, Hohai University.
************************************************************************************/

#include "Class_Functions.h"
#include "Header_Parameters.h"
#include "Header_Option.h"

/*TS transformation by Yao Yangping and SMP criterion by Matsuoka-Nakai*/
void clStraStre_Fun::TSSMP_before(Particle *pPar, double alpha)
{

	double q, qc, sx[6];
	double I1, I2, I3;
	double temp;

	qc = 0.0;
	pPar->tsqc = 0.0;

	//invariants
	sx[0] = -pPar->sig[0] + pPar->vsig;
	sx[1] = -pPar->sig[1] + pPar->vsig;
	sx[2] = -pPar->sig[2] + pPar->vsig;
	sx[3] = -pPar->sig[3];
	sx[4] = -pPar->sig[4];
	sx[5] = -pPar->sig[5];
	q = pPar->divq;
	I1 = -(pPar->strep[0] + pPar->strep[1] + pPar->strep[2]);
	I2 = pPar->strep[0] * pPar->strep[1] + pPar->strep[1] * pPar->strep[2] + pPar->strep[2] * pPar->strep[1];
	I3 = -pPar->strep[0] * pPar->strep[1] * pPar->strep[2];

	//calculate qc with nonlinear criterion
	if ((I1 * I1 - 3.0 * I2) >= 1.0e-7)
	{
		temp = 3.0 * sqrt((I1 * I2 - I3) / (I1 * I2 - 9.0 * I3)) - 1.0;
		if (fabs(temp) >= 1.0e-7)
			qc = alpha * sqrt(I1 * I1 - 3.0 * I2) + 2.0 * (1.0 - alpha) * I1 / (3.0 * sqrt((I1 * I2 - I3) / (I1 * I2 - 9.0 * I3)) - 1.0);
		else
			qc = 0.0;
	}
	else
		qc = 0.0;
	pPar->tsqc = qc;

	pPar->tssig[0] = pPar->sig[0];
	pPar->tssig[1] = pPar->sig[1];
	pPar->tssig[2] = pPar->sig[2];
	pPar->tssig[3] = pPar->sig[3];
	pPar->tssig[4] = pPar->sig[4];
	pPar->tssig[5] = pPar->sig[5];

	//calculate transformed stress tensor
	if (q != 0.0)
	{
		pPar->sig[0] = -I1 * 0.33333333 - qc / q * sx[0];
		pPar->sig[1] = -I1 * 0.33333333 - qc / q * sx[1];
		pPar->sig[2] = -I1 * 0.33333333 - qc / q * sx[2];
		pPar->sig[3] = -qc / q * sx[3];
		pPar->sig[4] = -qc / q * sx[4];
		pPar->sig[5] = -qc / q * sx[5];
	}
}

void clStraStre_Fun::TSSMP_after(Particle *pPar)
{
	double q, qc, dvsig, temp[6];

	q = pPar->divq;
	qc = pPar->tsqc;

	dvsig = (pPar->dsig[0] + pPar->dsig[1] + pPar->dsig[2]) / 3.0;

	if (qc != 0.0)
	{
		pPar->dsig[0] = dvsig + q / qc * (pPar->dsig[0] - dvsig);
		pPar->dsig[1] = dvsig + q / qc * (pPar->dsig[1] - dvsig);
		pPar->dsig[2] = dvsig + q / qc * (pPar->dsig[2] - dvsig);
		pPar->dsig[3] = q / qc * pPar->dsig[3];
		pPar->dsig[4] = q / qc * pPar->dsig[4];
		pPar->dsig[5] = q / qc * pPar->dsig[5];
	}

	temp[0] = pPar->tssig[0];
	temp[1] = pPar->tssig[1];
	temp[2] = pPar->tssig[2];
	temp[3] = pPar->tssig[3];
	temp[4] = pPar->tssig[4];
	temp[5] = pPar->tssig[5];

	pPar->tssig[0] = pPar->sig[0];
	pPar->tssig[1] = pPar->sig[1];
	pPar->tssig[2] = pPar->sig[2];
	pPar->tssig[3] = pPar->sig[3];
	pPar->tssig[4] = pPar->sig[4];
	pPar->tssig[5] = pPar->sig[5];

	pPar->sig[0] = temp[0];
	pPar->sig[1] = temp[1];
	pPar->sig[2] = temp[2];
	pPar->sig[3] = temp[3];
	pPar->sig[4] = temp[4];
	pPar->sig[5] = temp[5];
}
