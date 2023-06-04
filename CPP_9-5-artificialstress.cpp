/************************************************************************************
Multi-Phase and Parallelized SPH program
--parallelized and optimized for SMP (Symmetric multiprocessing ) architecture.
--Copyright (c) 2018-2023, Weijie ZHANG, GeoHohai, Hohai University.
************************************************************************************/

#include <cmath>
#include "Header_Parameters.h"
#include "Header_Option.h"
#include "Class_Functions.h"

//multiply the matrix C[][]=a[][]*b[][]
inline void matrix_multiply(double (*a)[3], double (*b)[3], double (*c)[3])
{
	c[0][0] = a[0][0] * b[0][0] + a[0][1] * b[1][0] + a[0][2] * b[2][0];
	c[0][1] = a[0][0] * b[0][1] + a[0][1] * b[1][1] + a[0][2] * b[2][1];
	c[0][2] = a[0][0] * b[0][2] + a[0][1] * b[1][2] + a[0][2] * b[2][2];

	c[1][0] = a[1][0] * b[0][0] + a[1][1] * b[1][0] + a[1][2] * b[2][0];
	c[1][1] = a[1][0] * b[0][1] + a[1][1] * b[1][1] + a[1][2] * b[2][1];
	c[1][2] = a[1][0] * b[0][2] + a[1][1] * b[1][2] + a[1][2] * b[2][2];

	c[2][0] = a[2][0] * b[0][0] + a[2][1] * b[1][0] + a[2][2] * b[2][0];
	c[2][1] = a[2][0] * b[0][1] + a[2][1] * b[1][1] + a[2][2] * b[2][1];
	c[2][2] = a[2][0] * b[0][2] + a[2][1] * b[1][2] + a[2][2] * b[2][2];
}
//artificial stress for the 2D problem
void clAcce_Fun::artificial_stress_2d(Particle *pPar, Par_Cell *pParCell, const Para_Pro &pPPro, 
	clRaij* Raij, int cn)
{
	int i, j, k, nc, nj;
	int ntotal, ndim;
	double tempup[3], hl, dst, q;
	double fijn, Rabi[3][3];
	double sita, sigmaxx, sigmayy, rxx, ryy;
	double xishu;

	xishu = pPPro.art_stre_coe;
	ntotal = pPPro.ntotal;
	ndim = pPPro.ndim;

#pragma omp parallel for schedule(static) private(k, Rabi, sita, sigmaxx, sigmayy, rxx, ryy)
	for (i = 0; i < ntotal; i++)
	{
		if (pPar[i].type == 2)
		{
			for (k = 0; k < 3; k++)
			{
				Rabi[k][0] = 0.0;
				Rabi[k][1] = 0.0;
				Rabi[k][2] = 0.0;
			}

			//solve the Rabi
			if (fabs(pPar[i].sig[0] - pPar[i].sig[1]) < 1.0e-6)
			{
				if (pPar[i].sig[3] > 0.0) sita = -pi * 0.78539816325;
				else if (pPar[i].sig[3] < 0.0) sita = pi * 0.78539816325;
				else sita = 0.0;
			}
			else
				sita = 0.5 * atan(2.0 * pPar[i].sig[3] / (pPar[i].sig[0] - pPar[i].sig[1]));

			sigmaxx = cos(sita) * cos(sita) * pPar[i].sig[0] + 2.0 * sin(sita) * cos(sita) * pPar[i].sig[3] + sin(sita) * sin(sita) * pPar[i].sig[1];
			sigmayy = sin(sita) * sin(sita) * pPar[i].sig[0] - 2.0 * sin(sita) * cos(sita) * pPar[i].sig[3] + cos(sita) * cos(sita) * pPar[i].sig[1];

			if (sigmaxx > 0.0)
				rxx = -xishu * sigmaxx / pPar[i].rho / pPar[i].rho;
			else
				rxx = 0.0;
			if (sigmayy > 0.0)
				ryy = -xishu * sigmayy / pPar[i].rho / pPar[i].rho;
			else
				ryy = 0.0;

			Rabi[0][0] = rxx * cos(sita) * cos(sita) + ryy * sin(sita) * sin(sita);
			Rabi[1][1] = rxx * sin(sita) * sin(sita) + ryy * cos(sita) * cos(sita);
			Rabi[0][1] = (rxx - ryy) * sin(sita) * cos(sita);
			Rabi[1][0] = Rabi[0][1];

			for (k = 0; k < 3; k++)
			{
				Raij[i].Rij[k][0] = Rabi[k][0];
				Raij[i].Rij[k][1] = Rabi[k][1];
				Raij[i].Rij[k][2] = Rabi[k][2];
			}
		}
	}

#pragma omp parallel for schedule(static) private(j, k, nc, nj, tempup, hl, dst, \
																q, fijn)
	for (i = 0; i < ntotal; i++)
	{
		if (pPar[i].type == 2)
		{
			tempup[0] = 0.0;
			tempup[1] = 0.0;
			tempup[2] = 0.0;

			nc = pParCell[i].ninflu;
			for (nj = 1; nj <= nc; nj++)
			{
				j = pParCell[i].influ[nj];
				if (pPar[i].type == pPar[j].type)
				{
					//distance between two particles
					dst = 0.0;
					for (k = 0; k < 3; k++)
					{
						dst = dst + (pPar[i].xp[k] - pPar[j].xp[k]) * (pPar[i].xp[k] - pPar[j].xp[k]);
					}
					dst = sqrt(dst);
					hl = (pPar[i].hl + pPar[j].hl) * 0.5;
					q = dst / hl;

					//solve fijn
					if (q >= 0.0 && q < 1.0)
					{
						fijn = (0.666666666666667 - q * q + 0.5 * q * q * q) * 3.82300884955752;
						fijn = pow(fijn, 2.55);
					}
					else if (q >= 1.0 && q < 2.0)
					{
						fijn = (2.0 - q) * (2.0 - q) * (2.0 - q) * 0.629737609329447;
						fijn = pow(fijn, 2.55);
					}
					else
						fijn = 0.0;

					//sovle the acceleration
					tempup[0] = tempup[0] + pPar[j].mass * fijn * (Raij[i].Rij[0][0] + Raij[j].Rij[0][0]) * pParCell[i].wij[nj][0]
						+ pPar[j].mass * fijn * (Raij[i].Rij[0][1] + Raij[j].Rij[0][1]) * pParCell[i].wij[nj][1];
					tempup[1] = tempup[1] + pPar[j].mass * fijn * (Raij[i].Rij[1][0] + Raij[j].Rij[1][0]) * pParCell[i].wij[nj][0]
						+ pPar[j].mass * fijn * (Raij[i].Rij[1][1] + Raij[j].Rij[1][1]) * pParCell[i].wij[nj][1];
				}
			}

			for (k = 0; k < ndim; k++)
				pPar[i].ax[k] = pPar[i].ax[k] + tempup[k];
		}
	}
}

//calculating the artificial stress for the 3D problem
void clAcce_Fun::artificial_stress_3d(Particle *pPar, Par_Cell *pParCell, const Para_Pro &pPPro, 
	clRaij* Raij, int cn)
{
	int i, j, k, nc, nj;
	int ntotal, ndim;
	double tempup[3], hl, dst, q;
	double fijn, Rabi[3][3], stre[3][3];
	double beta_in[3][3], beta_ints[3][3], temp[3][3];
	double sigmaxx, sigmayy, sigmazz, rxx, ryy, rzz;
	double xishu;

	//initialiing of variables
	xishu = pPPro.art_stre_coe;
	ntotal = pPPro.ntotal;
	ndim = pPPro.ndim;

#pragma omp parallel for schedule(static) private(k, Rabi,sigmaxx, sigmayy, sigmazz, \
rxx, ryy, rzz, beta_in, beta_ints, temp)
	for (i = 0; i < ntotal; i++)
	{
		if (pPar[i].type == 2)
		{
			for (k = 0; k < 3; k++)
			{
				Rabi[k][0] = 0.0;
				Rabi[k][1] = 0.0;
				Rabi[k][2] = 0.0;

				temp[k][0] = 0.0;
				temp[k][1] = 0.0;
				temp[k][2] = 0.0;

				beta_in[k][0] = 0.0;
				beta_in[k][1] = 0.0;
				beta_in[k][2] = 0.0;

				beta_ints[k][0] = 0.0;
				beta_ints[k][1] = 0.0;
				beta_ints[k][2] = 0.0;
			}

			//sigmaxx=sigma1; sigmaxx=sigma2; sigmaxx=sigma3
			sigmaxx = pPar[i].strep[0];
			sigmayy = pPar[i].strep[1];
			sigmazz = pPar[i].strep[2];

			//checking if the principle stress is positive
			if (sigmaxx > 0.0)
				rxx = -xishu * sigmaxx / pPar[i].rho / pPar[i].rho;
			else
				rxx = 0.0;
			if (sigmayy > 0.0)
				ryy = -xishu * sigmayy / pPar[i].rho / pPar[i].rho;
			else
				ryy = 0.0;
			if (sigmazz > 0.0)
				rzz = -xishu * sigmazz / pPar[i].rho / pPar[i].rho;
			else
				rzz = 0.0;

			//transforming the principle stress to ordinary space
			//inverse of beta
			inverse_mat(pPar[i].beta, beta_in);

			//inverse of beta_ts
			inverse_mat(pPar[i].beta_ts, beta_ints);

			//forming matrix stre
			for (k = 0; k < 3; k++)
			{
				stre[k][0] = 0.0;
				stre[k][1] = 0.0;
				stre[k][2] = 0.0;
			}
			stre[0][0] = rxx;
			stre[1][1] = ryy;
			stre[2][2] = rzz;

			//[temp]=inverse(beta)*[stre]
			matrix_multiply(beta_in, stre, temp);

			//[Rabi]=[temp]*inverse(beta_ts)
			matrix_multiply(temp, beta_ints, Rabi);

			for (k = 0; k < 3; k++)
			{
				Raij[i].Rij[k][0] = Rabi[k][0];
				Raij[i].Rij[k][1] = Rabi[k][1];
				Raij[i].Rij[k][2] = Rabi[k][2];
			}
		}
	}

#pragma omp parallel for schedule(static) private(j, k, nc, nj, tempup, hl, dst, \
																q, fijn)
	for (i = 0; i < ntotal; i++)
	{
		if (pPar[i].type == 2)
		{
			tempup[0] = 0.0;
			tempup[1] = 0.0;
			tempup[2] = 0.0;

			nc = pParCell[i].ninflu;
			for (nj = 1; nj <= nc; nj++)
			{
				j = pParCell[i].influ[nj];
				if (pPar[i].type == pPar[j].type)
				{
					//distance between two particles
					dst = 0.0;
					for (k = 0; k < 3; k++)
					{
						dst = dst + (pPar[i].xp[k] - pPar[j].xp[k]) * (pPar[i].xp[k] - pPar[j].xp[k]);
					}
					dst = sqrt(dst);
					hl = (pPar[i].hl + pPar[j].hl) / 2;
					q = dst / hl;

					//solve fijn
					if (q >= 0.0 && q < 1.0)
					{
						fijn = (0.666666666666667 - q * q + 0.5 * q * q * q) * 3.82300884955752;
						fijn = pow(fijn, 2.55);
					}
					else if (q >= 1.0 && q < 2.0)
					{
						fijn = (2.0 - q) * (2.0 - q) * (2.0 - q) * 0.629737609329447;
						fijn = pow(fijn, 2.55);
					}
					else
						fijn = 0.0;

					//sovle the acceleration
					tempup[0] = tempup[0] + pPar[j].mass * fijn * (Raij[i].Rij[0][0] + Raij[j].Rij[0][0]) * pParCell[i].wij[nj][0]
						+ pPar[j].mass * fijn * (Raij[i].Rij[0][1] + Raij[j].Rij[0][1]) * pParCell[i].wij[nj][1]
						+ pPar[j].mass * fijn * (Raij[i].Rij[0][2] + Raij[j].Rij[0][2]) * pParCell[i].wij[nj][2];
					tempup[1] = tempup[1] + pPar[j].mass * fijn * (Raij[i].Rij[1][0] + Raij[j].Rij[1][0]) * pParCell[i].wij[nj][0]
						+ pPar[j].mass * fijn * (Raij[i].Rij[1][1] + Raij[j].Rij[1][1]) * pParCell[i].wij[nj][1]
						+ pPar[j].mass * fijn * (Raij[i].Rij[1][2] + Raij[j].Rij[1][2]) * pParCell[i].wij[nj][2];
					tempup[2] = tempup[2] + pPar[j].mass * fijn * (Raij[i].Rij[2][0] + Raij[j].Rij[2][0]) * pParCell[i].wij[nj][0]
						+ pPar[j].mass * fijn * (Raij[i].Rij[2][1] + Raij[j].Rij[2][1]) * pParCell[i].wij[nj][1]
						+ pPar[j].mass * fijn * (Raij[i].Rij[2][2] + Raij[j].Rij[2][2]) * pParCell[i].wij[nj][2];
				}
			}

			for (k = 0; k < ndim; k++)
				pPar[i].ax[k] = pPar[i].ax[k] + tempup[k];
		}
	}
}
//get the inverse matrix of a[3][3]
void clAcce_Fun::inverse_mat(double (*A)[3], double (*A_inv)[3])
{
	double det;

	det = A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) - A[1][0] * (A[0][1] * A[2][2] - A[0][2] * A[2][1]) + A[2][0] * (A[0][1] * A[1][2] - A[0][2] * A[1][1]);

	if (fabs(det) > 1.0e-5)
	{
		A_inv[0][0] = (A[1][1] * A[2][2] - A[1][2] * A[2][1]) / det;
		A_inv[0][1] = (A[0][2] * A[2][1] - A[0][1] * A[2][2]) / det;
		A_inv[0][2] = (A[0][1] * A[1][2] - A[0][2] * A[1][1]) / det;

		A_inv[1][0] = (A[1][2] * A[2][0] - A[1][0] * A[2][2]) / det;
		A_inv[1][1] = (A[0][0] * A[2][2] - A[0][2] * A[2][0]) / det;
		A_inv[1][2] = (A[1][0] * A[0][2] - A[0][0] * A[1][2]) / det;

		A_inv[2][0] = (A[1][0] * A[2][1] - A[1][1] * A[2][0]) / det;
		A_inv[2][1] = (A[0][1] * A[2][0] - A[0][0] * A[2][1]) / det;
		A_inv[2][2] = (A[0][0] * A[1][1] - A[1][0] * A[0][1]) / det;
	}
}
