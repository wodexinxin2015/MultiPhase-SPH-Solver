/************************************************************************************
Multi-Phase and Parallelized SPH program
--parallelized and optimized for SMP (Symmetric multiprocessing ) architecture.
--Copyright (c) 2018-2023, Weijie ZHANG, GeoHohai, Hohai University.
************************************************************************************/
#include <math.h>
#include <map>
#include <utility>
#include <cstring>
#include "Class_Functions.h"
#include "Header_Option.h"

double moisturesuct(double satu, double dsatu, const Para_Soil &Soil);
double moisturesr(double satu, double dsatu, const Para_Soil &Soil);
int eejcb(double a[], int n, double v[], double eps, int jt);

double moisturesuct(double satu, double dsatu, const Para_Soil &Soil);
double moisturesr(double satu, double dsatu, const Para_Soil &Soil);

void clStraStre_Fun::strain_noreg(Particle *pPar, Par_Cell *pParCell, 
	cl_bndy_pair* pBndy_Pair, clVar_Boundary* pParti_VariBndy, 
	const Para_Pro &pPPro, int bndytype)
{
	int i, j, k, nc, nj;
	int ntotal, dim;
	double dvx[3], tempup[9], vxj[3], tempdown[3];

	ntotal = pPPro.ntotal;
	dim = pPPro.ndim;

#pragma omp parallel for schedule(static) private(j, k, nc, nj, dvx, tempup, tempdown, vxj)
	for (i = 0; i < ntotal; i++)
	{

		if (pPPro.l == pPPro.inip)
		{
			for (j = 0; j < 6; j++)
				pPar[i].eps[j] = 0.0;
		}

		for (j = 0; j < 6; j++)
		{
			pPar[i].deps[j] = 0.0;
			pPar[i].dsig[j] = 0.0;
		}
		for (j = 0; j < 3; j++)
		{
			pPar[i].weps[j][0] = 0.0;
			pPar[i].weps[j][1] = 0.0;
			pPar[i].weps[j][2] = 0.0;
			tempup[j] = 0.0;
			tempup[j + 3] = 0.0;
			tempup[j + 6] = 0.0;

			pPar[i].divx[j][0] = 0.0;
			pPar[i].divx[j][1] = 0.0;
			pPar[i].divx[j][2] = 0.0;

			tempdown[j] = 0;
		}

		nc = pParCell[i].ninflu;
		for (nj = 1; nj <= nc; nj++)
		{
			for (k = 0; k < 3; k++)
				dvx[k] = 0.0;
			j = pParCell[i].influ[nj];

			if (pPar[i].type == pPar[j].type && pPar[i].type != 7 && pPar[i].type != 0)
			{
				for (k = 0; k < 3; k++)
					dvx[k] = pPar[i].vxp[k] - pPar[j].vxp[k];

				tempup[0] = tempup[0] + pPar[j].mass / pPar[j].rho * dvx[0] * pParCell[i].wij[nj][0];
				tempup[1] = tempup[1] + pPar[j].mass / pPar[j].rho * dvx[0] * pParCell[i].wij[nj][1];
				tempup[2] = tempup[2] + pPar[j].mass / pPar[j].rho * dvx[0] * pParCell[i].wij[nj][2];

				tempup[3] = tempup[3] + pPar[j].mass / pPar[j].rho * dvx[1] * pParCell[i].wij[nj][0];
				tempup[4] = tempup[4] + pPar[j].mass / pPar[j].rho * dvx[1] * pParCell[i].wij[nj][1];
				tempup[5] = tempup[5] + pPar[j].mass / pPar[j].rho * dvx[1] * pParCell[i].wij[nj][2];

				tempup[6] = tempup[6] + pPar[j].mass / pPar[j].rho * dvx[2] * pParCell[i].wij[nj][0];
				tempup[7] = tempup[7] + pPar[j].mass / pPar[j].rho * dvx[2] * pParCell[i].wij[nj][1];
				tempup[8] = tempup[8] + pPar[j].mass / pPar[j].rho * dvx[2] * pParCell[i].wij[nj][2];

				if (pPar[i].type == 2 || pPar[i].type == 4) {
					tempdown[0] = tempdown[0] + pPar[j].mass / pPar[j].rho * (pPar[i].xp[0] - pPar[j].xp[0]) * pParCell[i].wij[nj][0];
					tempdown[1] = tempdown[1] + pPar[j].mass / pPar[j].rho * (pPar[i].xp[1] - pPar[j].xp[1]) * pParCell[i].wij[nj][1];
					tempdown[2] = tempdown[2] + pPar[j].mass / pPar[j].rho * (pPar[i].xp[2] - pPar[j].xp[2]) * pParCell[i].wij[nj][2];
				}
			}
			else if ((pPar[i].type == 1 || pPar[i].type == 2 || pPar[i].type == 3 || pPar[i].type == 4) && pPar[j].type == 0) {
				//settting virtual velocity
				vxj[0] = 0.0;
				vxj[1] = 0.0;
				vxj[2] = 0.0;
				if (bndytype == 1) {
					vxj[0] = pPar[i].vxp[0];
					vxj[1] = pPar[i].vxp[1];
					vxj[2] = pPar[i].vxp[2];
				}
				else if(bndytype == 2 || bndytype == 3 || bndytype == 4 || bndytype == 13 || bndytype == 14){
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

				for (k = 0; k < 3; k++)
					dvx[k] = pPar[i].vxp[k] - vxj[k];

				tempup[0] = tempup[0] + pPar[j].mass / pPar[j].rho * dvx[0] * pParCell[i].wij[nj][0];
				tempup[1] = tempup[1] + pPar[j].mass / pPar[j].rho * dvx[0] * pParCell[i].wij[nj][1];
				tempup[2] = tempup[2] + pPar[j].mass / pPar[j].rho * dvx[0] * pParCell[i].wij[nj][2];

				tempup[3] = tempup[3] + pPar[j].mass / pPar[j].rho * dvx[1] * pParCell[i].wij[nj][0];
				tempup[4] = tempup[4] + pPar[j].mass / pPar[j].rho * dvx[1] * pParCell[i].wij[nj][1];
				tempup[5] = tempup[5] + pPar[j].mass / pPar[j].rho * dvx[1] * pParCell[i].wij[nj][2];

				tempup[6] = tempup[6] + pPar[j].mass / pPar[j].rho * dvx[2] * pParCell[i].wij[nj][0];
				tempup[7] = tempup[7] + pPar[j].mass / pPar[j].rho * dvx[2] * pParCell[i].wij[nj][1];
				tempup[8] = tempup[8] + pPar[j].mass / pPar[j].rho * dvx[2] * pParCell[i].wij[nj][2];

				if (pPar[i].type == 2 || pPar[i].type == 4) {
					tempdown[0] = tempdown[0] + pPar[j].mass / pPar[j].rho * (pPar[i].xp[0] - pPar[j].xp[0]) * pParCell[i].wij[nj][0];
					tempdown[1] = tempdown[1] + pPar[j].mass / pPar[j].rho * (pPar[i].xp[1] - pPar[j].xp[1]) * pParCell[i].wij[nj][1];
					tempdown[2] = tempdown[2] + pPar[j].mass / pPar[j].rho * (pPar[i].xp[2] - pPar[j].xp[2]) * pParCell[i].wij[nj][2];
				}
			}
		}

		if (pPar[i].type == 1 || pPar[i].type == 3)
		{
			pPar[i].divx[0][0] = tempup[0];
			pPar[i].divx[1][0] = tempup[3];
			pPar[i].divx[2][0] = tempup[6];

			pPar[i].divx[0][1] = tempup[1];
			pPar[i].divx[1][1] = tempup[4];
			pPar[i].divx[2][1] = tempup[7];

			pPar[i].divx[0][2] = tempup[2];
			pPar[i].divx[1][2] = tempup[5];
			pPar[i].divx[2][2] = tempup[8];
		}
		else
		{
			if (fabs(tempdown[0]) > 1.0e-7) {
				pPar[i].divx[0][0] = tempup[0] / tempdown[0];
				pPar[i].divx[1][0] = tempup[3] / tempdown[0];
				pPar[i].divx[2][0] = tempup[6] / tempdown[0];
			}

			if (fabs(tempdown[1]) > 1.0e-7) {
				pPar[i].divx[0][1] = tempup[1] / tempdown[1];
				pPar[i].divx[1][1] = tempup[4] / tempdown[1];
				pPar[i].divx[2][1] = tempup[7] / tempdown[1];
			}

			if (fabs(tempdown[2]) > 1.0e-7) {
				pPar[i].divx[0][2] = tempup[2] / tempdown[2];
				pPar[i].divx[1][2] = tempup[5] / tempdown[2];
				pPar[i].divx[2][2] = tempup[8] / tempdown[2];
			}
		}

		//strain rate tensor
		pPar[i].deps[0] = pPar[i].divx[0][0];
		pPar[i].deps[1] = pPar[i].divx[1][1];
		pPar[i].deps[2] = pPar[i].divx[2][2];
		pPar[i].deps[3] = pPar[i].divx[0][1] + pPar[i].divx[1][0];
		pPar[i].deps[4] = pPar[i].divx[1][2] + pPar[i].divx[2][1];
		pPar[i].deps[5] = pPar[i].divx[2][0] + pPar[i].divx[0][2];

		if (dim == 2)
		{
			pPar[i].deps[2] = 0.0;
			pPar[i].deps[4] = 0.0;
			pPar[i].deps[5] = 0.0;
		}

		//spin rate tensor
		pPar[i].weps[0][0] = 0.0;
		pPar[i].weps[0][1] = 0.5 * (pPar[i].divx[0][1] - pPar[i].divx[1][0]);
		pPar[i].weps[0][2] = 0.5 * (pPar[i].divx[0][2] - pPar[i].divx[2][0]);

		pPar[i].weps[1][0] = 0.5 * (pPar[i].divx[1][0] - pPar[i].divx[0][1]);
		pPar[i].weps[1][1] = 0.0;
		pPar[i].weps[1][2] = 0.5 * (pPar[i].divx[1][2] - pPar[i].divx[2][1]);

		pPar[i].weps[2][0] = 0.5 * (pPar[i].divx[2][0] - pPar[i].divx[0][2]);
		pPar[i].weps[2][1] = 0.5 * (pPar[i].divx[2][1] - pPar[i].divx[1][2]);
		pPar[i].weps[2][2] = 0.0;

		if (dim == 2)
		{
			pPar[i].weps[0][2] = 0.0;
			pPar[i].weps[1][2] = 0.0;
			pPar[i].weps[2][0] = 0.0;
			pPar[i].weps[2][1] = 0.0;
			pPar[i].weps[2][2] = 0.0;
		}
	}
}

void clInterFor_Fun::saturation(Particle *pPar, Par_Cell *pParCell, const Para_Soil *pParti_ConsPara, const Para_Pro &pPPro, int cn)
{

	int i, j, nc, nj;
	int ntotal = pPPro.ntotal;
	double cof, psatu, tempup, tempdown, sitas, sitar;

#pragma omp parallel for schedule(static) private(j, nc, nj, psatu, tempup, tempdown, cof, sitas, sitar)
	for (i = 0; i < ntotal; i++)
	{
		cof = pParti_ConsPara[i].sitas - pParti_ConsPara[i].sitar;
		sitas = pParti_ConsPara[i].sitas;
		sitar = pParti_ConsPara[i].sitar;

		tempup = 0.0;
		tempdown = 0.0;

		nc = pParCell[i].ninflu;
		for (nj = 1; nj <= nc; nj++)
		{
			j = pParCell[i].influ[nj];
			if (pPar[i].type != pPar[j].type && (pPar[i].type == 2 && pPar[j].type == 1))
			{
				tempup = tempup + pPar[j].mass / pPar[j].rho * cof * pParCell[i].wij[nj][3];
				tempdown = tempdown + pPar[j].mass / pPar[j].rho * pParCell[i].wij[nj][3];
			}
		}
		if (pPar[i].type == 2)
		{
			psatu = pPar[i].satu;
			if (tempdown != 0.0)
				pPar[i].satu = tempup / tempdown + sitar;
			else
				pPar[i].satu = tempup + sitar;

			pPar[i].dsatu = pPar[i].satu - psatu;

			if (pPar[i].dsatu < 0)
			{ //wetting
				pPar[i].satu = pPar[i].satu - pPar[i].dsatu;
				pPar[i].dsatu = 0.0;
			}

			/*if (pPar[i].dsatu > 0){ //drying
			pPar[i].satu = pPar[i].satu - pPar[i].dsatu;
			pPar[i].dsatu = 0.0;
			}*/

			//correct the saturation
			if (pPar[i].satu > sitas)
				pPar[i].satu = sitas;
			else if (pPar[i].satu < sitar)
				pPar[i].satu = sitar;

			/*Suction*/
			pPar[i].suct = moisturesuct(pPar[i].satu, pPar[i].dsatu, pParti_ConsPara[i]);
		}
	}
}

//calculate the principle stresses
void clStraStre_Fun::Stress_eigValue(Particle *pPar, const Para_Pro &pPPro, int cn)
{
	int i, k, err;
	double stre[9], eig_vector[9];
	int ntotal = pPPro.ntotal;

#pragma omp parallel for schedule(static) private(k, stre, eig_vector, err)
	for (i = 0; i < ntotal; i++)
	{
		//initializing of vector
		for (k = 0; k < 9; k++)
		{
			eig_vector[k] = 0.0;
		}

		//0-xx; 1-yy; 2-zz; 3-xy; 4-yz; 5-zx.
		stre[0] = pPar[i].sig[0]; //xx
		stre[1] = pPar[i].sig[3]; //xy
		stre[2] = pPar[i].sig[5]; //zx

		stre[3] = pPar[i].sig[3]; //xy
		stre[4] = pPar[i].sig[1]; //yy
		stre[5] = pPar[i].sig[4]; //yz

		stre[6] = pPar[i].sig[5]; //zx
		stre[7] = pPar[i].sig[4]; //yz
		stre[8] = pPar[i].sig[2]; //zz

		err = eejcb(stre, 3, eig_vector, 0.001, 7);

		//eigen values
		for (k = 0; k < 3; k++)
		{
			int index = k * 3 + k;
			pPar[i].strep[k] = stre[index];
		}

		//transforming matrix and its transposed matrix
		for (k = 0; k < 3; k++)
		{
			pPar[i].beta_ts[k][0] = eig_vector[k * 3];
			pPar[i].beta_ts[k][1] = eig_vector[k * 3 + 1];
			pPar[i].beta_ts[k][2] = eig_vector[k * 3 + 2];

			pPar[i].beta[0][k] = eig_vector[k * 3];
			pPar[i].beta[1][k] = eig_vector[k * 3 + 1];
			pPar[i].beta[2][k] = eig_vector[k * 3 + 2];
		}
	}
}

//Jacobi method for solving the eigen values and its corresponding vectors
//return > 0 means normal; else means not meeting the demanded precision
//a-array of n*n, input or output the eigen values
//n-dimension of matrix
//v-array of n*n, output the vectors
//eps-precision
//jt-iteration loops
int eejcb(double a[], int n, double v[], double eps, int jt)
{
	int i, j, p, q, u, w, t, s, l;
	double fm, cn, sn, omega, x, y, d;

	p = 0;
	q = 0;

	l = 1;
	for (i = 0; i <= n - 1; i++)
	{
		v[i * n + i] = 1.0;
		for (j = 0; j <= n - 1; j++)
		{
			if (i != j)
			{
				v[i * n + j] = 0.0;
			}
		}
	}
	while (1 == 1)
	{
		fm = 0.0;
		for (i = 0; i <= n - 1; i++)
		{
			for (j = 0; j <= n - 1; j++)
			{
				d = fabs(a[i * n + j]);
				if ((i != j) && (d > fm))
				{
					fm = d;
					p = i;
					q = j;
				}
			}
		}
		if (fm < eps)
		{
			return (1);
		}
		if (l > jt)
		{
			return (-1);
		}
		l = l + 1;
		u = p * n + q;
		w = p * n + p;
		t = q * n + p;
		s = q * n + q;
		x = -a[u];
		y = (a[s] - a[w]) / 2.0;
		omega = x / sqrt(x * x + y * y);
		if (y < 0.0)
		{
			omega = -omega;
		}
		sn = 1.0 + sqrt(1.0 - omega * omega);
		sn = omega / sqrt(2.0 * sn);
		cn = sqrt(1.0 - sn * sn);
		fm = a[w];
		a[w] = fm * cn * cn + a[s] * sn * sn + a[u] * omega;
		a[s] = fm * sn * sn + a[s] * cn * cn - a[u] * omega;
		a[u] = 0.0;
		a[t] = 0.0;
		for (j = 0; j <= n - 1; j++)
		{
			if ((j != p) && (j != q))
			{
				u = p * n + j;
				w = q * n + j;
				fm = a[u];
				a[u] = fm * cn + a[w] * sn;
				a[w] = -fm * sn + a[w] * cn;
			}
		}
		for (i = 0; i <= n - 1; i++)
		{
			if ((i != p) && (i != q))
			{
				u = i * n + p;
				w = i * n + q;
				fm = a[u];
				a[u] = fm * cn + a[w] * sn;
				a[w] = -fm * sn + a[w] * cn;
			}
		}
		for (i = 0; i <= n - 1; i++)
		{
			u = i * n + p;
			w = i * n + q;
			fm = v[u];
			v[u] = fm * cn + v[w] * sn;
			v[w] = -fm * sn + v[w] * cn;
		}
	}
	return (1);
}
