/************************************************************************************
Multi-Phase and Parallelized SPH program
--parallelized and optimized for SMP (Symmetric multiprocessing ) architecture.
--Copyright (c) 2018-2023, Weijie ZHANG, GeoHohai, Hohai University.
************************************************************************************/

#define _CRT_SECURE_NO_DEPRECATE
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "Class_Functions.h"
#include "Header_Parameters.h"
#include "Header_Option.h"

/* Smoothing kernel function using 3rd-order B-spline function*/
inline void weight(double dst, double ra, double dx[], double *w1ij, int ndim)
{
	double a1, s;
	int k;
	s = dst / ra;
	a1 = 0.0;
	if (ndim == 2)
	{
		a1 = 0.682092613250980 / ra / ra;
		w1ij[ndim] = 0.0;
	}
	if (ndim == 3)
		a1 = 0.477464829275686 / ra / ra / ra;
	if (s == 0.0)
	{
		for (k = 0; k < ndim; k++)
			w1ij[k] = 0.0;
		w1ij[3] = a1 * (0.66666667 - s * s + 0.5 * s * s * s);
	}
	if (s > 0.0 && s <= 1.0)
	{
		for (k = 0; k < ndim; k++)
			w1ij[k] = a1 * dx[k] * (-2 * s + 1.5 * s * s) / (ra * dst);
		w1ij[3] = a1 * (0.66666667 - s * s + 0.5 * s * s * s);
	}
	else if (s > 1.0 && s <= 2.0)
	{
		for (k = 0; k < ndim; k++)
			w1ij[k] = a1 * (-0.5) * dx[k] * (2.0 - s) * (2.0 - s) / (ra * dst);
		w1ij[3] = a1 * (2.0 - s) * (2.0 - s) * (2.0 - s) / 6.0;
	}
}
/* Smoothing kernel function by Yang and Liu 2012*/
// 2016-7-25
inline void weight2(double dst, double ra, double dx[], double *w1ij, int ndim)
{
	double a1, s;
	int k;
	s = dst / ra;
	a1 = 0.0;
	if (ndim == 2)
	{
		a1 = 0.106103295394597 / (ra * ra);
		w1ij[ndim] = 0.0;
	}
	if (ndim == 3)
		a1 = 0.077010456334788 / (ra * ra * ra);
	if (s == 0.0)
	{
		for (k = 0; k < ndim; k++)
			w1ij[k] = 0.0;
		w1ij[3] = a1 * (6.0 - 6.0 * s + 1.0 * s * s * s);
	}
	if (s > 0.0 && s <= 1.0)
	{
		for (k = 0; k < ndim; k++)
			w1ij[k] = a1 * dx[k] * (-6.0 + 3.0 * s * s) / (ra * dst);
		w1ij[3] = a1 * (6.0 - 6.0 * s + 1.0 * s * s * s);
	}
	else if (s > 1.0 && s <= 2.0)
	{
		for (k = 0; k < ndim; k++)
			w1ij[k] = a1 * (-3.0) * dx[k] * (2.0 - s) * (2.0 - s) / (ra * dst);
		w1ij[3] = a1 * (2.0 - s) * (2.0 - s) * (2.0 - s);
	}
}

// multiply the matrix c[]=a[][]*b[]
inline void mat_multiply(double (*a)[3], double(*b), double(*c))
{
	c[0] = a[0][0] * b[0] + a[0][1] * b[1] + a[0][2] * b[2];
	c[1] = a[1][0] * b[0] + a[1][1] * b[1] + a[1][2] * b[2];
	c[2] = a[2][0] * b[0] + a[2][1] * b[1] + a[2][2] * b[2];
}

void weight(double dst, double ra, double dx[], double *w1ij, int ndim);
void weight2(double dst, double ra, double dx[], double *w1ij, int ndim);

/*initializing for cell number and problem domain*/
void clNNPS_Fun::celllist_ini1(Cell_Con *pCellc, const Para_Pro &pPPro)
{
	double xm, ym, zm;
	double dr = 2.0 * pPPro.dr;

	if (pPPro.ndim == 2)
	{
		xm = (pCellc->xmax - pCellc->xmin) / dr + 3.0;
		ym = (pCellc->ymax - pCellc->ymin) / dr + 3.0;
		zm = 1.0;
		pCellc->cxm = (int)(xm);
		pCellc->cym = (int)(ym);
		pCellc->czm = (int)(zm);
	}
	else if (pPPro.ndim == 3)
	{
		xm = (pCellc->xmax - pCellc->xmin) / dr + 3.0;
		ym = (pCellc->ymax - pCellc->ymin) / dr + 3.0;
		zm = (pCellc->zmax - pCellc->zmin) / dr + 3.0;
		pCellc->cxm = (int)(xm);
		pCellc->cym = (int)(ym);
		pCellc->czm = (int)(zm);
	}

	pCellc->ctotal = pCellc->cxm * pCellc->cym * pCellc->czm;
}

/* initializing of the cell link */
void clNNPS_Fun::celllist_ini2(cell_link *pCell_Link, const Cell_Con &pCellc, int dim, int cn)
{
	int p, i, kx, ky, kz;
	int nXm = pCellc.cxm;
	int nYm = pCellc.cym;
	int nZm = pCellc.czm;
	int m = pCellc.ctotal;

	if (dim == 2)
	{
#pragma omp parallel for schedule(static) private(kx, ky, p)
		for (i = 0; i <= m; i++)
		{
			p = 0;
			for (ky = 1; ky <= 3; ky++)
			{
				for (kx = 1; kx <= 3; kx++)
				{
					pCell_Link[i].nncell[p] = (kx - 2) * nYm + ky - 2 + i;
					if (pCell_Link[i].nncell[p] < 0 || pCell_Link[i].nncell[p] > m)
						pCell_Link[i].nncell[p] = 0;
					p = p + 1;
				}
			}
		}
	}
	else if (dim == 3)
	{
#pragma omp parallel for schedule(static) private(kx, ky, kz, p)
		for (i = 0; i <= m; i++)
		{
			p = 0;
			for (kz = 1; kz <= 3; kz++)
			{
				for (kx = 1; kx <= 3; kx++)
				{
					for (ky = 1; ky <= 3; ky++)
					{
						pCell_Link[i].nncell[p] = (ky - 2) * nXm * nZm + (kx - 2) * nZm + kz - 2 + i;
						if (pCell_Link[i].nncell[p] < 0 || pCell_Link[i].nncell[p] > m)
							pCell_Link[i].nncell[p] = 0;
						p = p + 1;
					}
				}
			}
		}
	}
}

/* particles searching method*/
// fully parallelized
int clNNPS_Fun::parttocell(Par_Cell *pParCell, cell_info *pCell_Info, const Cell_Con &pCellc,
						   parti_cellid *pParti_Cell_Sorted, parti_cellid *pParti_Cell_Temp, 
						   const Para_Pro &pPPro, Particle *pPar, int dim, int cn, int size_m, int k)
{
	int nx_t[3], m;
	int i, j, ntotal;
	double nx[3], xmin[3], dr;
	int ntd, id, nXm, nYm, nZm;
	int t = 0;

	m = pCellc.ctotal;

	xmin[0] = pCellc.xmin;
	xmin[1] = pCellc.ymin;
	xmin[2] = pCellc.zmin;
	nXm = pCellc.cxm;
	nYm = pCellc.cym;
	nZm = pCellc.czm;
	dr = 2.0 * pPPro.dr;
	ntotal = pPPro.ntotal;

#pragma omp parallel for schedule(static)
	for (i = 0; i < m; i += 1)
	{
		pCell_Info[i].start = 10000000;
		pCell_Info[i].end = -10000000;
	}

	if (dim == 2)
	{
#pragma omp parallel for schedule(static) private(id, i, j, nx, nx_t)
		/* get the cell number for each particle and put particles into cell for velocity particles */
		for (ntd = 0; ntd < cn; ntd++)
		{
			id = ntd;
			for (i = id; i < ntotal; i += cn)
			{
				if (pPar[i].type != 7)
				{
					for (j = 0; j < 2; j++)
					{
						nx[j] = fabs(pPar[i].xp[j] - xmin[j]) / dr + 2.0;
						nx_t[j] = (int)(nx[j]);
						if (nx_t[j] < 2)
							nx_t[j] = 2;
					}
					pParti_Cell_Sorted[i].cell_id = nYm * (nx_t[0] - 1) + nx_t[1];
					pParti_Cell_Sorted[i].parti_id = i;
				}
			}
		}
	}
	else if (dim == 3)
	{
#pragma omp parallel for schedule(static) private(id, i, j, nx, nx_t)
		/* get the cell number for each particle and put particles into cell for velocity particles */
		for (ntd = 0; ntd < cn; ntd++)
		{
			id = ntd;
			for (i = id; i < ntotal; i += cn)
			{
				if (pPar[i].type != 7)
				{
					for (j = 0; j < 3; j++)
					{
						nx[j] = fabs(pPar[i].xp[j] - xmin[j]) / dr + 2.0;
						nx_t[j] = (int)(nx[j]);
						if (nx_t[j] < 2)
							nx_t[j] = 2;
					}
					pParti_Cell_Sorted[i].cell_id = nXm * nZm * (nx_t[1] - 1) + nZm * (nx_t[0] - 1) + nx_t[2];
					pParti_Cell_Sorted[i].parti_id = i;
				}
			}
		}
	}

	// cell id information to the pParcell
#pragma omp parallel for schedule(static)
	for (i = 0; i < ntotal; i++)
	{
		if (pParti_Cell_Sorted[i].cell_id > 0 && pParti_Cell_Sorted[i].cell_id <= m)
			pParCell[i].cell_id = pParti_Cell_Sorted[i].cell_id;
		else {
			t = 3;
		}	
	}

	if (t > 0) return t;

	//put particle id into cell information by the binary sorting method
	Binary_Sort(pParti_Cell_Sorted, pParti_Cell_Temp, pCell_Info, ntotal, m, size_m, k);

	return t;
}

int clNNPS_Fun::partisearching(Par_Cell *pParCell, cell_info *pCell_Info, cell_link *pCell_Link, parti_cellid * pParti_Cell_Sorted,
							   const Cell_Con &pCellc, const Para_Pro &pPPro, Particle *pPar, int cn, int flknl)
{
	int p, temp, t;
	int i, j, k, ntotal, nc, nj, nl, np, start, end, dim;
	double dst, dx[3], ra;

	dim = pPPro.ndim;
	ntotal = pPPro.ntotal;
	temp = (dim == 2) ? 9 : 27;
	t = 0;

#pragma omp parallel for schedule(static)
	for (i = 0; i < ntotal; i++) {
		pParCell[i].initial();
	}

	if (flknl == 1 || flknl == 3)
	{
		/*searching particles and calculation of kernel function*/
#pragma omp parallel for schedule(static) private(j, k, nc, nl, nj, np, start, end, dst, dx, ra, p)
		for (i = 0; i < ntotal; i++)
		{
			if (pPar[i].type != 7)
			{
				nc = pParCell[i].cell_id;
				for (nj = 0; nj < temp; nj++)
				{
					nl = pCell_Link[nc].nncell[nj];
					if (nl > 0)
					{
						start = pCell_Info[nl].start;
						end = pCell_Info[nl].end;
						for (np = start; np <= end; np++)
						{
							j = pParti_Cell_Sorted[np].parti_id;
							dst = 0.0;
							for (k = 0; k < dim; k++)
							{
								dst = dst + (pPar[i].xp[k] - pPar[j].xp[k]) * (pPar[i].xp[k] - pPar[j].xp[k]);
								dx[k] = pPar[i].xp[k] - pPar[j].xp[k];
							}
							dst = sqrt(dst);
							ra = (pPar[i].hl + pPar[j].hl) / 2;
							if (dst <= 2 * ra && i != j)
							{
								pParCell[i].ninflu += 1;
								p = pParCell[i].ninflu;
								if (p <= 99)
								{
									pParCell[i].influ[p] = j;
									weight(dst, ra, dx, pParCell[i].wij[p], dim);
								}
								else
									t = 5;
							}
						}
					}
				}
			}
		}
	}
	else
	{
		/*searching particles and calculation of kernel function*/
#pragma omp parallel for schedule(static) private(j, k, nc, nl, nj, np, start, end, dst, dx, ra, p)
		for (i = 0; i < ntotal; i++)
		{
			if (pPar[i].type != 7)
			{
				nc = pParCell[i].cell_id;
				for (nj = 0; nj < temp; nj++)
				{
					nl = pCell_Link[nc].nncell[nj];
					if (nl > 0)
					{
						start = pCell_Info[nl].start;
						end = pCell_Info[nl].end;
						for (np = start; np <= end; np++)
						{
							j = pParti_Cell_Sorted[np].parti_id;
							dst = 0.0;
							for (k = 0; k < dim; k++)
							{
								dst = dst + (pPar[i].xp[k] - pPar[j].xp[k]) * (pPar[i].xp[k] - pPar[j].xp[k]);
								dx[k] = pPar[i].xp[k] - pPar[j].xp[k];
							}
							dst = sqrt(dst);
							ra = (pPar[i].hl + pPar[j].hl) / 2;
							if (dst <= 2 * ra && i != j)
							{
								pParCell[i].ninflu += 1;
								p = pParCell[i].ninflu;
								if (p <= 99)
								{
									pParCell[i].influ[p] = j;
									weight2(dst, ra, dx, pParCell[i].wij[p], dim);
								}
								else
									t = 5;
							}
						}
					}
				}
			}
		}
	}

	return t;
}

void clNNPS_Fun::kernel_gradient_correction(Particle *pPar, Par_Cell *pParCell, const Para_Pro &pPPro)
{
	int i, j, k, l, nc, nj, ntotal, ndim;
	double xij[3], l_inv_mat[3][3], wij[3];
	double l_soil[3][3], l_struct[3][3];
	double(*l_mat)[3];
	ntotal = pPPro.ntotal;
	ndim = pPPro.ndim;

#pragma omp parallel for schedule(static) private(j, k, l, nc, nj, xij, \
												  l_inv_mat, wij, l_soil, l_struct, l_mat)
	for (i = 0; i < ntotal; i++)
	{
		if (pPar[i].type == 2 || pPar[i].type == 4)
		{
			// initialization of l_mat and l_inv_mat
			for (k = 0; k < 3; k++)
			{
				for (l = 0; l < 3; l++)
				{
					l_soil[k][l] = 0.0;
					l_struct[k][l] = 0.0;
					l_inv_mat[k][l] = 0.0;
				}
			}
			// calculateing l_mat
			nc = pParCell[i].ninflu;
			for (nj = 1; nj <= nc; nj++)
			{
				j = pParCell[i].influ[nj];
				// soil or structure
				if (pPar[i].type == 2 && (pPar[j].type == 2 || pPar[j].type == 0 || pPar[j].type == 4))
				{
					l_mat = l_soil;
					// xij
					for (k = 0; k < ndim; k++)
						xij[k] = pPar[j].xp[k] - pPar[i].xp[k];
					// 2-dimensional problem
					if (ndim == 2)
					{
						l_mat[0][0] = l_mat[0][0] + pPar[j].mass / pPar[j].rho * xij[0] * pParCell[i].wij[nj][0];
						l_mat[0][1] = l_mat[0][1] + pPar[j].mass / pPar[j].rho * xij[0] * pParCell[i].wij[nj][1];
						l_mat[1][0] = l_mat[1][0] + pPar[j].mass / pPar[j].rho * xij[1] * pParCell[i].wij[nj][0];
						l_mat[1][1] = l_mat[1][1] + pPar[j].mass / pPar[j].rho * xij[1] * pParCell[i].wij[nj][1];
					}
					// 3-dimensional problem
					else if (ndim == 3)
					{
						l_mat[0][0] = l_mat[0][0] + pPar[j].mass / pPar[j].rho * xij[0] * pParCell[i].wij[nj][0];
						l_mat[0][1] = l_mat[0][1] + pPar[j].mass / pPar[j].rho * xij[0] * pParCell[i].wij[nj][1];
						l_mat[0][2] = l_mat[0][2] + pPar[j].mass / pPar[j].rho * xij[0] * pParCell[i].wij[nj][2];
						l_mat[1][0] = l_mat[1][0] + pPar[j].mass / pPar[j].rho * xij[1] * pParCell[i].wij[nj][0];
						l_mat[1][1] = l_mat[1][1] + pPar[j].mass / pPar[j].rho * xij[1] * pParCell[i].wij[nj][1];
						l_mat[1][2] = l_mat[1][2] + pPar[j].mass / pPar[j].rho * xij[1] * pParCell[i].wij[nj][2];
						l_mat[2][0] = l_mat[2][0] + pPar[j].mass / pPar[j].rho * xij[2] * pParCell[i].wij[nj][0];
						l_mat[2][1] = l_mat[2][1] + pPar[j].mass / pPar[j].rho * xij[2] * pParCell[i].wij[nj][1];
						l_mat[2][2] = l_mat[2][2] + pPar[j].mass / pPar[j].rho * xij[2] * pParCell[i].wij[nj][2];
					}
				}
				else if (pPar[i].type == 4 && (pPar[j].type == 4 || pPar[j].type == 0 || pPar[j].type == 2))
				{
					l_mat = l_struct;
					// xij
					for (k = 0; k < ndim; k++)
						xij[k] = pPar[j].xp[k] - pPar[i].xp[k];
					// 2-dimensional problem
					if (ndim == 2)
					{
						l_mat[0][0] = l_mat[0][0] + pPar[j].mass / pPar[j].rho * xij[0] * pParCell[i].wij[nj][0];
						l_mat[0][1] = l_mat[0][1] + pPar[j].mass / pPar[j].rho * xij[0] * pParCell[i].wij[nj][1];
						l_mat[1][0] = l_mat[1][0] + pPar[j].mass / pPar[j].rho * xij[1] * pParCell[i].wij[nj][0];
						l_mat[1][1] = l_mat[1][1] + pPar[j].mass / pPar[j].rho * xij[1] * pParCell[i].wij[nj][1];
					}
					// 3-dimensional problem
					else if (ndim == 3)
					{
						l_mat[0][0] = l_mat[0][0] + pPar[j].mass / pPar[j].rho * xij[0] * pParCell[i].wij[nj][0];
						l_mat[0][1] = l_mat[0][1] + pPar[j].mass / pPar[j].rho * xij[0] * pParCell[i].wij[nj][1];
						l_mat[0][2] = l_mat[0][2] + pPar[j].mass / pPar[j].rho * xij[0] * pParCell[i].wij[nj][2];
						l_mat[1][0] = l_mat[1][0] + pPar[j].mass / pPar[j].rho * xij[1] * pParCell[i].wij[nj][0];
						l_mat[1][1] = l_mat[1][1] + pPar[j].mass / pPar[j].rho * xij[1] * pParCell[i].wij[nj][1];
						l_mat[1][2] = l_mat[1][2] + pPar[j].mass / pPar[j].rho * xij[1] * pParCell[i].wij[nj][2];
						l_mat[2][0] = l_mat[2][0] + pPar[j].mass / pPar[j].rho * xij[2] * pParCell[i].wij[nj][0];
						l_mat[2][1] = l_mat[2][1] + pPar[j].mass / pPar[j].rho * xij[2] * pParCell[i].wij[nj][1];
						l_mat[2][2] = l_mat[2][2] + pPar[j].mass / pPar[j].rho * xij[2] * pParCell[i].wij[nj][2];
					}
				}
			}
			// calculating normalized kernel gradient
			nc = pParCell[i].ninflu;
			for (nj = 1; nj <= nc; nj++)
			{
				// soil or structure
				if (pPar[i].type == 2)
				{
					l_mat = l_soil;
					// calculateing l_inv_mat
					if (ndim == 2)
						inverse_mat_2D(l_mat, l_inv_mat);
					else if (ndim == 3)
						inverse_mat_3D(l_mat, l_inv_mat);
					mat_multiply(l_inv_mat, pParCell[i].wij[nj], wij);
					for (k = 0; k < 3; k++)
						pParCell[i].wij[nj][k] = wij[k];
				}
				else if (pPar[i].type == 4)
				{
					l_mat = l_struct;
					// calculateing l_inv_mat
					if (ndim == 2)
						inverse_mat_2D(l_mat, l_inv_mat);
					else if (ndim == 3)
						inverse_mat_3D(l_mat, l_inv_mat);
					mat_multiply(l_inv_mat, pParCell[i].wij[nj], wij);
					for (k = 0; k < 3; k++)
						pParCell[i].wij[nj][k] = wij[k];
				}
			}
		}
	}
}

// get the inverse matrix of a[2][2]
void clNNPS_Fun::inverse_mat_2D(double (*A)[3], double (*A_inv)[3])
{
	double det;

	det = A[0][0] * A[1][1] - A[0][1] * A[1][0];

	if (fabs(det) > 1.0e-5)
	{
		A_inv[0][0] = A[1][1] / det;
		A_inv[0][1] = -A[0][1] / det;

		A_inv[1][0] = -A[1][0] / det;
		A_inv[1][1] = A[0][0] / det;
	}
	else
	{
		A_inv[0][0] = 1.0;
		A_inv[0][1] = 0.0;

		A_inv[1][0] = 0.0;
		A_inv[1][1] = 1.0;
	}
}

// get the inverse matrix of a[3][3]
void clNNPS_Fun::inverse_mat_3D(double (*A)[3], double (*A_inv)[3])
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
	else
	{
		A_inv[0][0] = 1.0;
		A_inv[0][1] = 0.0;
		A_inv[0][2] = 0.0;

		A_inv[1][0] = 0.0;
		A_inv[1][1] = 1.0;
		A_inv[1][2] = 0.0;

		A_inv[2][0] = 0.0;
		A_inv[2][1] = 0.0;
		A_inv[2][2] = 1.0;
	}
}

// merge all the elements by the Cell id
void clNNPS_Fun::Merge_CellID(parti_cellid *arr_A, parti_cellid *arr_temp, int start, int mid_index, int end)
{
	int i, j, k;

	for (i = start; i <= end; i++)
	{ // transport selected part to temp array
		arr_temp[i] = arr_A[i];
	}

	// merge the elements
	i = start;
	j = mid_index + 1;
	k = start;
	while (i < (mid_index + 1) && j <= end)
	{
		if (arr_temp[i].cell_id < arr_temp[j].cell_id)
		{
			arr_A[k++] = arr_temp[i++];
		}
		else
		{
			arr_A[k++] = arr_temp[j++];
		}
	}
	// merge the rest elements
	while (i < (mid_index + 1))
	{
		arr_A[k++] = arr_temp[i++];
	}
	while (j <= end)
	{
		arr_A[k++] = arr_temp[j++];
	}
}

// divide the array
void clNNPS_Fun::Binary_Sort(parti_cellid *pParti_Cell_Sorted, parti_cellid *pParti_Cell_Temp,
							 cell_info *pCell_Info, int ntotal, int ctotal, int size_m, int k)
{
	// sorting by the first number
	int i_k = 1;
	int interval;			   // interval for each group
	int group_num;			   // group number
	int group_id;			   // group id
	int start, mid_index, end; // start, middle and end id of particles
	while (i_k <= k)
	{
		interval = 1 << i_k;
		group_num = size_m / interval;
#pragma omp parallel for private(start, end, mid_index)
		for (group_id = 0; group_id < group_num; group_id++)
		{ // loop for the group_id in the range of group_num
			start = group_id * interval;
			end = (group_id + 1) * interval - 1;
			mid_index = (start + end) >> 1;
			Merge_CellID(pParti_Cell_Sorted, pParti_Cell_Temp, start, mid_index, end);
		}
		i_k++;
	}

	// generating the cell infotmation
	for (int i = 0; i < ntotal; i++)
	{
		int cell_id = pParti_Cell_Sorted[i].cell_id;
		if (i < pCell_Info[cell_id].start)
			pCell_Info[cell_id].start = i;
		if (i > pCell_Info[cell_id].end)
			pCell_Info[cell_id].end = i;
	}
}