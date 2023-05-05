/************************************************************************************
Multi-Phase and Parallelized SPH program
--parallelized and optimized for SMP (Symmetric multiprocessing ) architecture.
--Copyright (c) 2018-2023, Weijie ZHANG, GeoHohai, Hohai University.
************************************************************************************/

#include <math.h>
#include <omp.h>
#include "Vector3.h"
#include "Header_Option.h"
#include "Class_Functions.h"

inline void set_vibration_velocity(double* vx, double* vbx) {
	for (int i = 0; i < 3; i++) {
		vx[i] =  vbx[i];
	}
}

/* check if the particle is approching the problem domain and resist it*/
void clBndy_Fun::check_domain(Particle* pPar, const Cell_Con& pCellc,
	const Para_Boundary& pBou, const Para_Pro& pPPro) {

	double dx1, dx2, dy1, dy2, dz1, dz2;
	int i, flag;
	double dr = pPPro.dr;
	int ntotal = pPPro.ntotal;
	int dim = pPPro.ndim;
	double bmax[3], bmin[3], vbx[3];
	double a = 0.0;
	double b = 0.0;

	bmax[0] = pCellc.xmax;
	bmax[1] = pCellc.ymax;
	bmax[2] = pCellc.zmax;
	bmin[0] = pCellc.xmin;
	bmin[1] = pCellc.ymin;
	bmin[2] = pCellc.zmin;

	vbx[0] = 0.0;
	vbx[1] = 0.0;
	vbx[2] = 0.0;

	flag = 0;
#pragma omp parallel for schedule(static) private(dx1,dx2,dy1,dy2,dz1,dz2)
	for (i = 0; i < ntotal; i++) {
		//velocity particles
		if (pPar[i].type != 0 && pPar[i].type != 7) {
			dx1 = pPar[i].xp[0] - bmax[0];
			dx2 = pPar[i].xp[0] - bmin[0];
			dy1 = pPar[i].xp[1] - bmax[1];
			dy2 = pPar[i].xp[1] - bmin[1];
			dz1 = pPar[i].xp[2] - bmax[2];
			dz2 = pPar[i].xp[2] - bmin[2];
			if (dx1 >= -0.5 * dr) set_vibration_velocity(pPar[i].vx, vbx);
			if (dx2 <= 0.5 * dr)  set_vibration_velocity(pPar[i].vx, vbx);
			if (dy1 >= -0.5 * dr)  set_vibration_velocity(pPar[i].vx, vbx);
			if (dy2 <= 0.5 * dr)  set_vibration_velocity(pPar[i].vx, vbx);
			if (dz1 >= -0.5 * dr && dim == 3)  set_vibration_velocity(pPar[i].vx, vbx);
			if (dz2 <= 0.5 * dr && dim == 3)  set_vibration_velocity(pPar[i].vx, vbx);
		}
		else if (pPar[i].type == 0 && pPar[i].matype > 2 && flag == 0) {
			flag = 1;
			vbx[0] = pPar[i].interfss[0];
			vbx[1] = pPar[i].interfss[1];
			vbx[2] = pPar[i].interfss[2];
		}
	}
}

/* check if the particle is approching the problem domain and resist it*/
void clBndy_Fun::velocity_set_domain(Particle* pPar, Par_Cell* pParCell, cl_bndy_pair* pBndy_Pair, const Para_Pro& pPPro, double coe_bndy) {
	int i, j, k, nc, nj;
	double inprod, dst, beta, db;
	double a, b, tempup1[3], vx2[3];
	int ntotal = pPPro.ntotal;
	int ndim = pPPro.ndim;
	double dr = pPPro.dr;

	a = coe_bndy;
	b = coe_bndy;

#pragma omp parallel for schedule(static) private(j,k,nc,nj,dst,db,inprod,beta,tempup1,vx2)
	for (i = 0; i < ntotal; i++) {
		if (pPar[i].type == 1 || pPar[i].type == 3) {

			for (k = 0; k < 3; k++) {
				tempup1[k] = 0.0;
			}

			//velocity particles
			nc = pParCell[i].ninflu;
			for (nj = 1; nj <= nc; nj++) {

				j = pParCell[i].influ[nj];
				if (pPar[j].type == 0) {
					inprod = 0.0;
					dst = 0.0;
					for (k = 0; k < ndim; k++) {
						inprod = inprod + (pPar[i].vx[k] - pPar[j].interfss[k])
							* (pPar[i].xp[k] - pPar[j].xp[k]);
						dst = dst + (pPar[i].xp[k] - pPar[j].xp[k]) * (pPar[i].xp[k] - pPar[j].xp[k]);
					}
					if (inprod < 0.0) {
						dst = sqrt(dst);
						db = a * dr;
						beta = 1 + db / (dst - b * dr);
						beta = fmin(1.5, beta);

						for (k = 0; k < ndim; k++)
							tempup1[k] = tempup1[k]
							+ (1 - beta) * (pPar[i].xp[k] - pPar[j].xp[k]) * inprod / dst / dst;
					}
				}
			}

			for (k = 0; k < ndim; k++) {
				vx2[k] = pPar[i].vx[k] + tempup1[k];
				pPar[i].vx[k] = vx2[k];
			}
		}
	}
}

/* boundary effect:  using the Takeda & Morris and no permeability*/
void clBndy_Fun::tmboundary(Particle* pPar, Par_Cell* pParCell, cl_bndy_pair* pBndy_Pair, const Para_Pro& pPPro) {

	int i, j, k, nc, nj;
	double inprod, dst, beta, db;
	double a, b, tempup1[3], vx2[3];
	int ntotal = pPPro.ntotal;
	int ndim = pPPro.ndim;
	double dr = pPPro.dr;

	a = 0.30;
	b = 0.30;

#pragma omp parallel for schedule(static) private(j,k,nc,nj,dst,db,inprod,beta,tempup1,vx2)
	for (i = 0; i < ntotal; i++) {
		if (pPar[i].type != 0 && pPar[i].type != 7) {

			for (k = 0; k < 3; k++) {
				tempup1[k] = 0.0;
			}

			//velocity particles
			nc = pParCell[i].ninflu;
			for (nj = 1; nj <= nc; nj++) {

				j = pParCell[i].influ[nj];
				if (pPar[j].type == 0) {
					inprod = 0.0;
					dst = 0.0;
					for (k = 0; k < ndim; k++) {
						inprod = inprod + (pPar[i].vx[k] - pPar[j].interfss[k])
							* (pPar[i].xp[k] - pPar[j].xp[k]);
						dst = dst + (pPar[i].xp[k] - pPar[j].xp[k]) * (pPar[i].xp[k] - pPar[j].xp[k]);
					}
					if (inprod < 0.0) {
						dst = sqrt(dst);
						db = a * dr;
						beta = 1 + db / (dst - b * dr);
						beta = fmin(1.5, beta);

						for (k = 0; k < ndim; k++)
							tempup1[k] = tempup1[k]
							+ (1 - beta) * (pPar[i].xp[k] - pPar[j].xp[k]) * inprod / dst / dst;
					}
				}
			}

			for (k = 0; k < ndim; k++) {
				vx2[k]= pPar[i].vx[k] + tempup1[k];
				pPar[i].vx[k] = vx2[k];
			}
		}
	}
}

void clBndy_Fun::tm_boundary_2d(Particle *pPar, Par_Cell *pParCell, cl_bndy_pair *pBndy_Pair, clVar_Boundary* pParti_VariBndy,
								const Para_Pro &pPPro, int cn)
{
	//variables and their initialization
	int ntotal = pPPro.ntotal;
	int ndim = pPPro.ndim;
	double dr = pPPro.dr;
	int i, j, k, nc, nj, id;
	double dst, da, db, beta;
	int id1, id2;
	double vel_b[3];
	double tempup_vel[4][3], tempdown[4];

	// class variables
	clParti_Pair pPartiPair_SS;

	//geometry variables
	clPoint p1, p2, p4;
	clLine line1;

	//calculating relative distance between cumputing particle and structure particle
#pragma omp parallel for schedule(static) private(j, nc, nj, k, pPartiPair_SS, dst, da, db, beta, id1, id2, \
p1, p2, p4, line1, vel_b, id)
	for (i = 0; i < ntotal; i++)
	{
		if (pPar[i].type != 0 && pPar[i].type != 7)
		{
			//particle information initialization
			pPartiPair_SS.info_reset();

			//searching contacting boundary particles
			k = 0;
			id = 0;
			nc = pParCell[i].ninflu;
			for (nj = 1; nj <= nc; nj++)
			{
				j = pParCell[i].influ[nj];
				if (pPar[j].type == 0)
				{
					pPartiPair_SS.parti_id[k] = j;
					dst = (pPar[i].xp[0] - pPar[j].xp[0]) * (pPar[i].xp[0] - pPar[j].xp[0]) 
						+ (pPar[i].xp[1] - pPar[j].xp[1]) * (pPar[i].xp[1] - pPar[j].xp[1]) 
						+ (pPar[i].xp[2] - pPar[j].xp[2]) * (pPar[i].xp[2] - pPar[j].xp[2]);
					dst = sqrt(dst);
					pPartiPair_SS.parti_dst[k] = dst;
					k += 1;
					if (pPar[j].matype > 2) id = j;
				}
			}
			pPartiPair_SS.total = k;
			pPartiPair_SS.dist_sorting();

			//converting particle information to point and vector
			if (pPartiPair_SS.total >= 2)
			{
				converse_point(pPar[i], &p4);

				//converting particle information to point and vector
				id1 = pPartiPair_SS.parti_id[0];
				converse_point(pPar[id1], &p1);

				id2 = pPartiPair_SS.parti_id[1];
				converse_point(pPar[id2], &p2);

				//calculating line between p1 and p2
				cal_line(p1, p2, &line1);

				//calculating distance between moving particle and boundary line
				cal_dpline(p4, line1, &dst);

				//calculating da
				if (dst > 0.6 * dr)
					da = dst - 0.5* dr;
				else
					da = 0.1* dr;
				db = 0.5 * dr;

				//calculating beta
				if (fabs(da) > 1.0e-6)
				{
					beta = 1.0 + db / da; //calculating beta and compare beta and beta_max
					beta = fmin(1.5, beta);
				}
				else
					beta = 1.5; //setting bete=beta_max

				//calculating virtual velocity
				vel_b[0] = (1.0 - beta) * pPar[i].vx[0] + beta * pPar[id].interfss[0];
				vel_b[1] = (1.0 - beta) * pPar[i].vx[1] + beta * pPar[id].interfss[1];
				vel_b[2] = 0.0;

				//setting virtual velocity
				pBndy_Pair[i].vel_bndy[0] = vel_b[0];
				pBndy_Pair[i].vel_bndy[1] = vel_b[1];
				pBndy_Pair[i].vel_bndy[2] = vel_b[2];
			}
			else {
				pBndy_Pair[i].vel_bndy[0] = pPar[i].vx[0];
				pBndy_Pair[i].vel_bndy[1] = pPar[i].vx[1];
				pBndy_Pair[i].vel_bndy[2] = pPar[i].vx[2];
			}
		}
	}
#pragma omp parallel for schedule(static) private(j, nc, nj, tempup_vel, tempdown)
	for (i = 0; i < ntotal; i++)
	{
		if (pPar[i].type == 0)
		{
			tempup_vel[0][0] = 0.0;
			tempup_vel[0][1] = 0.0;
			tempup_vel[0][2] = 0.0;
			tempup_vel[1][0] = 0.0;
			tempup_vel[1][1] = 0.0;
			tempup_vel[1][2] = 0.0;
			tempup_vel[2][0] = 0.0;
			tempup_vel[2][1] = 0.0;
			tempup_vel[2][2] = 0.0;
			tempup_vel[3][0] = 0.0;
			tempup_vel[3][1] = 0.0;
			tempup_vel[3][2] = 0.0;

			tempdown[0] = 0.0;
			tempdown[1] = 0.0;
			tempdown[2] = 0.0;
			tempdown[3] = 0.0;

			nc = pParCell[i].ninflu;
			//calculating every db and vb for each boundary particles in supporting domain
			for (nj = 1; nj <= nc; nj++)
			{
				j = pParCell[i].influ[nj];
				if (pPar[j].type == 1) {
					tempup_vel[0][0] = tempup_vel[0][0] + pPar[j].mass * pBndy_Pair[j].vel_bndy[0] * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempup_vel[0][1] = tempup_vel[0][0] + pPar[j].mass * pBndy_Pair[j].vel_bndy[1] * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempup_vel[0][2] = tempup_vel[0][0] + pPar[j].mass * pBndy_Pair[j].vel_bndy[2] * pParCell[i].wij[nj][3] / pPar[j].rho;
					
					tempdown[0] = tempdown[0] + pPar[j].mass * pParCell[i].wij[nj][3] / pPar[j].rho;
				}
				else if (pPar[j].type == 2) {
					tempup_vel[1][0] = tempup_vel[1][0] + pPar[j].mass * pBndy_Pair[j].vel_bndy[0] * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempup_vel[1][1] = tempup_vel[1][0] + pPar[j].mass * pBndy_Pair[j].vel_bndy[1] * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempup_vel[1][2] = tempup_vel[1][0] + pPar[j].mass * pBndy_Pair[j].vel_bndy[2] * pParCell[i].wij[nj][3] / pPar[j].rho;

					tempdown[1] = tempdown[1] + pPar[j].mass * pParCell[i].wij[nj][3] / pPar[j].rho;
				}
				else if (pPar[j].type == 3) {
					tempup_vel[2][0] = tempup_vel[2][0] + pPar[j].mass * pBndy_Pair[j].vel_bndy[0] * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempup_vel[2][1] = tempup_vel[2][0] + pPar[j].mass * pBndy_Pair[j].vel_bndy[1] * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempup_vel[2][2] = tempup_vel[2][0] + pPar[j].mass * pBndy_Pair[j].vel_bndy[2] * pParCell[i].wij[nj][3] / pPar[j].rho;

					tempdown[2] = tempdown[2] + pPar[j].mass * pParCell[i].wij[nj][3] / pPar[j].rho;
				}
				else if (pPar[j].type == 4) {
					tempup_vel[3][0] = tempup_vel[3][0] + pPar[j].mass * pBndy_Pair[j].vel_bndy[0] * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempup_vel[3][1] = tempup_vel[3][0] + pPar[j].mass * pBndy_Pair[j].vel_bndy[1] * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempup_vel[3][2] = tempup_vel[3][0] + pPar[j].mass * pBndy_Pair[j].vel_bndy[2] * pParCell[i].wij[nj][3] / pPar[j].rho;

					tempdown[3] = tempdown[3] + pPar[j].mass * pParCell[i].wij[nj][3] / pPar[j].rho;
				}
			}

			//calculating velocity  of boundary particles
			if (fabs(tempdown[0]) > 1.0e-6) { //water
				pParti_VariBndy[i].vx_water[0] = tempup_vel[0][0] / tempdown[0];
				pParti_VariBndy[i].vx_water[1] = tempup_vel[0][1] / tempdown[0];
				pParti_VariBndy[i].vx_water[2] = tempup_vel[0][2] / tempdown[0];
			}
			else {
				pParti_VariBndy[i].vx_water[0] = 0.0;
				pParti_VariBndy[i].vx_water[1] = 0.0;
				pParti_VariBndy[i].vx_water[2] = 0.0;
			}

			if (fabs(tempdown[1]) > 1.0e-6) { //soil
				pParti_VariBndy[i].vx_soil[0] = tempup_vel[1][0] / tempdown[1];
				pParti_VariBndy[i].vx_soil[1] = tempup_vel[1][1] / tempdown[1];
				pParti_VariBndy[i].vx_soil[2] = tempup_vel[1][2] / tempdown[1];
			}
			else {
				pParti_VariBndy[i].vx_soil[0] = 0.0;
				pParti_VariBndy[i].vx_soil[1] = 0.0;
				pParti_VariBndy[i].vx_soil[2] = 0.0;
			}

			if (fabs(tempdown[2]) > 1.0e-6) { //air
				pParti_VariBndy[i].vx_air[0] = tempup_vel[2][0] / tempdown[2];
				pParti_VariBndy[i].vx_air[1] = tempup_vel[2][1] / tempdown[2];
				pParti_VariBndy[i].vx_air[2] = tempup_vel[2][2] / tempdown[2];
			}
			else {
				pParti_VariBndy[i].vx_air[0] = 0.0;
				pParti_VariBndy[i].vx_air[1] = 0.0;
				pParti_VariBndy[i].vx_air[2] = 0.0;
			}

			if (fabs(tempdown[3]) > 1.0e-6) { //structure
				pParti_VariBndy[i].vx_struct[0] = tempup_vel[3][0] / tempdown[3];
				pParti_VariBndy[i].vx_struct[1] = tempup_vel[3][1] / tempdown[3];
				pParti_VariBndy[i].vx_struct[2] = tempup_vel[3][2] / tempdown[3];
			}
			else {
				pParti_VariBndy[i].vx_struct[0] = 0.0;
				pParti_VariBndy[i].vx_struct[1] = 0.0;
				pParti_VariBndy[i].vx_struct[2] = 0.0;
			}
		}
	}
}

void clBndy_Fun::tm_boundary_3d(Particle *pPar, Par_Cell *pParCell, cl_bndy_pair *pBndy_Pair, clVar_Boundary* pParti_VariBndy,
								const Para_Pro &pPPro, int cn)
{
	//variables and their initialization
	int ntotal = pPPro.ntotal;
	int ndim = pPPro.ndim;
	double dr = pPPro.dr;
	int i, j, k, nc, nj, id;
	double dst, da, db;
	int id1, id2, id3;
	double beta, vel_b[3];
	double cos_sita;
	double tempup_vel[4][3], tempdown[4];

	// class variables
	clParti_Pair pPartiPair_SS;

	//geometry variables
	clPoint p1, p2, p3, p4;
	clPlane plane1;
	clPlane plane_c1, plane_c2;

	//calculating relative distance between cumputing particle and structure particle
#pragma omp parallel for schedule(static) private(j, nc, nj, k, dst, da, db, pPartiPair_SS, \
id1, id2, id3, p1, p2, p3, p4, plane1, plane_c1, plane_c2, beta, vel_b, cos_sita, id)
	for (i = 0; i < ntotal; i++)
	{
		if (pPar[i].type != 0 && pPar[i].type != 7)
		{
			//particle information initialization
			pPartiPair_SS.info_reset();

			//searching contacting boundary particles
			k = 0;
			id = 0;
			nc = pParCell[i].ninflu;
			for (nj = 1; nj <= nc; nj++)
			{
				j = pParCell[i].influ[nj];
				if (pPar[j].type == 0)
				{
					pPartiPair_SS.parti_id[k] = j;
					dst = (pPar[i].xp[0] - pPar[j].xp[0]) * (pPar[i].xp[0] - pPar[j].xp[0])
						+ (pPar[i].xp[1] - pPar[j].xp[1]) * (pPar[i].xp[1] - pPar[j].xp[1])
						+ (pPar[i].xp[2] - pPar[j].xp[2]) * (pPar[i].xp[2] - pPar[j].xp[2]);
					dst = sqrt(dst);
					pPartiPair_SS.parti_dst[k] = dst;
					k += 1;
					if (pPar[j].matype > 2) id = j;
				}
			}
			pPartiPair_SS.total = k;
			pPartiPair_SS.dist_sorting();

			//converting particle information to point and vector
			if (pPartiPair_SS.total >= 3)
			{
				converse_point(pPar[i], &p4);

				//converting particle information to point and vector
				id1 = pPartiPair_SS.parti_id[0];
				converse_point(pPar[id1], &p1);

				id2 = pPartiPair_SS.parti_id[1];
				converse_point(pPar[id2], &p2);

				id3 = pPartiPair_SS.parti_id[2];
				converse_point(pPar[id3], &p3);

				//calculating plane between p1, p2, and p3
				cal_plane(p1, p2, p4, &plane_c1);
				cal_plane(p1, p3, p4, &plane_c2);

				cos_sita = plane_c1.A * plane_c2.A + plane_c1.B * plane_c2.B + plane_c1.C * plane_c2.C;
				cos_sita = cos_sita / sqrt(plane_c1.A * plane_c1.A + plane_c1.B * plane_c1.B + plane_c1.C * plane_c1.C)
					/ sqrt(plane_c2.A * plane_c2.A + plane_c2.B * plane_c2.B + plane_c2.C * plane_c2.C);

				if (fabs(cos_sita) > 0.999999)
				{
					//three nerest boundary particles in a straight line
					plane_two_points(p4, p1, &plane1);
				}
				else
					cal_plane(p1, p2, p3, &plane1);

				//calculating distance between moving particle and boundary plane
				cal_dpplane(p4, plane1, &dst);

				//calculating da
				if (dst > 0.6 * dr)
					da = dst - 0.5* dr;
				else
					da = 0.1* dr;
				db = 0.5 * dr;

				//calculating beta
				if (fabs(da) > 1.0e-6)
				{
					beta = 1.0 + db / da; //calculating beta and compare beta and beta_max
					beta = fmin(1.5, beta);
				}
				else
					beta = 1.5; //setting bete=beta_max

				//calculating virtual velocity
				vel_b[0] = (1.0 - beta) * pPar[i].vx[0] + beta * pPar[id].interfss[0];
				vel_b[1] = (1.0 - beta) * pPar[i].vx[1] + beta * pPar[id].interfss[1];
				vel_b[2] = (1.0 - beta) * pPar[i].vx[2] + beta * pPar[id].interfss[2];

				//setting virtual velocity
				pBndy_Pair[i].vel_bndy[0] = vel_b[0];
				pBndy_Pair[i].vel_bndy[1] = vel_b[1];
				pBndy_Pair[i].vel_bndy[2] = vel_b[2];
			}
			else {
				pBndy_Pair[i].vel_bndy[0] = pPar[i].vx[0];
				pBndy_Pair[i].vel_bndy[1] = pPar[i].vx[1];
				pBndy_Pair[i].vel_bndy[2] = pPar[i].vx[2];
			}
		}
	}

#pragma omp parallel for schedule(static) private(j, nc, nj, tempup_vel, tempdown)
	for (i = 0; i < ntotal; i++)
	{
		if (pPar[i].type == 0)
		{
			tempup_vel[0][0] = 0.0;
			tempup_vel[0][1] = 0.0;
			tempup_vel[0][2] = 0.0;
			tempup_vel[1][0] = 0.0;
			tempup_vel[1][1] = 0.0;
			tempup_vel[1][2] = 0.0;
			tempup_vel[2][0] = 0.0;
			tempup_vel[2][1] = 0.0;
			tempup_vel[2][2] = 0.0;
			tempup_vel[3][0] = 0.0;
			tempup_vel[3][1] = 0.0;
			tempup_vel[3][2] = 0.0;

			tempdown[0] = 0.0;
			tempdown[1] = 0.0;
			tempdown[2] = 0.0;
			tempdown[3] = 0.0;

			nc = pParCell[i].ninflu;
			//calculating every db and vb for each boundary particles in supporting domain
			for (nj = 1; nj <= nc; nj++)
			{
				j = pParCell[i].influ[nj];
				if (pPar[j].type == 1) {
					tempup_vel[0][0] = tempup_vel[0][0] + pPar[j].mass * pBndy_Pair[j].vel_bndy[0] * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempup_vel[0][1] = tempup_vel[0][0] + pPar[j].mass * pBndy_Pair[j].vel_bndy[1] * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempup_vel[0][2] = tempup_vel[0][0] + pPar[j].mass * pBndy_Pair[j].vel_bndy[2] * pParCell[i].wij[nj][3] / pPar[j].rho;
					
					tempdown[0] = tempdown[0] + pPar[j].mass * pParCell[i].wij[nj][3] / pPar[j].rho;
				}
				else if (pPar[j].type == 2) {
					tempup_vel[1][0] = tempup_vel[1][0] + pPar[j].mass * pBndy_Pair[j].vel_bndy[0] * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempup_vel[1][1] = tempup_vel[1][0] + pPar[j].mass * pBndy_Pair[j].vel_bndy[1] * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempup_vel[1][2] = tempup_vel[1][0] + pPar[j].mass * pBndy_Pair[j].vel_bndy[2] * pParCell[i].wij[nj][3] / pPar[j].rho;

					tempdown[1] = tempdown[1] + pPar[j].mass * pParCell[i].wij[nj][3] / pPar[j].rho;
				}
				else if (pPar[j].type == 3) {
					tempup_vel[2][0] = tempup_vel[2][0] + pPar[j].mass * pBndy_Pair[j].vel_bndy[0] * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempup_vel[2][1] = tempup_vel[2][0] + pPar[j].mass * pBndy_Pair[j].vel_bndy[1] * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempup_vel[2][2] = tempup_vel[2][0] + pPar[j].mass * pBndy_Pair[j].vel_bndy[2] * pParCell[i].wij[nj][3] / pPar[j].rho;

					tempdown[2] = tempdown[2] + pPar[j].mass * pParCell[i].wij[nj][3] / pPar[j].rho;
				}
				else if (pPar[j].type == 4) {
					tempup_vel[3][0] = tempup_vel[3][0] + pPar[j].mass * pBndy_Pair[j].vel_bndy[0] * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempup_vel[3][1] = tempup_vel[3][0] + pPar[j].mass * pBndy_Pair[j].vel_bndy[1] * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempup_vel[3][2] = tempup_vel[3][0] + pPar[j].mass * pBndy_Pair[j].vel_bndy[2] * pParCell[i].wij[nj][3] / pPar[j].rho;

					tempdown[3] = tempdown[3] + pPar[j].mass * pParCell[i].wij[nj][3] / pPar[j].rho;
				}
			}

			//calculating velocity  of boundary particles
			if (fabs(tempdown[0]) > 1.0e-6) { //water
				pParti_VariBndy[i].vx_water[0] = tempup_vel[0][0] / tempdown[0];
				pParti_VariBndy[i].vx_water[1] = tempup_vel[0][1] / tempdown[0];
				pParti_VariBndy[i].vx_water[2] = tempup_vel[0][2] / tempdown[0];
			}
			else {
				pParti_VariBndy[i].vx_water[0] = 0.0;
				pParti_VariBndy[i].vx_water[1] = 0.0;
				pParti_VariBndy[i].vx_water[2] = 0.0;
			}

			if (fabs(tempdown[1]) > 1.0e-6) { //soil
				pParti_VariBndy[i].vx_soil[0] = tempup_vel[1][0] / tempdown[1];
				pParti_VariBndy[i].vx_soil[1] = tempup_vel[1][1] / tempdown[1];
				pParti_VariBndy[i].vx_soil[2] = tempup_vel[1][2] / tempdown[1];
			}
			else {
				pParti_VariBndy[i].vx_soil[0] = 0.0;
				pParti_VariBndy[i].vx_soil[1] = 0.0;
				pParti_VariBndy[i].vx_soil[2] = 0.0;
			}

			if (fabs(tempdown[2]) > 1.0e-6) { //air
				pParti_VariBndy[i].vx_air[0] = tempup_vel[2][0] / tempdown[2];
				pParti_VariBndy[i].vx_air[1] = tempup_vel[2][1] / tempdown[2];
				pParti_VariBndy[i].vx_air[2] = tempup_vel[2][2] / tempdown[2];
			}
			else {
				pParti_VariBndy[i].vx_air[0] = 0.0;
				pParti_VariBndy[i].vx_air[1] = 0.0;
				pParti_VariBndy[i].vx_air[2] = 0.0;
			}

			if (fabs(tempdown[3]) > 1.0e-6) { //structure
				pParti_VariBndy[i].vx_struct[0] = tempup_vel[3][0] / tempdown[3];
				pParti_VariBndy[i].vx_struct[1] = tempup_vel[3][1] / tempdown[3];
				pParti_VariBndy[i].vx_struct[2] = tempup_vel[3][2] / tempdown[3];
			}
			else {
				pParti_VariBndy[i].vx_struct[0] = 0.0;
				pParti_VariBndy[i].vx_struct[1] = 0.0;
				pParti_VariBndy[i].vx_struct[2] = 0.0;
			}
		}
	}
}

//calculating distance of a point to a line
void clBndy_Fun::cal_dpline(clPoint ref_point, clLine inLine, double *dp)
{
	double temp, dst;

	dst = fabs(inLine.A * ref_point.x + inLine.B * ref_point.y + inLine.C);
	temp = sqrt(inLine.A * inLine.A + inLine.B * inLine.B);

	if (temp > ep_err)
		*dp = dst / temp;
	else
		*dp = 0.0;
}

//calculating distance of a point to a plane
void clBndy_Fun::cal_dpplane(clPoint ref_point, clPlane inPlane, double *dp)
{
	double temp, dst;

	dst = fabs(inPlane.A * ref_point.x + inPlane.B * ref_point.y + inPlane.C * ref_point.z + inPlane.D);
	temp = sqrt(inPlane.A * inPlane.A + inPlane.B * inPlane.B + inPlane.C * inPlane.C);

	if (temp > ep_err)
		*dp = dst / temp;
	else
		*dp = 0.0;
}

//calculating line parameters Ax+By+C=0
void clBndy_Fun::cal_line(clPoint p1, clPoint p2, clLine *outLine)
{
	double dst = p2.y * p1.x - p1.y * p2.x;
	outLine->A = (p1.y - p2.y) / dst;	
	outLine->B = (p2.x - p1.x) / dst;
}

//calculating plane parameters Ax+By+Cz+D=0
void clBndy_Fun::cal_plane(clPoint p1, clPoint p2, clPoint p3, clPlane *outPlane)
{
	double dst =
		-p1.x * (p2.y * p3.z - p3.y * p2.z) - p2.x * (p3.y * p1.z - p1.y * p3.z) - p3.x * (p1.y * p2.z - p2.y * p1.z);

	outPlane->A = p1.y * (p2.z - p3.z) + p2.y * (p3.z - p1.z) + p3.y * (p1.z - p2.z);
	outPlane->B = p1.z * (p2.x - p3.x) + p2.z * (p3.x - p1.x) + p3.z * (p1.x - p2.x);
	outPlane->C = p1.x * (p2.y - p3.y) + p2.x * (p3.y - p1.y) + p3.x * (p1.y - p2.y);

	outPlane->A = outPlane->A / dst;
	outPlane->B = outPlane->B / dst;
	outPlane->C = outPlane->C / dst;
}

//converse the particle information to clPoint type
void clBndy_Fun::converse_point(Particle inPar, clPoint *outPoint)
{
	outPoint->x = inPar.xp[0];
	outPoint->y = inPar.xp[1];
	outPoint->z = inPar.xp[2];
}

//converse the particle information to Vector3 type
void clBndy_Fun::converse_vector(Particle inPar, Vector3 *outVect)
{
	outVect->x = inPar.xp[0];
	outVect->y = inPar.xp[1];
	outVect->z = inPar.xp[2];
}

//calculate plane from two points across p2
void clBndy_Fun::plane_two_points(clPoint p1, clPoint p2_across, clPlane *outplane)
{
	double temp;
	double a, b, c;

	a = p1.x - p2_across.x;
	b = p1.y - p2_across.y;
	c = p1.z - p2_across.z;

	temp = -a * p2_across.x - b * p2_across.y - c * p2_across.z;
	if (fabs(temp) > 1.0e-6)
	{
		outplane->A = a / temp;
		outplane->B = b / temp;
		outplane->C = c / temp;
		outplane->D = 1.0;
	}
	else
	{
		outplane->A = a;
		outplane->B = b;
		outplane->C = c;
		outplane->D = 0.0;
	}
}

//no-slip boundary by Tran et al. 2019 Computers and Geotechnics
void clBndy_Fun::no_slip_tran(Particle* pPar, Par_Cell* pParCell, clVar_Boundary* pParti_VariBndy,
	const Para_Pro& pPPro, int cn) {
	//local variables
	int ntotal = pPPro.ntotal;
	int ndim = pPPro.ndim;
	int i, j, k, nc, nj;
	double tempup_sig[4][6], tempup_vel[4][3], tempup_rho[4], tempdown[4];
	double satu;

	//codes
#pragma omp parallel for schedule(static) private(j, nc, nj, k, tempup_sig, tempup_vel,\
 tempup_rho, tempdown, satu)
	for (i = 0; i < ntotal; i++)
	{
		//initialization
		for (k = 0; k < 6; k++) {
			tempup_sig[0][k] = 0.0;
			tempup_sig[1][k] = 0.0;
			tempup_sig[2][k] = 0.0;
			tempup_sig[3][k] = 0.0;
		}
		tempup_vel[0][0] = 0.0;
		tempup_vel[0][1] = 0.0;
		tempup_vel[0][2] = 0.0;
		tempup_vel[1][0] = 0.0;
		tempup_vel[1][1] = 0.0;
		tempup_vel[1][2] = 0.0;
		tempup_vel[2][0] = 0.0;
		tempup_vel[2][1] = 0.0;
		tempup_vel[2][2] = 0.0;
		tempup_vel[3][0] = 0.0;
		tempup_vel[3][1] = 0.0;
		tempup_vel[3][2] = 0.0;

		tempup_rho[0] = 0.0;
		tempup_rho[1] = 0.0;
		tempup_rho[2] = 0.0;
		tempup_rho[3] = 0.0;

		tempdown[0] = 0.0;
		tempdown[1] = 0.0;
		tempdown[2] = 0.0;
		tempdown[3] = 0.0;

		if (pPar[i].type == 0) {
			//calculating tempup_sig[4][6], tempup_vel[4][3], tempup_rho[4], and tempdown[4];
			nc = pParCell[i].ninflu;
			for (nj = 1; nj <= nc; nj++) {
				j = pParCell[i].influ[nj];

				//water
				if (pPar[j].type == 1) {
					tempup_sig[0][0] = tempup_sig[0][0] + pPar[j].mass * pPar[j].sig[0] * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempup_sig[0][1] = tempup_sig[0][1] + pPar[j].mass * pPar[j].sig[1] * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempup_sig[0][2] = tempup_sig[0][2] + pPar[j].mass * pPar[j].sig[2] * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempup_sig[0][3] = tempup_sig[0][3] + pPar[j].mass * pPar[j].sig[3] * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempup_sig[0][4] = tempup_sig[0][4] + pPar[j].mass * pPar[j].sig[4] * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempup_sig[0][5] = tempup_sig[0][5] + pPar[j].mass * pPar[j].sig[5] * pParCell[i].wij[nj][3] / pPar[j].rho;

					tempup_vel[0][0] = tempup_vel[0][0] + pPar[j].mass * pPar[j].vx[0] * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempup_vel[0][1] = tempup_vel[0][1] + pPar[j].mass * pPar[j].vx[1] * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempup_vel[0][2] = tempup_vel[0][2] + pPar[j].mass * pPar[j].vx[2] * pParCell[i].wij[nj][3] / pPar[j].rho;

					tempup_rho[0] = tempup_rho[0] + pPar[j].mass * pPar[j].rho * pParCell[i].wij[nj][3] / pPar[j].rho;

					tempdown[0] = tempdown[0] + pPar[j].mass * pParCell[i].wij[nj][3] / pPar[j].rho;
				}

				//soil
				if (pPar[j].type == 2) {
					satu = pPar[j].satu;
					tempup_sig[1][0] = tempup_sig[1][0] + 
						pPar[j].mass * (pPar[j].sig[0] - satu * pPar[j].prew + (1 - satu) * pPar[j].prea) * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempup_sig[1][1] = tempup_sig[1][1] + 
						pPar[j].mass * (pPar[j].sig[1] - satu * pPar[j].prew + (1 - satu) * pPar[j].prea) * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempup_sig[1][2] = tempup_sig[1][2] + 
						pPar[j].mass * (pPar[j].sig[2] - satu * pPar[j].prew + (1 - satu) * pPar[j].prea) * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempup_sig[1][3] = tempup_sig[1][3] + pPar[j].mass * pPar[j].sig[3] * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempup_sig[1][4] = tempup_sig[1][4] + pPar[j].mass * pPar[j].sig[4] * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempup_sig[1][5] = tempup_sig[1][5] + pPar[j].mass * pPar[j].sig[5] * pParCell[i].wij[nj][3] / pPar[j].rho;

					tempup_vel[1][0] = tempup_vel[1][0] + pPar[j].mass * pPar[j].vx[0] * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempup_vel[1][1] = tempup_vel[1][1] + pPar[j].mass * pPar[j].vx[1] * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempup_vel[1][2] = tempup_vel[1][2] + pPar[j].mass * pPar[j].vx[2] * pParCell[i].wij[nj][3] / pPar[j].rho;

					tempup_rho[1] = tempup_rho[1] + pPar[j].mass * pPar[j].rho * pParCell[i].wij[nj][3] / pPar[j].rho;

					tempdown[1] = tempdown[1] + pPar[j].mass * pParCell[i].wij[nj][3] / pPar[j].rho;
				}

				//air
				if (pPar[j].type == 3) {
					tempup_sig[2][0] = tempup_sig[2][0] + pPar[j].mass * pPar[j].sig[0] * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempup_sig[2][1] = tempup_sig[2][1] + pPar[j].mass * pPar[j].sig[1] * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempup_sig[2][2] = tempup_sig[2][2] + pPar[j].mass * pPar[j].sig[2] * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempup_sig[2][3] = tempup_sig[2][3] + pPar[j].mass * pPar[j].sig[3] * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempup_sig[2][4] = tempup_sig[2][4] + pPar[j].mass * pPar[j].sig[4] * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempup_sig[2][5] = tempup_sig[2][5] + pPar[j].mass * pPar[j].sig[5] * pParCell[i].wij[nj][3] / pPar[j].rho;

					tempup_vel[2][0] = tempup_vel[2][0] + pPar[j].mass * pPar[j].vx[0] * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempup_vel[2][1] = tempup_vel[2][1] + pPar[j].mass * pPar[j].vx[1] * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempup_vel[2][2] = tempup_vel[2][2] + pPar[j].mass * pPar[j].vx[2] * pParCell[i].wij[nj][3] / pPar[j].rho;

					tempup_rho[2] = tempup_rho[2] + pPar[j].mass * pPar[j].rho * pParCell[i].wij[nj][3] / pPar[j].rho;

					tempdown[2] = tempdown[2] + pPar[j].mass * pParCell[i].wij[nj][3] / pPar[j].rho;
				}

				//structure
				if (pPar[j].type == 4) {
					tempup_sig[3][0] = tempup_sig[3][0] + pPar[j].mass * pPar[j].sig[0] * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempup_sig[3][1] = tempup_sig[3][1] + pPar[j].mass * pPar[j].sig[1] * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempup_sig[3][2] = tempup_sig[3][2] + pPar[j].mass * pPar[j].sig[2] * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempup_sig[3][3] = tempup_sig[3][3] + pPar[j].mass * pPar[j].sig[3] * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempup_sig[3][4] = tempup_sig[3][4] + pPar[j].mass * pPar[j].sig[4] * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempup_sig[3][5] = tempup_sig[3][5] + pPar[j].mass * pPar[j].sig[5] * pParCell[i].wij[nj][3] / pPar[j].rho;

					tempup_vel[3][0] = tempup_vel[3][0] + pPar[j].mass * pPar[j].vx[0] * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempup_vel[3][1] = tempup_vel[3][1] + pPar[j].mass * pPar[j].vx[1] * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempup_vel[3][2] = tempup_vel[3][2] + pPar[j].mass * pPar[j].vx[2] * pParCell[i].wij[nj][3] / pPar[j].rho;

					tempup_rho[3] = tempup_rho[3] + pPar[j].mass * pPar[j].rho * pParCell[i].wij[nj][3] / pPar[j].rho;

					tempdown[3] = tempdown[3] + pPar[j].mass * pParCell[i].wij[nj][3] / pPar[j].rho;
				}

			}

			//calculating stress tensor, velocity and density of boundary particles
			if (fabs(tempdown[0]) > 1.0e-6) { //water
				pParti_VariBndy[i].sig_water[0] = tempup_sig[0][0] / tempdown[0];
				pParti_VariBndy[i].sig_water[1] = tempup_sig[0][1] / tempdown[0];
				pParti_VariBndy[i].sig_water[2] = tempup_sig[0][2] / tempdown[0];
				pParti_VariBndy[i].sig_water[3] = tempup_sig[0][3] / tempdown[0];
				pParti_VariBndy[i].sig_water[4] = tempup_sig[0][4] / tempdown[0];
				pParti_VariBndy[i].sig_water[5] = tempup_sig[0][5] / tempdown[0];

				pParti_VariBndy[i].vx_water[0] = pPar[i].interfss[0] - tempup_vel[0][0] / tempdown[0];
				pParti_VariBndy[i].vx_water[1] = pPar[i].interfss[1] - tempup_vel[0][1] / tempdown[0];
				pParti_VariBndy[i].vx_water[2] = pPar[i].interfss[2] - tempup_vel[0][2] / tempdown[0];

				pParti_VariBndy[i].rho_water = tempup_rho[0] / tempdown[0];
			}
			else {
				pParti_VariBndy[i].sig_water[0] = 0.0;
				pParti_VariBndy[i].sig_water[1] = 0.0;
				pParti_VariBndy[i].sig_water[2] = 0.0;
				pParti_VariBndy[i].sig_water[3] = 0.0;
				pParti_VariBndy[i].sig_water[4] = 0.0;
				pParti_VariBndy[i].sig_water[5] = 0.0;

				pParti_VariBndy[i].vx_water[0] = pPar[i].interfss[0];
				pParti_VariBndy[i].vx_water[1] = pPar[i].interfss[1];
				pParti_VariBndy[i].vx_water[2] = pPar[i].interfss[2];

				pParti_VariBndy[i].rho_water = 0.0;
			}

			if (fabs(tempdown[1]) > 1.0e-6) { //soil
				pParti_VariBndy[i].sig_soil[0] = tempup_sig[1][0] / tempdown[1];
				pParti_VariBndy[i].sig_soil[1] = tempup_sig[1][1] / tempdown[1];
				pParti_VariBndy[i].sig_soil[2] = tempup_sig[1][2] / tempdown[1];
				pParti_VariBndy[i].sig_soil[3] = tempup_sig[1][3] / tempdown[1];
				pParti_VariBndy[i].sig_soil[4] = tempup_sig[1][4] / tempdown[1];
				pParti_VariBndy[i].sig_soil[5] = tempup_sig[1][5] / tempdown[1];

				pParti_VariBndy[i].vx_soil[0] = pPar[i].interfss[0] - tempup_vel[1][0] / tempdown[1];
				pParti_VariBndy[i].vx_soil[1] = pPar[i].interfss[1] - tempup_vel[1][1] / tempdown[1];
				pParti_VariBndy[i].vx_soil[2] = pPar[i].interfss[2] - tempup_vel[1][2] / tempdown[1];

				pParti_VariBndy[i].rho_soil = tempup_rho[1] / tempdown[1];
			}
			else {
				pParti_VariBndy[i].sig_soil[0] = 0.0;
				pParti_VariBndy[i].sig_soil[1] = 0.0;
				pParti_VariBndy[i].sig_soil[2] = 0.0;
				pParti_VariBndy[i].sig_soil[3] = 0.0;
				pParti_VariBndy[i].sig_soil[4] = 0.0;
				pParti_VariBndy[i].sig_soil[5] = 0.0;

				pParti_VariBndy[i].vx_soil[0] = pPar[i].interfss[0];
				pParti_VariBndy[i].vx_soil[1] = pPar[i].interfss[1];
				pParti_VariBndy[i].vx_soil[2] = pPar[i].interfss[2];

				pParti_VariBndy[i].rho_soil = 0.0;
			}

			if (fabs(tempdown[2]) > 1.0e-6) { //air
				pParti_VariBndy[i].sig_air[0] = tempup_sig[2][0] / tempdown[2];
				pParti_VariBndy[i].sig_air[1] = tempup_sig[2][1] / tempdown[2];
				pParti_VariBndy[i].sig_air[2] = tempup_sig[2][2] / tempdown[2];
				pParti_VariBndy[i].sig_air[3] = tempup_sig[2][3] / tempdown[2];
				pParti_VariBndy[i].sig_air[4] = tempup_sig[2][4] / tempdown[2];
				pParti_VariBndy[i].sig_air[5] = tempup_sig[2][5] / tempdown[2];

				pParti_VariBndy[i].vx_air[0] = pPar[i].interfss[0] - tempup_vel[2][0] / tempdown[2];
				pParti_VariBndy[i].vx_air[1] = pPar[i].interfss[1] - tempup_vel[2][1] / tempdown[2];
				pParti_VariBndy[i].vx_air[2] = pPar[i].interfss[2] - tempup_vel[2][2] / tempdown[2];

				pParti_VariBndy[i].rho_air = tempup_rho[2] / tempdown[2];
			}
			else {
				pParti_VariBndy[i].sig_air[0] = 0.0;
				pParti_VariBndy[i].sig_air[1] = 0.0;
				pParti_VariBndy[i].sig_air[2] = 0.0;
				pParti_VariBndy[i].sig_air[3] = 0.0;
				pParti_VariBndy[i].sig_air[4] = 0.0;
				pParti_VariBndy[i].sig_air[5] = 0.0;

				pParti_VariBndy[i].vx_air[0] = pPar[i].interfss[0];
				pParti_VariBndy[i].vx_air[1] = pPar[i].interfss[1];
				pParti_VariBndy[i].vx_air[2] = pPar[i].interfss[2];

				pParti_VariBndy[i].rho_air = 0.0;
			}

			if (fabs(tempdown[3]) > 1.0e-6) { //structure
				pParti_VariBndy[i].sig_struct[0] = tempup_sig[3][0] / tempdown[3];
				pParti_VariBndy[i].sig_struct[1] = tempup_sig[3][1] / tempdown[3];
				pParti_VariBndy[i].sig_struct[2] = tempup_sig[3][2] / tempdown[3];
				pParti_VariBndy[i].sig_struct[3] = tempup_sig[3][3] / tempdown[3];
				pParti_VariBndy[i].sig_struct[4] = tempup_sig[3][4] / tempdown[3];
				pParti_VariBndy[i].sig_struct[5] = tempup_sig[3][5] / tempdown[3];

				pParti_VariBndy[i].vx_struct[0] = pPar[i].interfss[0] - tempup_vel[3][0] / tempdown[3];
				pParti_VariBndy[i].vx_struct[1] = pPar[i].interfss[1] - tempup_vel[3][1] / tempdown[3];
				pParti_VariBndy[i].vx_struct[2] = pPar[i].interfss[2] - tempup_vel[3][2] / tempdown[3];

				pParti_VariBndy[i].rho_struct = tempup_rho[3] / tempdown[3];
			}
			else {
				pParti_VariBndy[i].sig_struct[0] = 0.0;
				pParti_VariBndy[i].sig_struct[1] = 0.0;
				pParti_VariBndy[i].sig_struct[2] = 0.0;
				pParti_VariBndy[i].sig_struct[3] = 0.0;
				pParti_VariBndy[i].sig_struct[4] = 0.0;
				pParti_VariBndy[i].sig_struct[5] = 0.0;

				pParti_VariBndy[i].vx_struct[0] = pPar[i].interfss[0];
				pParti_VariBndy[i].vx_struct[1] = pPar[i].interfss[1];
				pParti_VariBndy[i].vx_struct[2] = pPar[i].interfss[2];

				pParti_VariBndy[i].rho_struct = 0.0;
			}
		}
	}
}

//free-slip boundary by Tran et al. 2019 Computers and Geotechnics
void clBndy_Fun::free_slip_tran(Particle* pPar, Par_Cell* pParCell, clVar_Boundary* pParti_VariBndy,
	const Para_Pro& pPPro, int cn) {
	//local variables
	int ntotal = pPPro.ntotal;
	int ndim = pPPro.ndim;
	int i, j, k, nc, nj;
	double tempup_sig[4][6], tempup_rho[4], tempdown[4];
	double satu, vect_1[3], vect_2[4][3], dst, inprod;
	double dr = pPPro.dr;

	//codes
#pragma omp parallel for schedule(static) private(j, nc, nj, k, tempup_sig, \
tempup_rho, tempdown, satu, vect_1, vect_2, dst, inprod)
	for (i = 0; i < ntotal; i++)
	{
		//initialization
		for (k = 0; k < 6; k++) {
			tempup_sig[0][k] = 0.0;
			tempup_sig[1][k] = 0.0;
			tempup_sig[2][k] = 0.0;
			tempup_sig[3][k] = 0.0;
		}

		tempup_rho[0] = 0.0;
		tempup_rho[1] = 0.0;
		tempup_rho[2] = 0.0;
		tempup_rho[3] = 0.0;

		tempdown[0] = 0.0;
		tempdown[1] = 0.0;
		tempdown[2] = 0.0;
		tempdown[3] = 0.0;

		if (pPar[i].type == 0) {

			vect_2[0][0] = 0.0;
			vect_2[0][1] = 0.0;
			vect_2[0][2] = 0.0;
			vect_2[1][0] = 0.0;
			vect_2[1][1] = 0.0;
			vect_2[1][2] = 0.0;
			vect_2[2][0] = 0.0;
			vect_2[2][1] = 0.0;
			vect_2[2][2] = 0.0;
			vect_2[3][0] = 0.0;
			vect_2[3][1] = 0.0;
			vect_2[3][2] = 0.0;

			//calculating tempup_sig[4][6], tempup_vel[4][3], tempup_rho[4], and tempdown[4];
			nc = pParCell[i].ninflu;
			for (nj = 1; nj <= nc; nj++) {
				j = pParCell[i].influ[nj];
				//vector
				vect_1[0] = pPar[i].xp[0] - pPar[j].xp[0];
				vect_1[1] = pPar[i].xp[1] - pPar[j].xp[1];
				vect_1[2] = pPar[i].xp[2] - pPar[j].xp[2];

				dst = sqrt(vect_1[0] * vect_1[0] + vect_1[1] * vect_1[1] + vect_1[2] * vect_1[2]);
				inprod = pPar[j].vx[0] * vect_1[0] + pPar[j].vx[1] * vect_1[1] + pPar[j].vx[2] * vect_1[2];
				//water
				if (pPar[j].type == 1) {
					tempup_sig[0][0] = tempup_sig[0][0] + pPar[j].mass * pPar[j].sig[0] * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempup_sig[0][1] = tempup_sig[0][1] + pPar[j].mass * pPar[j].sig[1] * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempup_sig[0][2] = tempup_sig[0][2] + pPar[j].mass * pPar[j].sig[2] * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempup_sig[0][3] = tempup_sig[0][3] + pPar[j].mass * pPar[j].sig[3] * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempup_sig[0][4] = tempup_sig[0][4] + pPar[j].mass * pPar[j].sig[4] * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempup_sig[0][5] = tempup_sig[0][5] + pPar[j].mass * pPar[j].sig[5] * pParCell[i].wij[nj][3] / pPar[j].rho;

					if(dst > 0.1 *dr) {
						vect_2[0][0] = vect_2[0][0] + 
						pPar[j].mass * (pPar[j].vx[0] - inprod * vect_1[0] / dst / dst) * pParCell[i].wij[nj][3] / pPar[j].rho;
						vect_2[0][1] = vect_2[0][1] + 
						pPar[j].mass * (pPar[j].vx[1] - inprod * vect_1[1] / dst / dst) * pParCell[i].wij[nj][3] / pPar[j].rho;
						vect_2[0][2] = vect_2[0][2] + 
						pPar[j].mass * (pPar[j].vx[2] - inprod * vect_1[2] / dst / dst) * pParCell[i].wij[nj][3] / pPar[j].rho;
					}

					tempup_rho[0] = tempup_rho[0] + pPar[j].mass * pPar[j].rho * pParCell[i].wij[nj][3] / pPar[j].rho;

					tempdown[0] = tempdown[0] + pPar[j].mass * pParCell[i].wij[nj][3] / pPar[j].rho;
				}

				//soil
				if (pPar[j].type == 2) {
					satu = pPar[j].satu;
					tempup_sig[1][0] = tempup_sig[1][0] +
						pPar[j].mass * (pPar[j].sig[0] - satu * pPar[j].prew + (1 - satu) * pPar[j].prea) * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempup_sig[1][1] = tempup_sig[1][1] +
						pPar[j].mass * (pPar[j].sig[1] - satu * pPar[j].prew + (1 - satu) * pPar[j].prea) * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempup_sig[1][2] = tempup_sig[1][2] +
						pPar[j].mass * (pPar[j].sig[2] - satu * pPar[j].prew + (1 - satu) * pPar[j].prea) * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempup_sig[1][3] = tempup_sig[1][3] + pPar[j].mass * pPar[j].sig[3] * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempup_sig[1][4] = tempup_sig[1][4] + pPar[j].mass * pPar[j].sig[4] * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempup_sig[1][5] = tempup_sig[1][5] + pPar[j].mass * pPar[j].sig[5] * pParCell[i].wij[nj][3] / pPar[j].rho;

					if(dst > 0.1 *dr) {
						vect_2[1][0] = vect_2[1][0] +
						pPar[j].mass * (pPar[j].vx[0] - inprod * vect_1[0] / dst / dst) * pParCell[i].wij[nj][3] / pPar[j].rho;
						vect_2[1][1] = vect_2[1][1] +
						pPar[j].mass * (pPar[j].vx[1] - inprod * vect_1[1] / dst / dst) * pParCell[i].wij[nj][3] / pPar[j].rho;
						vect_2[1][2] = vect_2[1][2] +
						pPar[j].mass * (pPar[j].vx[2] - inprod * vect_1[2] / dst / dst) * pParCell[i].wij[nj][3] / pPar[j].rho;
					}

					tempup_rho[1] = tempup_rho[1] + pPar[j].mass * pPar[j].rho * pParCell[i].wij[nj][3] / pPar[j].rho;

					tempdown[1] = tempdown[1] + pPar[j].mass * pParCell[i].wij[nj][3] / pPar[j].rho;
				}

				//air
				if (pPar[j].type == 3) {
					tempup_sig[2][0] = tempup_sig[2][0] + pPar[j].mass * pPar[j].sig[0] * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempup_sig[2][1] = tempup_sig[2][1] + pPar[j].mass * pPar[j].sig[1] * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempup_sig[2][2] = tempup_sig[2][2] + pPar[j].mass * pPar[j].sig[2] * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempup_sig[2][3] = tempup_sig[2][3] + pPar[j].mass * pPar[j].sig[3] * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempup_sig[2][4] = tempup_sig[2][4] + pPar[j].mass * pPar[j].sig[4] * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempup_sig[2][5] = tempup_sig[2][5] + pPar[j].mass * pPar[j].sig[5] * pParCell[i].wij[nj][3] / pPar[j].rho;

					if(dst > 0.1 *dr) {
						vect_2[2][0] = vect_2[2][0] +
						pPar[j].mass * (pPar[j].vx[0] - inprod * vect_1[0] / dst / dst) * pParCell[i].wij[nj][3] / pPar[j].rho;
						vect_2[2][1] = vect_2[2][1] +
						pPar[j].mass * (pPar[j].vx[1] - inprod * vect_1[1] / dst / dst) * pParCell[i].wij[nj][3] / pPar[j].rho;
						vect_2[2][2] = vect_2[2][2] +
						pPar[j].mass * (pPar[j].vx[2] - inprod * vect_1[2] / dst / dst) * pParCell[i].wij[nj][3] / pPar[j].rho;
					}

					tempup_rho[2] = tempup_rho[2] + pPar[j].mass * pPar[j].rho * pParCell[i].wij[nj][3] / pPar[j].rho;

					tempdown[2] = tempdown[2] + pPar[j].mass * pParCell[i].wij[nj][3] / pPar[j].rho;
				}

				//structure
				if (pPar[j].type == 4) {
					tempup_sig[3][0] = tempup_sig[3][0] + pPar[j].mass * pPar[j].sig[0] * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempup_sig[3][1] = tempup_sig[3][1] + pPar[j].mass * pPar[j].sig[1] * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempup_sig[3][2] = tempup_sig[3][2] + pPar[j].mass * pPar[j].sig[2] * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempup_sig[3][3] = tempup_sig[3][3] + pPar[j].mass * pPar[j].sig[3] * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempup_sig[3][4] = tempup_sig[3][4] + pPar[j].mass * pPar[j].sig[4] * pParCell[i].wij[nj][3] / pPar[j].rho;
					tempup_sig[3][5] = tempup_sig[3][5] + pPar[j].mass * pPar[j].sig[5] * pParCell[i].wij[nj][3] / pPar[j].rho;

					if(dst > 0.1 *dr) {
						vect_2[3][0] = vect_2[3][0] +
						pPar[j].mass * (pPar[j].vx[0] - inprod * vect_1[0] / dst / dst) * pParCell[i].wij[nj][3] / pPar[j].rho;
						vect_2[3][1] = vect_2[3][1] +
						pPar[j].mass * (pPar[j].vx[1] - inprod * vect_1[1] / dst / dst) * pParCell[i].wij[nj][3] / pPar[j].rho;
						vect_2[3][2] = vect_2[3][2] +
						pPar[j].mass * (pPar[j].vx[2] - inprod * vect_1[2] / dst / dst) * pParCell[i].wij[nj][3] / pPar[j].rho;
					}
					
					tempup_rho[3] = tempup_rho[3] + pPar[j].mass * pPar[j].rho * pParCell[i].wij[nj][3] / pPar[j].rho;

					tempdown[3] = tempdown[3] + pPar[j].mass * pParCell[i].wij[nj][3] / pPar[j].rho;
				}

			}

			//calculating stress tensor, velocity and density of boundary particles
			if (fabs(tempdown[0]) > 1.0e-6) { //water
				pParti_VariBndy[i].sig_water[0] = tempup_sig[0][0] / tempdown[0];
				pParti_VariBndy[i].sig_water[1] = tempup_sig[0][1] / tempdown[0];
				pParti_VariBndy[i].sig_water[2] = tempup_sig[0][2] / tempdown[0];
				pParti_VariBndy[i].sig_water[3] = tempup_sig[0][3] / tempdown[0];
				pParti_VariBndy[i].sig_water[4] = tempup_sig[0][4] / tempdown[0];
				pParti_VariBndy[i].sig_water[5] = tempup_sig[0][5] / tempdown[0];

				pParti_VariBndy[i].vx_water[0] = pPar[i].interfss[0] + vect_2[0][0] / tempdown[0];
				pParti_VariBndy[i].vx_water[1] = pPar[i].interfss[1] + vect_2[0][1] / tempdown[0];
				pParti_VariBndy[i].vx_water[2] = pPar[i].interfss[2] + vect_2[0][2] / tempdown[0];

				pParti_VariBndy[i].rho_water = tempup_rho[0] / tempdown[0];
			}
			else {
				pParti_VariBndy[i].sig_water[0] = 0.0;
				pParti_VariBndy[i].sig_water[1] = 0.0;
				pParti_VariBndy[i].sig_water[2] = 0.0;
				pParti_VariBndy[i].sig_water[3] = 0.0;
				pParti_VariBndy[i].sig_water[4] = 0.0;
				pParti_VariBndy[i].sig_water[5] = 0.0;

				pParti_VariBndy[i].vx_water[0] = pPar[i].interfss[0];
				pParti_VariBndy[i].vx_water[1] = pPar[i].interfss[1];
				pParti_VariBndy[i].vx_water[2] = pPar[i].interfss[2];

				pParti_VariBndy[i].rho_water = 0.0;
			}

			if (fabs(tempdown[1]) > 1.0e-6) { //soil
				pParti_VariBndy[i].sig_soil[0] = tempup_sig[1][0] / tempdown[1];
				pParti_VariBndy[i].sig_soil[1] = tempup_sig[1][1] / tempdown[1];
				pParti_VariBndy[i].sig_soil[2] = tempup_sig[1][2] / tempdown[1];
				pParti_VariBndy[i].sig_soil[3] = tempup_sig[1][3] / tempdown[1];
				pParti_VariBndy[i].sig_soil[4] = tempup_sig[1][4] / tempdown[1];
				pParti_VariBndy[i].sig_soil[5] = tempup_sig[1][5] / tempdown[1];

				pParti_VariBndy[i].vx_soil[0] = pPar[i].interfss[0] + vect_2[1][0] / tempdown[1];
				pParti_VariBndy[i].vx_soil[1] = pPar[i].interfss[1] + vect_2[1][1] / tempdown[1];
				pParti_VariBndy[i].vx_soil[2] = pPar[i].interfss[2] + vect_2[1][2] / tempdown[1];

				pParti_VariBndy[i].rho_soil = tempup_rho[1] / tempdown[1];
			}
			else {
				pParti_VariBndy[i].sig_soil[0] = 0.0;
				pParti_VariBndy[i].sig_soil[1] = 0.0;
				pParti_VariBndy[i].sig_soil[2] = 0.0;
				pParti_VariBndy[i].sig_soil[3] = 0.0;
				pParti_VariBndy[i].sig_soil[4] = 0.0;
				pParti_VariBndy[i].sig_soil[5] = 0.0;

				pParti_VariBndy[i].vx_soil[0] = pPar[i].interfss[0];
				pParti_VariBndy[i].vx_soil[1] = pPar[i].interfss[1];
				pParti_VariBndy[i].vx_soil[2] = pPar[i].interfss[2];

				pParti_VariBndy[i].rho_soil = 0.0;
			}

			if (fabs(tempdown[2]) > 1.0e-6) { //air
				pParti_VariBndy[i].sig_air[0] = tempup_sig[2][0] / tempdown[2];
				pParti_VariBndy[i].sig_air[1] = tempup_sig[2][1] / tempdown[2];
				pParti_VariBndy[i].sig_air[2] = tempup_sig[2][2] / tempdown[2];
				pParti_VariBndy[i].sig_air[3] = tempup_sig[2][3] / tempdown[2];
				pParti_VariBndy[i].sig_air[4] = tempup_sig[2][4] / tempdown[2];
				pParti_VariBndy[i].sig_air[5] = tempup_sig[2][5] / tempdown[2];

				pParti_VariBndy[i].vx_air[0] = pPar[i].interfss[0] + vect_2[2][0] / tempdown[2];
				pParti_VariBndy[i].vx_air[1] = pPar[i].interfss[1] + vect_2[2][1] / tempdown[2];
				pParti_VariBndy[i].vx_air[2] = pPar[i].interfss[2] + vect_2[2][2] / tempdown[2];

				pParti_VariBndy[i].rho_air = tempup_rho[2] / tempdown[2];
			}
			else {
				pParti_VariBndy[i].sig_air[0] = 0.0;
				pParti_VariBndy[i].sig_air[1] = 0.0;
				pParti_VariBndy[i].sig_air[2] = 0.0;
				pParti_VariBndy[i].sig_air[3] = 0.0;
				pParti_VariBndy[i].sig_air[4] = 0.0;
				pParti_VariBndy[i].sig_air[5] = 0.0;

				pParti_VariBndy[i].vx_air[0] = pPar[i].interfss[0];
				pParti_VariBndy[i].vx_air[1] = pPar[i].interfss[1];
				pParti_VariBndy[i].vx_air[2] = pPar[i].interfss[2];

				pParti_VariBndy[i].rho_air = 0.0;
			}

			if (fabs(tempdown[3]) > 1.0e-6) { //structure
				pParti_VariBndy[i].sig_struct[0] = tempup_sig[3][0] / tempdown[3];
				pParti_VariBndy[i].sig_struct[1] = tempup_sig[3][1] / tempdown[3];
				pParti_VariBndy[i].sig_struct[2] = tempup_sig[3][2] / tempdown[3];
				pParti_VariBndy[i].sig_struct[3] = tempup_sig[3][3] / tempdown[3];
				pParti_VariBndy[i].sig_struct[4] = tempup_sig[3][4] / tempdown[3];
				pParti_VariBndy[i].sig_struct[5] = tempup_sig[3][5] / tempdown[3];

				pParti_VariBndy[i].vx_struct[0] = pPar[i].interfss[0] + vect_2[3][0] / tempdown[3];
				pParti_VariBndy[i].vx_struct[1] = pPar[i].interfss[1] + vect_2[3][1] / tempdown[3];
				pParti_VariBndy[i].vx_struct[2] = pPar[i].interfss[2] + vect_2[3][2] / tempdown[3];

				pParti_VariBndy[i].rho_struct = tempup_rho[3] / tempdown[3];
			}
			else {
				pParti_VariBndy[i].sig_struct[0] = 0.0;
				pParti_VariBndy[i].sig_struct[1] = 0.0;
				pParti_VariBndy[i].sig_struct[2] = 0.0;
				pParti_VariBndy[i].sig_struct[3] = 0.0;
				pParti_VariBndy[i].sig_struct[4] = 0.0;
				pParti_VariBndy[i].sig_struct[5] = 0.0;

				pParti_VariBndy[i].vx_struct[0] = pPar[i].interfss[0];
				pParti_VariBndy[i].vx_struct[1] = pPar[i].interfss[1];
				pParti_VariBndy[i].vx_struct[2] = pPar[i].interfss[2];

				pParti_VariBndy[i].rho_struct = 0.0;
			}
		}
	}
}


//boundary treatment using the interaction forces (normal force and tangential frction force)
void clBndy_Fun::boundary_free_solid(Particle* pPar, Par_Cell* pParCell, const Para_Pro& pPPro, Para_SSInter pSSInter, int cn) {

	int i, j, k, nc, nj, err;
	int id1, id2, id3;
	double dst, mag_fn, mag_ft;
	clPoint p1, p2, p3, p4, outp, sectp;
	Vector3 v1, v2, v3, v4;
	clLine line1;
	clLine3D line31, line32;
	clPlane plane1;
	double temp_ms, temp_vs[3];
	bool checker;

	//intial setting
	double dt = pPPro.dt;
	double dr = pPPro.dr;
	int ntotal = pPPro.ntotal;
	int ndim = pPPro.ndim;
	double myu = pSSInter.myu;
	double zeta = pSSInter.zeta;

	//variables
	clPair_SS pPair_SS;
	clInterf_SS pInterf_SS;
	clParti_Pair pPartiPair_SS;

	//2D problem
	if (ndim == 2) {
#pragma omp parallel for schedule(static) private(j, nc, nj, k, dst, id1, id2, p1, p2, p4, outp, \
v1, v2, v3, line1, mag_fn, mag_ft, checker, pPair_SS, pInterf_SS, pPartiPair_SS)
		for (i = 0; i < ntotal; i++)
		{
			if (pPar[i].type == 1 || pPar[i].type == 2 || pPar[i].type == 3 || pPar[i].type == 4) {
				//variable initialization
				pPartiPair_SS.info_reset();
				pPair_SS.info_reset();
				pInterf_SS.info_reset();
				//searching contacting boundary particles and sorting by distance
				k = 0;
				nc = pParCell[i].ninflu;
				for (nj = 1; nj <= nc; nj++)
				{
					j = pParCell[i].influ[nj];
					if (pPar[j].type == 0)
					{
						pPartiPair_SS.parti_id[k] = j;
						dst = (pPar[i].xp[0] - pPar[j].xp[0]) * (pPar[i].xp[0] - pPar[j].xp[0]) + (pPar[i].xp[1] - pPar[j].xp[1]) * (pPar[i].xp[1] - pPar[j].xp[1]) + (pPar[i].xp[2] - pPar[j].xp[2]) * (pPar[i].xp[2] - pPar[j].xp[2]);
						dst = sqrt(dst);
						pPartiPair_SS.parti_dst[k] = dst;
						k += 1;
					}
				}
				pPartiPair_SS.total = k;
				pPartiPair_SS.dist_sorting();
				//two dimensional: calculating contacting forces between soil and boundary
				if (pPartiPair_SS.total >= 2) {
					converse_point(pPar[i], &p4);

					//converting particle information to point and vector
					pPair_SS.p1id = pPartiPair_SS.parti_id[0];
					id1 = pPartiPair_SS.parti_id[0];
					converse_point(pPar[id1], &p1);
					converse_vector(pPar[id1], &v1);

					pPair_SS.p2id = pPartiPair_SS.parti_id[1];
					id2 = pPartiPair_SS.parti_id[1];
					converse_point(pPar[id2], &p2);
					converse_vector(pPar[id2], &v2);

					//calculating line between p1 and p2
					cal_line(p1, p2, &line1);

					//calculating distance between point and line
					cal_dpline(p4, line1, &pPair_SS.dp);

					//calculating reflecting point on the line
					cal_mappoint_line(p4, line1, &outp);

					//check if the mapping point is between p1 and p2
					v3.x = outp.x;
					v3.y = outp.y;
					v3.z = 0.0;

					checker = check_ifinside_line(v1, v2, v3);

					if (checker && pPair_SS.dp < 0.5 * dr)
					{
						//calculating the unit normal  vector of line
						line1.set_normal();
						if (check_ifoutward2D(p4, outp, line1) == false)
						{
							line1.n_x = -line1.n_x;
							line1.n_y = -line1.n_y;
						}

						//interpolate velocity and mass on the mapping point
						linear_interpolate(p1, p2, outp, pPar[id1].mass, pPar[id2].mass,
							pPar[id1].interfss, pPar[id2].interfss, &pPair_SS.ms, pPair_SS.vs);
						pPair_SS.vs[0] -= pPar[i].vxp[0];
						pPair_SS.vs[1] -= pPar[i].vxp[1];
						decomposition_velocity2D(&pPair_SS, line1);

						//calculating the contacting force
						pInterf_SS.f_n[0] = pPar[i].mass * pPair_SS.vn[0] / dt;
						pInterf_SS.f_n[1] = pPar[i].mass * pPair_SS.vn[1] / dt;

						pInterf_SS.f_t[0] = pPar[i].mass * pPair_SS.vt[0] / dt;
						pInterf_SS.f_t[1] = pPar[i].mass * pPair_SS.vt[1] / dt;

						//telling if the frictional force exceeds u*fn and correction
						mag_fn = sqrt(pInterf_SS.f_n[0] * pInterf_SS.f_n[0] +
							pInterf_SS.f_n[1] * pInterf_SS.f_n[1]);

						mag_ft = sqrt(pInterf_SS.f_t[0] * pInterf_SS.f_t[0] +
							pInterf_SS.f_t[1] * pInterf_SS.f_t[1]);

						if (mag_ft > (mag_fn * myu))
						{
							pPair_SS.vt[0] = pPair_SS.vt[0] * myu * mag_fn / mag_ft;
							pPair_SS.vt[1] = pPair_SS.vt[1] * myu * mag_fn / mag_ft;

							pInterf_SS.f_t[0] = pPar[i].mass * pPair_SS.vt[0] / dt;
							pInterf_SS.f_t[1] = pPar[i].mass * pPair_SS.vt[1] / dt;
						}

						//update the acceleration
						pPar[i].ax[0] += zeta * (pInterf_SS.f_t[0] + pInterf_SS.f_n[0]) / pPar[i].mass;
						pPar[i].ax[1] += zeta * (pInterf_SS.f_t[1] + pInterf_SS.f_n[1]) / pPar[i].mass;
						pPar[i].ax[2] += 0.0;
					}
					else
						continue;
				}
				else
					continue;
			}
		}
	}
	//3D problem
	else if (ndim == 3) {
#pragma omp parallel for schedule(static) private(j, nc, nj, k, dst, id1, id2, id3, p1, p2, p3, p4, outp, \
sectp, v1, v2, v3, v4, line31, line32, plane1, temp_ms, temp_vs, mag_fn, mag_ft, checker, err, pPair_SS, pInterf_SS, pPartiPair_SS)
		for (i = 0; i < ntotal; i++)
		{
			//variable initialization
			pPartiPair_SS.info_reset();
			pPair_SS.info_reset();
			pInterf_SS.info_reset();
			//searching contacting boundary particles and sorting by distance
			k = 0;
			nc = pParCell[i].ninflu;
			for (nj = 1; nj <= nc; nj++)
			{
				j = pParCell[i].influ[nj];
				if (pPar[j].type == 0)
				{
					pPartiPair_SS.parti_id[k] = j;
					dst = (pPar[i].xp[0] - pPar[j].xp[0]) * (pPar[i].xp[0] - pPar[j].xp[0]) + (pPar[i].xp[1] - pPar[j].xp[1]) * (pPar[i].xp[1] - pPar[j].xp[1]) + (pPar[i].xp[2] - pPar[j].xp[2]) * (pPar[i].xp[2] - pPar[j].xp[2]);
					dst = sqrt(dst);
					pPartiPair_SS.parti_dst[k] = dst;
					k += 1;
				}
			}
			pPartiPair_SS.total = k;
			pPartiPair_SS.dist_sorting();
			//three dimensional: calculating contacting forces between soil and boundary
			if (pPartiPair_SS.total >= 3) {
				converse_point(pPar[i], &p4);
				//converting particle information to point and vector
				pPair_SS.p1id = pPartiPair_SS.parti_id[0];
				id1 = pPartiPair_SS.parti_id[0];
				converse_point(pPar[id1], &p1);
				converse_vector(pPar[id1], &v1);

				pPair_SS.p2id = pPartiPair_SS.parti_id[1];
				id2 = pPartiPair_SS.parti_id[1];
				converse_point(pPar[id2], &p2);
				converse_vector(pPar[id2], &v2);

				pPair_SS.p3id = pPartiPair_SS.parti_id[2];
				id3 = pPartiPair_SS.parti_id[2];
				converse_point(pPar[id3], &p3);
				converse_vector(pPar[id3], &v3);

				//calculating plane between p1, p2, and p3
				cal_plane(p1, p2, p3, &plane1);

				//calculating distance between point and plane
				cal_dpplane(p4, plane1, &pPair_SS.dp);

				//calculating reflecting point on the line
				cal_mappoint_plane(p4, plane1, &outp);

				//check if the mapping point is between p1 and p2
				v4.x = outp.x;
				v4.y = outp.y;
				v4.z = outp.z;

				checker = check_ifinside_triangle(v1, v2, v3, v4);

				if (checker && pPair_SS.dp < 0.5 * dr)
				{
					//calculating the unit normal  vector of line
					plane1.set_normal();
					if (check_ifoutward3D(p4, outp, plane1) == false)
					{
						plane1.n_x = -plane1.n_x;
						plane1.n_y = -plane1.n_y;
						plane1.n_z = -plane1.n_z;
					}

					//calculating line1 between p1 and p2
					cal_line3D(p1, p2, &line31);
					//calculating line2 between p3 and mapping point
					cal_line3D(p3, outp, &line32);

					//calculating intersecting pbetween line1 and line2
					err = cal_intersect_lines3D(line31, line32, &sectp);
					if (err != 0)
						continue;

					//interpolate velocity and mass on line of p1 and p2
					linear_interpolate(p1, p2, sectp, pPar[id1].mass, pPar[id2].mass,
						pPar[id1].interfss, pPar[id2].interfss, &temp_ms, temp_vs);

					//interpolate velocity and mass on line of p3 and mapping point
					linear_interpolate(p3, sectp, outp, pPar[id3].mass, temp_ms,
						pPar[id3].interfss, temp_vs, &pPair_SS.ms, pPair_SS.vs);
					pPair_SS.vs[0] -= pPar[i].vxp[0];
					pPair_SS.vs[1] -= pPar[i].vxp[1];
					pPair_SS.vs[2] -= pPar[i].vxp[2];

					//velocity decomposing
					decomposition_velocity3D(&pPair_SS, plane1);

					//calculating the contacting force
					pInterf_SS.f_n[0] = pPar[i].mass * pPair_SS.vn[0] / dt;
					pInterf_SS.f_n[1] = pPar[i].mass * pPair_SS.vn[1] / dt;
					pInterf_SS.f_n[2] = pPar[i].mass * pPair_SS.vn[2] / dt;

					pInterf_SS.f_t[0] = pPar[i].mass * pPair_SS.vt[0] / dt;
					pInterf_SS.f_t[1] = pPar[i].mass * pPair_SS.vt[1] / dt;
					pInterf_SS.f_t[2] = pPar[i].mass * pPair_SS.vt[2] / dt;

					//telling if the frictional force exceeds u*fn and correction
					mag_fn = sqrt(pInterf_SS.f_n[0] * pInterf_SS.f_n[0] +
						pInterf_SS.f_n[1] * pInterf_SS.f_n[1] +
						pInterf_SS.f_n[2] * pInterf_SS.f_n[2]);

					mag_ft = sqrt(pInterf_SS.f_t[0] * pInterf_SS.f_t[0] +
						pInterf_SS.f_t[1] * pInterf_SS.f_t[1] +
						pInterf_SS.f_t[2] * pInterf_SS.f_t[2]);

					if (mag_ft > (mag_fn * myu))
					{
						pPair_SS.vt[0] = pPair_SS.vt[0] * myu * mag_fn / mag_ft;
						pPair_SS.vt[1] = pPair_SS.vt[1] * myu * mag_fn / mag_ft;
						pPair_SS.vt[2] = pPair_SS.vt[2] * myu * mag_fn / mag_ft;

						pInterf_SS.f_t[0] = pPar[i].mass * pPair_SS.vt[0] / dt;
						pInterf_SS.f_t[1] = pPar[i].mass * pPair_SS.vt[1] / dt;
						pInterf_SS.f_t[2] = pPar[i].mass * pPair_SS.vt[2] / dt;
					}

					//update the acceleration
					pPar[i].ax[0] += zeta * (pInterf_SS.f_t[0] + pInterf_SS.f_n[0]) / pPar[i].mass;
					pPar[i].ax[1] += zeta * (pInterf_SS.f_t[1] + pInterf_SS.f_n[1]) / pPar[i].mass;
					pPar[i].ax[2] += zeta * (pInterf_SS.f_t[2] + pInterf_SS.f_n[2]) / pPar[i].mass;
				}
				else
					continue;
			}
			else
				continue;
		}
	}
}

//calculating a line by two points in 3D space
void clBndy_Fun::cal_line3D(clPoint p1, clPoint p2, clLine3D* outLine)
{
	outLine->x0 = p1.x;
	outLine->y0 = p1.y;
	outLine->z0 = p1.z;

	double dst = 0.0;
	dst = (p2.x - p1.x) * (p2.x - p1.x) + (p2.y - p1.y) * (p2.y - p1.y) + (p2.z - p1.z) * (p2.z - p1.z);
	dst = sqrt(dst);

	outLine->a = (p2.x - p1.x) / dst;
	outLine->b = (p2.y - p1.y) / dst;
	outLine->c = (p2.z - p1.z) / dst;
}

//calculating the interacting point for 3D lines
int clBndy_Fun::cal_intersect_lines3D(clLine3D line1, clLine3D line2, clPoint* outPoint)
{
	double t, s;
	double flag, tempup;

	flag = line1.a * line2.b - line2.a * line1.b;

	if (fabs(flag) > ep_err)
	{
		tempup = line2.a * (line1.y0 - line2.y0) - line2.b * (line1.x0 - line2.x0);
		t = tempup / flag;
		tempup = line1.a * (line1.y0 - line2.y0) - line1.b * (line1.x0 - line2.x0);
		s = tempup / flag;

		outPoint->x = line1.x0 + line1.a * t;
		outPoint->y = line1.y0 + line1.b * t;
		outPoint->z = line1.z0 + line1.c * t;
	}
	else
		return 1;

	return 0;
}

//Check if a vector is outward for 2D case
bool clBndy_Fun::check_ifoutward2D(clPoint base_p, clPoint reflect_p, clLine inLine)
{
	Vector3 v1, v2;

	v1.x = reflect_p.x - base_p.x;
	v1.y = reflect_p.y - base_p.y;
	v1.z = 0.0;

	v2.x = inLine.n_x;
	v2.y = inLine.n_y;
	v2.z = 0.0;

	double dot = v1 * v2;
	if (dot <= 0.0)
	{
		return true;
	}
	else
		return false;
}

//Check if a vector is outward for 3D case
bool clBndy_Fun::check_ifoutward3D(clPoint base_p, clPoint reflect_p, clPlane inPlane)
{
	Vector3 v1, v2;

	v1.x = reflect_p.x - base_p.x;
	v1.y = reflect_p.y - base_p.y;
	v1.z = reflect_p.z - base_p.z;

	v2.x = inPlane.n_x;
	v2.y = inPlane.n_y;
	v2.z = inPlane.n_z;

	double dot = v1 * v2;
	if (dot <= 0.0)
	{
		return true;
	}
	else
		return false;
}

//check if a point is in a line
bool clBndy_Fun::check_ifinside_line(Vector3 A, Vector3 B, Vector3 P)
{
	Vector3 v0 = B - A;
	Vector3 v1 = P - A;

	double dot00 = v0 * v0;
	double dot01 = v0 * v1;

	double inver = 1.0 / dot00;

	double u = inver * dot01;

	if (u < 0 || u > 1) // if u out of range, return directly
	{
		return false;
	}
	else
		return true;
}

//check if a point is in a triangle
bool clBndy_Fun::check_ifinside_triangle(Vector3 A, Vector3 B, Vector3 C, Vector3 P)
{
	Vector3 v0 = C - A;
	Vector3 v1 = B - A;
	Vector3 v2 = P - A;

	double dot00 = v0 * v0;
	double dot01 = v0 * v1;
	double dot02 = v0 * v2;
	double dot11 = v1 * v1;
	double dot12 = v1 * v2;

	double inverDeno = 1 / (dot00 * dot11 - dot01 * dot01);

	double u = (dot11 * dot02 - dot01 * dot12) * inverDeno;
	if (u < 0 || u > 1) // if u out of range, return directly
	{
		return false;
	}

	double v = (dot00 * dot12 - dot01 * dot02) * inverDeno;
	if (v < 0 || v > 1) // if v out of range, return directly
	{
		return false;
	}
	return u + v <= 1;
}

//calcualting the referencing mass and velocity of particle P in literature
void clBndy_Fun::linear_interpolate(clPoint P1, clPoint P2, clPoint ref_p,
	double m1, double m2, double v1[3], double v2[3], double* ref_m, double* ref_v)
{
	double xk1, xk2, xp;

	//calculating distance
	xk1 = 0.0;
	xk2 = sqrt((P2.x - P1.x) * (P2.x - P1.x) + (P2.y - P1.y) * (P2.y - P1.y) + (P2.z - P1.z) * (P2.z - P1.z));
	xp = sqrt((ref_p.x - P1.x) * (ref_p.x - P1.x) + (ref_p.y - P1.y) * (ref_p.y - P1.y) + (ref_p.z - P1.z) * (ref_p.z - P1.z));

	*ref_m = (m1 * fabs(xp - xk2) + m2 * fabs(xp - xk1)) / fabs(xk2 - xk1);
	ref_v[0] = (v1[0] * fabs(xp - xk2) + v2[0] * fabs(xp - xk1)) / fabs(xk2 - xk1);
	ref_v[1] = (v1[1] * fabs(xp - xk2) + v2[1] * fabs(xp - xk1)) / fabs(xk2 - xk1);
	ref_v[2] = (v1[2] * fabs(xp - xk2) + v2[2] * fabs(xp - xk1)) / fabs(xk2 - xk1);
}

//interpolate soil structure force on the structure particles
void clBndy_Fun::linear_interpolate_force(clPoint P1, clPoint P2, clPoint ref_p,
	double ss_force[3], double* outforce1, double* outforce2)
{
	double xk1, xk2, xp;

	//calculating distance
	xk1 = 0.0;
	xk2 = sqrt((P2.x - P1.x) * (P2.x - P1.x) + (P2.y - P1.y) * (P2.y - P1.y) + (P2.z - P1.z) * (P2.z - P1.z));
	xp = sqrt((ref_p.x - P1.x) * (ref_p.x - P1.x) + (ref_p.y - P1.y) * (ref_p.y - P1.y) + (ref_p.z - P1.z) * (ref_p.z - P1.z));

	outforce1[0] = -fabs(xp - xk2) / fabs(xk1 - xk2) * ss_force[0];
	outforce1[1] = -fabs(xp - xk2) / fabs(xk1 - xk2) * ss_force[1];
	outforce1[2] = -fabs(xp - xk2) / fabs(xk1 - xk2) * ss_force[2];

	outforce2[0] = -fabs(xp - xk1) / fabs(xk1 - xk2) * ss_force[0];
	outforce2[1] = -fabs(xp - xk1) / fabs(xk1 - xk2) * ss_force[1];
	outforce2[2] = -fabs(xp - xk1) / fabs(xk1 - xk2) * ss_force[2];
}

//calculating the mapping point on line
void clBndy_Fun::cal_mappoint_line(clPoint ref_point, clLine Line, clPoint* OutPoint)
{
	double d_A, d_B; //direction vector
	clPoint P1, P2;	 //extended points
	double dp, temp; //relative distance

	//setting direction vector
	d_A = Line.A;
	d_B = Line.B;

	//calculating relative distance
	cal_dpline(ref_point, Line, &dp);

	//extend points
	temp = sqrt(d_A * d_A + d_B * d_B);

	P1.x = ref_point.x + dp * d_A / temp;
	P1.y = ref_point.y + dp * d_B / temp;

	P2.x = ref_point.x - dp * d_A / temp;
	P2.y = ref_point.y - dp * d_B / temp;

	//telling which of P1 and P2 is on the line
	cal_dpline(P1, Line, &dp);
	if (dp < ep_err)
	{
		*OutPoint = P1;
	}
	else
		*OutPoint = P2;
}

//calculating the mapping point on plane
void clBndy_Fun::cal_mappoint_plane(clPoint ref_point, clPlane Plane, clPoint* OutPoint)
{
	double d_A, d_B, d_C; //direction vector
	clPoint P1, P2;		  //extended points
	double dp, temp;	  //relative distance

	//setting direction vector
	d_A = Plane.A;
	d_B = Plane.B;
	d_C = Plane.C;

	//calculating relative distance
	cal_dpplane(ref_point, Plane, &dp);

	//extend points
	temp = sqrt(d_A * d_A + d_B * d_B + d_C * d_C);
	P1.x = ref_point.x + dp * d_A / temp;
	P1.y = ref_point.y + dp * d_B / temp;
	P1.z = ref_point.z + dp * d_C / temp;

	P2.x = ref_point.x - dp * d_A / temp;
	P2.y = ref_point.y - dp * d_B / temp;
	P2.z = ref_point.z - dp * d_C / temp;

	//telling which of P1 and P2 is on the plane
	cal_dpplane(P1, Plane, &dp);
	if (dp < ep_err)
	{
		*OutPoint = P1;
	}
	else
		*OutPoint = P2;
}

//decompose velocity
void clBndy_Fun::decomposition_velocity2D(clPair_SS* pPair_SS, clLine inLine)
{
	double v;
	v = pPair_SS->vs[0] * inLine.n_x + pPair_SS->vs[1] * inLine.n_y;
	pPair_SS->vn[0] = v * inLine.n_x;
	pPair_SS->vn[1] = v * inLine.n_y;

	pPair_SS->vt[0] = pPair_SS->vs[0] - pPair_SS->vn[0];
	pPair_SS->vt[1] = pPair_SS->vs[1] - pPair_SS->vn[1];
}

//decompose velocity
void clBndy_Fun::decomposition_velocity3D(clPair_SS* pPair_SS, clPlane inPlane)
{
	double v;
	v = pPair_SS->vs[0] * inPlane.n_x + pPair_SS->vs[1] * inPlane.n_y + pPair_SS->vs[2] * inPlane.n_z;

	pPair_SS->vn[0] = v * inPlane.n_x;
	pPair_SS->vn[1] = v * inPlane.n_y;
	pPair_SS->vn[2] = v * inPlane.n_z;

	pPair_SS->vt[0] = pPair_SS->vs[0] - pPair_SS->vn[0];
	pPair_SS->vt[1] = pPair_SS->vs[1] - pPair_SS->vn[1];
	pPair_SS->vt[2] = pPair_SS->vs[2] - pPair_SS->vn[2];
}