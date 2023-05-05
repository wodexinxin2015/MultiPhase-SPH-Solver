/************************************************************************************
Multi-Phase and Parallelized SPH program
--parallelized and optimized for SMP (Symmetric multiprocessing ) architecture.
--Copyright (c) 2018-2023, Weijie ZHANG, GeoHohai, Hohai University.
************************************************************************************/

#include <math.h>
#include <omp.h>
#include "Vector3.h"
#include "Header_Parameters.h"
#include "Header_Option.h"
#include "Class_Functions.h"

//calculating the interaction force for soil and structure frictional interaction
//J. Wang and D. Chan 2014
//rigid structure
int clSoilStruct_Fun::interaction_ss_rigid(Particle *pPar, Par_Cell *pParCell,
										   clPair_SS *pPair_SS, clInterf_SS **pInterf_SS, clParti_Pair *pPartiPair_SS,
										   const Para_SSInter &pSSInter, const Para_Pro &pPPro, int cn)
{

	int i, j, k, nc, nj, err;
	int id1, id2, id3, tid;
	double dst, mag_fn, mag_ft;
	clPoint p1, p2, p3, p4, outp, sectp;
	Vector3 v1, v2, v3, v4;
	clLine line1;
	clLine3D line31, line32;
	clPlane plane1;
	double temp_ms, temp_vs[3];
	double force1_n[3], force1_t[3], force2_n[3], force2_t[3],
		force3_n[3], force3_t[3];
	bool checker;

	//intial setting
	double dt = pPPro.dt;
	double dr = pPPro.dr;
	int ntotal = pPPro.ntotal;
	int ndim = pPPro.ndim;
	double myu = pSSInter.myu;
	double zeta = pSSInter.zeta;

	//calculating relative distance between cumputing particle and structure particle
#pragma omp parallel for schedule(static) private(j, nc, nj, k, dst)
	for (i = 0; i < ntotal; i++)
	{
		//initializing of particle information
		pPartiPair_SS[i].info_reset();
		pPair_SS[i].info_reset();
		for (j = 0; j < cn; j++)
			pInterf_SS[i][j].info_reset();
		k = 0;
		//searching contacting structure particles
		nc = pParCell[i].ninflu;
		for (nj = 1; nj <= nc; nj++)
		{
			j = pParCell[i].influ[nj];
			if ((pPar[i].type == 1 || pPar[i].type == 2) && pPar[j].type == 4)
			{
				pPartiPair_SS[i].parti_id[k] = j;
				dst = (pPar[i].xp[0] - pPar[j].xp[0]) * (pPar[i].xp[0] - pPar[j].xp[0]) + (pPar[i].xp[1] - pPar[j].xp[1]) * (pPar[i].xp[1] - pPar[j].xp[1]) + (pPar[i].xp[2] - pPar[j].xp[2]) * (pPar[i].xp[2] - pPar[j].xp[2]);
				dst = sqrt(dst);
				pPartiPair_SS[i].parti_dst[k] = dst;
				k += 1;
			}
		}
		pPartiPair_SS[i].total = k;
		pPartiPair_SS[i].dist_sorting();
	}

	//2 dimensional: calculating contacting forces between soil and structure
	if (ndim == 2)
	{
#pragma omp parallel for schedule(static) private(id1, id2, tid, p1, p2, p4, outp, \
																v1, v2, v3, line1, mag_fn, mag_ft, force1_n, force1_t, force2_n, force2_t, checker)
		for (i = 0; i < ntotal; i++)
		{
			tid = omp_get_thread_num();
			if ((pPar[i].type == 1 || pPar[i].type == 2) && pPartiPair_SS[i].total >= 2)
			{
				converse_point(pPar[i], &p4);

				//converting particle information to point and vector
				pPair_SS[i].p1id = pPartiPair_SS[i].parti_id[0];
				id1 = pPartiPair_SS[i].parti_id[0];
				converse_point(pPar[id1], &p1);
				converse_vector(pPar[id1], &v1);

				pPair_SS[i].p2id = pPartiPair_SS[i].parti_id[1];
				id2 = pPartiPair_SS[i].parti_id[1];
				converse_point(pPar[id2], &p2);
				converse_vector(pPar[id2], &v2);

				//calculating line between p1 and p2
				cal_line(p1, p2, &line1);

				//calculating distance between point and line
				cal_dpline(p4, line1, &pPair_SS[i].dp);

				//calculating reflecting point on the line
				cal_mappoint_line(p4, line1, &outp);

				//check if the mapping point is between p1 and p2
				v3.x = outp.x;
				v3.y = outp.y;
				v3.z = 0.0;

				checker = check_ifinside_line(v1, v2, v3);

				if (checker && pPair_SS[i].dp < 0.5 * dr)
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
									   pPar[id1].vxp, pPar[id2].vxp, &pPair_SS[i].ms, pPair_SS[i].vs);
					pPair_SS[i].vs[0] -= pPar[i].vxp[0];
					pPair_SS[i].vs[1] -= pPar[i].vxp[1];
					decomposition_velocity2D(&pPair_SS[i], line1);

					//calculating the contacting force
					pInterf_SS[i][tid].f_n[0] = pPar[i].mass * pPair_SS[i].vn[0] / dt;
					pInterf_SS[i][tid].f_n[1] = pPar[i].mass * pPair_SS[i].vn[1] / dt;

					pInterf_SS[i][tid].f_t[0] = pPar[i].mass * pPair_SS[i].vt[0] / dt;
					pInterf_SS[i][tid].f_t[1] = pPar[i].mass * pPair_SS[i].vt[1] / dt;

					//telling if the frictional force exceeds u*fn and correction
					mag_fn = sqrt(pInterf_SS[i][tid].f_n[0] * pInterf_SS[i][tid].f_n[0] +
								  pInterf_SS[i][tid].f_n[1] * pInterf_SS[i][tid].f_n[1]);

					mag_ft = sqrt(pInterf_SS[i][tid].f_t[0] * pInterf_SS[i][tid].f_t[0] +
								  pInterf_SS[i][tid].f_t[1] * pInterf_SS[i][tid].f_t[1]);

					if (mag_ft > (mag_fn * myu))
					{
						pPair_SS[i].vt[0] = pPair_SS[i].vt[0] * myu * mag_fn / mag_ft;
						pPair_SS[i].vt[1] = pPair_SS[i].vt[1] * myu * mag_fn / mag_ft;

						pInterf_SS[i][tid].f_t[0] = pPar[i].mass * pPair_SS[i].vt[0] / dt;
						pInterf_SS[i][tid].f_t[1] = pPar[i].mass * pPair_SS[i].vt[1] / dt;
					}

					//calculating contacting force on the structure particles
					linear_interpolate_force(p1, p2, outp, pInterf_SS[i][tid].f_n,
											 force1_n, force2_n);
					linear_interpolate_force(p1, p2, outp, pInterf_SS[i][tid].f_t,
											 force1_t, force2_t);

					pInterf_SS[id1][tid].f_n[0] += force1_n[0];
					pInterf_SS[id1][tid].f_n[1] += force1_n[1];
					pInterf_SS[id1][tid].f_n[2] += force1_n[2];

					pInterf_SS[id2][tid].f_n[0] += force2_n[0];
					pInterf_SS[id2][tid].f_n[1] += force2_n[1];
					pInterf_SS[id2][tid].f_n[2] += force2_n[2];

					pInterf_SS[id1][tid].f_t[0] += force1_t[0];
					pInterf_SS[id1][tid].f_t[1] += force1_t[1];
					pInterf_SS[id1][tid].f_t[2] += force1_t[2];

					pInterf_SS[id2][tid].f_t[0] += force2_t[0];
					pInterf_SS[id2][tid].f_t[1] += force2_t[1];
					pInterf_SS[id2][tid].f_t[2] += force2_t[2];
				}
				else
					continue;
			}
		}
	}
	//3 dimensional: calculating contacting forces between soil and structure
	else if (ndim == 3)
	{
#pragma omp parallel for schedule(static) private(id1, id2, id3, tid, p1, p2, p3, p4, outp, \
																sectp, v1, v2, v3, v4, line31, line32, plane1, temp_ms, temp_vs, mag_fn, mag_ft, force1_n, force1_t, force2_n, force2_t, force3_n, force3_t, checker, err)
		for (i = 0; i < ntotal; i++)
		{
			tid = omp_get_thread_num();
			if ((pPar[i].type == 1 || pPar[i].type == 2) && pPartiPair_SS[i].total >= 3)
			{
				converse_point(pPar[i], &p4);
				//converting particle information to point and vector
				pPair_SS[i].p1id = pPartiPair_SS[i].parti_id[0];
				id1 = pPartiPair_SS[i].parti_id[0];
				converse_point(pPar[id1], &p1);
				converse_vector(pPar[id1], &v1);

				pPair_SS[i].p2id = pPartiPair_SS[i].parti_id[1];
				id2 = pPartiPair_SS[i].parti_id[1];
				converse_point(pPar[id2], &p2);
				converse_vector(pPar[id2], &v2);

				pPair_SS[i].p3id = pPartiPair_SS[i].parti_id[2];
				id3 = pPartiPair_SS[i].parti_id[2];
				converse_point(pPar[id3], &p3);
				converse_vector(pPar[id3], &v3);

				//calculating plane between p1, p2, and p3
				cal_plane(p1, p2, p3, &plane1);

				//calculating distance between point and plane
				cal_dpplane(p4, plane1, &pPair_SS[i].dp);

				//calculating reflecting point on the line
				cal_mappoint_plane(p4, plane1, &outp);

				//check if the mapping point is between p1 and p2
				v4.x = outp.x;
				v4.y = outp.y;
				v4.z = outp.z;

				checker = check_ifinside_triangle(v1, v2, v3, v4);

				if (checker && pPair_SS[i].dp < 0.5 * dr)
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
									   pPar[id1].vxp, pPar[id2].vxp, &temp_ms, temp_vs);

					//interpolate velocity and mass on line of p3 and mapping point
					linear_interpolate(p3, sectp, outp, pPar[id3].mass, temp_ms,
									   pPar[id3].vxp, temp_vs, &pPair_SS[i].ms, pPair_SS[i].vs);
					pPair_SS[i].vs[0] -= pPar[i].vxp[0];
					pPair_SS[i].vs[1] -= pPar[i].vxp[1];
					pPair_SS[i].vs[2] -= pPar[i].vxp[2];

					//velocity decomposing
					decomposition_velocity3D(&pPair_SS[i], plane1);

					//calculating the contacting force
					pInterf_SS[i][tid].f_n[0] = pPar[i].mass * pPair_SS[i].vn[0] / dt;
					pInterf_SS[i][tid].f_n[1] = pPar[i].mass * pPair_SS[i].vn[1] / dt;
					pInterf_SS[i][tid].f_n[2] = pPar[i].mass * pPair_SS[i].vn[2] / dt;

					pInterf_SS[i][tid].f_t[0] = pPar[i].mass * pPair_SS[i].vt[0] / dt;
					pInterf_SS[i][tid].f_t[1] = pPar[i].mass * pPair_SS[i].vt[1] / dt;
					pInterf_SS[i][tid].f_t[2] = pPar[i].mass * pPair_SS[i].vt[2] / dt;

					//telling if the frictional force exceeds u*fn and correction
					mag_fn = sqrt(pInterf_SS[i][tid].f_n[0] * pInterf_SS[i][tid].f_n[0] +
								  pInterf_SS[i][tid].f_n[1] * pInterf_SS[i][tid].f_n[1] +
								  pInterf_SS[i][tid].f_n[2] * pInterf_SS[i][tid].f_n[2]);

					mag_ft = sqrt(pInterf_SS[i][tid].f_t[0] * pInterf_SS[i][tid].f_t[0] +
								  pInterf_SS[i][tid].f_t[1] * pInterf_SS[i][tid].f_t[1] +
								  pInterf_SS[i][tid].f_t[2] * pInterf_SS[i][tid].f_t[2]);

					if (mag_ft > (mag_fn * myu))
					{
						pPair_SS[i].vt[0] = pPair_SS[i].vt[0] * myu * mag_fn / mag_ft;
						pPair_SS[i].vt[1] = pPair_SS[i].vt[1] * myu * mag_fn / mag_ft;
						pPair_SS[i].vt[2] = pPair_SS[i].vt[2] * myu * mag_fn / mag_ft;

						pInterf_SS[i][tid].f_t[0] = pPar[i].mass * pPair_SS[i].vt[0] / dt;
						pInterf_SS[i][tid].f_t[1] = pPar[i].mass * pPair_SS[i].vt[1] / dt;
						pInterf_SS[i][tid].f_t[2] = pPar[i].mass * pPair_SS[i].vt[2] / dt;
					}

					//calculating contacting force on the structure particles
					linear_interpolate_force(p3, sectp, outp, pInterf_SS[i][tid].f_n,
											 force1_n, force3_n);
					linear_interpolate_force(p3, sectp, outp, pInterf_SS[i][tid].f_t,
											 force1_t, force3_t);

					pInterf_SS[id3][tid].f_n[0] += force1_n[0];
					pInterf_SS[id3][tid].f_n[1] += force1_n[1];
					pInterf_SS[id3][tid].f_n[2] += force1_n[2];

					pInterf_SS[id3][tid].f_t[0] += force1_t[0];
					pInterf_SS[id3][tid].f_t[1] += force1_t[1];
					pInterf_SS[id3][tid].f_t[2] += force1_t[2];

					linear_interpolate_force(p1, p2, sectp, force3_n,
											 force1_n, force2_n);
					linear_interpolate_force(p1, p2, sectp, force3_t,
											 force1_t, force2_t);

					pInterf_SS[id1][tid].f_n[0] -= force1_n[0];
					pInterf_SS[id1][tid].f_n[1] -= force1_n[1];
					pInterf_SS[id1][tid].f_n[2] -= force1_n[2];

					pInterf_SS[id2][tid].f_n[0] -= force2_n[0];
					pInterf_SS[id2][tid].f_n[1] -= force2_n[1];
					pInterf_SS[id2][tid].f_n[2] -= force2_n[2];

					pInterf_SS[id1][tid].f_t[0] -= force1_t[0];
					pInterf_SS[id1][tid].f_t[1] -= force1_t[1];
					pInterf_SS[id1][tid].f_t[2] -= force1_t[2];

					pInterf_SS[id2][tid].f_t[0] -= force2_t[0];
					pInterf_SS[id2][tid].f_t[1] -= force2_t[1];
					pInterf_SS[id2][tid].f_t[2] -= force2_t[2];
				}
				else
					continue;
			}
		}
	}

//summing contacting force of different threads into the particle information
#pragma omp parallel for schedule(static) private(j)
	for (i = 0; i < ntotal; i++)
	{
		if (pPar[i].type == 1 || pPar[i].type == 2 || pPar[i].type == 4)
		{
			pPar[i].interfss[0] = 0.0;
			pPar[i].interfss[1] = 0.0;
			pPar[i].interfss[2] = 0.0;

			for (j = 0; j < cn; j++)
			{
				pPar[i].interfss[0] += zeta * (pInterf_SS[i][j].f_t[0] + pInterf_SS[i][j].f_n[0]);
				pPar[i].interfss[1] += zeta * (pInterf_SS[i][j].f_t[1] + pInterf_SS[i][j].f_n[1]);
				pPar[i].interfss[2] += zeta * (pInterf_SS[i][j].f_t[2] + pInterf_SS[i][j].f_n[2]);
			}

			pPar[i].ax[0] += pPar[i].interfss[0] / pPar[i].mass;
			pPar[i].ax[1] += pPar[i].interfss[1] / pPar[i].mass;
			pPar[i].ax[2] += pPar[i].interfss[2] / pPar[i].mass;
		}
	}

	return 0;
}

//deformable structure
int clSoilStruct_Fun::interaction_ss_deformable(Particle *pPar, Par_Cell *pParCell,
												clPair_SS *pPair_SS, clInterf_SS **pInterf_SS, clParti_Pair *pPartiPair_SS,
												const Para_SSInter &pSSInter, const Para_Pro &pPPro, int cn)
{

	int i, j, k, nc, nj, err;
	int id1, id2, id3, tid;
	double dst, mag_fn, mag_ft;
	clPoint p1, p2, p3, p4, outp, sectp;
	Vector3 v1, v2, v3, v4;
	clLine line1;
	clLine3D line31, line32;
	clPlane plane1;
	double temp_ms, temp_vs[3];
	double force1_n[3], force1_t[3], force2_n[3], force2_t[3],
		force3_n[3], force3_t[3];
	bool checker;

	//intial setting
	double dt = pPPro.dt;
	double dr = pPPro.dr;
	int ntotal = pPPro.ntotal;
	int ndim = pPPro.ndim;
	double myu = pSSInter.myu;
	double zeta = pSSInter.zeta;

	//calculating relative distance between cumputing particle and structure particle
#pragma omp parallel for schedule(static) private(j, nc, nj, k, dst)
	for (i = 0; i < ntotal; i++)
	{
		//initializing of particle information
		pPartiPair_SS[i].info_reset();
		pPair_SS[i].info_reset();
		for (j = 0; j < cn; j++)
			pInterf_SS[i][j].info_reset();
		k = 0;
		//searching contacting structure particles
		nc = pParCell[i].ninflu;
		for (nj = 1; nj <= nc; nj++)
		{
			j = pParCell[i].influ[nj];
			if ((pPar[i].type == 1 || pPar[i].type == 2) && pPar[j].type == 4)
			{
				pPartiPair_SS[i].parti_id[k] = j;
				dst = (pPar[i].xp[0] - pPar[j].xp[0]) * (pPar[i].xp[0] - pPar[j].xp[0]) + (pPar[i].xp[1] - pPar[j].xp[1]) * (pPar[i].xp[1] - pPar[j].xp[1]) + (pPar[i].xp[2] - pPar[j].xp[2]) * (pPar[i].xp[2] - pPar[j].xp[2]);
				dst = sqrt(dst);
				pPartiPair_SS[i].parti_dst[k] = dst;
				k += 1;
			}
		}
		pPartiPair_SS[i].total = k;
		pPartiPair_SS[i].dist_sorting();
	}

	//2 dimensional: calculating contacting forces between soil and structure
	if (ndim == 2)
	{
#pragma omp parallel for schedule(static) private(id1, id2, tid, p1, p2, p4, outp, \
																v1, v2, v3, line1, mag_fn, mag_ft, force1_n, force1_t, force2_n, force2_t, checker)
		for (i = 0; i < ntotal; i++)
		{
			tid = omp_get_thread_num();
			if ((pPar[i].type == 1 || pPar[i].type == 2) && pPartiPair_SS[i].total >= 2)
			{
				converse_point(pPar[i], &p4);

				//converting particle information to point and vector
				pPair_SS[i].p1id = pPartiPair_SS[i].parti_id[0];
				id1 = pPartiPair_SS[i].parti_id[0];
				converse_point(pPar[id1], &p1);
				converse_vector(pPar[id1], &v1);

				pPair_SS[i].p2id = pPartiPair_SS[i].parti_id[1];
				id2 = pPartiPair_SS[i].parti_id[1];
				converse_point(pPar[id2], &p2);
				converse_vector(pPar[id2], &v2);

				//calculating line between p1 and p2
				cal_line(p1, p2, &line1);

				//calculating distance between point and line
				cal_dpline(p4, line1, &pPair_SS[i].dp);

				//calculating reflecting point on the line
				cal_mappoint_line(p4, line1, &outp);

				//check if the mapping point is between p1 and p2
				v3.x = outp.x;
				v3.y = outp.y;
				v3.z = 0.0;

				checker = check_ifinside_line(v1, v2, v3);

				if (checker && pPair_SS[i].dp < dr)
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
									   pPar[id1].vxp, pPar[id2].vxp, &pPair_SS[i].ms, pPair_SS[i].vs);

					pPair_SS[i].vs[0] -= pPar[i].vxp[0];
					pPair_SS[i].vs[1] -= pPar[i].vxp[1];
					pPair_SS[i].vs[0] *= pPair_SS[i].ms / (pPair_SS[i].ms + pPar[i].mass);
					pPair_SS[i].vs[1] *= pPair_SS[i].ms / (pPair_SS[i].ms + pPar[i].mass);
					decomposition_velocity2D(&pPair_SS[i], line1);

					//calculating the contacting force
					pInterf_SS[i][tid].f_n[0] = pPar[i].mass * pPair_SS[i].vn[0] / dt;
					pInterf_SS[i][tid].f_n[1] = pPar[i].mass * pPair_SS[i].vn[1] / dt;

					pInterf_SS[i][tid].f_t[0] = pPar[i].mass * pPair_SS[i].vt[0] / dt;
					pInterf_SS[i][tid].f_t[1] = pPar[i].mass * pPair_SS[i].vt[1] / dt;

					//telling if the frictional force exceeds u*fn and correction
					mag_fn = sqrt(pInterf_SS[i][tid].f_n[0] * pInterf_SS[i][tid].f_n[0] +
								  pInterf_SS[i][tid].f_n[1] * pInterf_SS[i][tid].f_n[1]);

					mag_ft = sqrt(pInterf_SS[i][tid].f_t[0] * pInterf_SS[i][tid].f_t[0] +
								  pInterf_SS[i][tid].f_t[1] * pInterf_SS[i][tid].f_t[1]);

					if (mag_ft > (mag_fn * myu))
					{
						pPair_SS[i].vt[0] = pPair_SS[i].vt[0] * myu * mag_fn / mag_ft;
						pPair_SS[i].vt[1] = pPair_SS[i].vt[1] * myu * mag_fn / mag_ft;

						pInterf_SS[i][tid].f_t[0] = pPar[i].mass * pPair_SS[i].vt[0] / dt;
						pInterf_SS[i][tid].f_t[1] = pPar[i].mass * pPair_SS[i].vt[1] / dt;
					}

					//calculating contacting force on the structure particles
					linear_interpolate_force(p1, p2, outp, pInterf_SS[i][tid].f_n,
											 force1_n, force2_n);
					linear_interpolate_force(p1, p2, outp, pInterf_SS[i][tid].f_t,
											 force1_t, force2_t);

					pInterf_SS[id1][tid].f_n[0] += force1_n[0];
					pInterf_SS[id1][tid].f_n[1] += force1_n[1];
					pInterf_SS[id1][tid].f_n[2] += force1_n[2];

					pInterf_SS[id2][tid].f_n[0] += force2_n[0];
					pInterf_SS[id2][tid].f_n[1] += force2_n[1];
					pInterf_SS[id2][tid].f_n[2] += force2_n[2];

					pInterf_SS[id1][tid].f_t[0] += force1_t[0];
					pInterf_SS[id1][tid].f_t[1] += force1_t[1];
					pInterf_SS[id1][tid].f_t[2] += force1_t[2];

					pInterf_SS[id2][tid].f_t[0] += force2_t[0];
					pInterf_SS[id2][tid].f_t[1] += force2_t[1];
					pInterf_SS[id2][tid].f_t[2] += force2_t[2];
				}
				else
					continue;
			}
		}
	}
	//3 dimensional: calculating contacting forces between soil and structure
	else if (ndim == 3)
	{
#pragma omp parallel for schedule(static) private(id1, id2, id3, tid, p1, p2, p3, p4, outp, \
																sectp, v1, v2, v3, v4, line31, line32, plane1, temp_ms, temp_vs, mag_fn, mag_ft, force1_n, force1_t, force2_n, force2_t, force3_n, force3_t, checker, err)
		for (i = 0; i < ntotal; i++)
		{
			tid = omp_get_thread_num();
			if ((pPar[i].type == 1 || pPar[i].type == 2) && pPartiPair_SS[i].total >= 3)
			{
				converse_point(pPar[i], &p4);
				//converting particle information to point and vector
				pPair_SS[i].p1id = pPartiPair_SS[i].parti_id[0];
				id1 = pPartiPair_SS[i].parti_id[0];
				converse_point(pPar[id1], &p1);
				converse_vector(pPar[id1], &v1);

				pPair_SS[i].p2id = pPartiPair_SS[i].parti_id[1];
				id2 = pPartiPair_SS[i].parti_id[1];
				converse_point(pPar[id2], &p2);
				converse_vector(pPar[id2], &v2);

				pPair_SS[i].p3id = pPartiPair_SS[i].parti_id[2];
				id3 = pPartiPair_SS[i].parti_id[2];
				converse_point(pPar[id3], &p3);
				converse_vector(pPar[id3], &v3);

				//calculating plane between p1, p2, and p3
				cal_plane(p1, p2, p3, &plane1);

				//calculating distance between point and plane
				cal_dpplane(p4, plane1, &pPair_SS[i].dp);

				//calculating reflecting point on the line
				cal_mappoint_plane(p4, plane1, &outp);

				//check if the mapping point is between p1 and p2
				v4.x = outp.x;
				v4.y = outp.y;
				v4.z = outp.z;

				checker = check_ifinside_triangle(v1, v2, v3, v4);

				if (checker && pPair_SS[i].dp < dr)
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
									   pPar[id1].vxp, pPar[id2].vxp, &temp_ms, temp_vs);

					//interpolate velocity and mass on line of p3 and mapping point
					linear_interpolate(p3, sectp, outp, pPar[id3].mass, temp_ms,
									   pPar[id3].vxp, temp_vs, &pPair_SS[i].ms, pPair_SS[i].vs);

					pPair_SS[i].vs[0] -= pPar[i].vxp[0];
					pPair_SS[i].vs[0] *= pPair_SS[i].ms / (pPair_SS[i].ms + pPar[i].mass);
					pPair_SS[i].vs[1] -= pPar[i].vxp[1];
					pPair_SS[i].vs[1] *= pPair_SS[i].ms / (pPair_SS[i].ms + pPar[i].mass);
					pPair_SS[i].vs[2] -= pPar[i].vxp[2];
					pPair_SS[i].vs[2] *= pPair_SS[i].ms / (pPair_SS[i].ms + pPar[i].mass);

					//velocity decomposing
					decomposition_velocity3D(&pPair_SS[i], plane1);

					//calculating the contacting force
					pInterf_SS[i][tid].f_n[0] = pPar[i].mass * pPair_SS[i].vn[0] / dt;
					pInterf_SS[i][tid].f_n[1] = pPar[i].mass * pPair_SS[i].vn[1] / dt;
					pInterf_SS[i][tid].f_n[2] = pPar[i].mass * pPair_SS[i].vn[2] / dt;

					pInterf_SS[i][tid].f_t[0] = pPar[i].mass * pPair_SS[i].vt[0] / dt;
					pInterf_SS[i][tid].f_t[1] = pPar[i].mass * pPair_SS[i].vt[1] / dt;
					pInterf_SS[i][tid].f_t[2] = pPar[i].mass * pPair_SS[i].vt[2] / dt;

					//telling if the frictional force exceeds u*fn and correction
					mag_fn = sqrt(pInterf_SS[i][tid].f_n[0] * pInterf_SS[i][tid].f_n[0] +
								  pInterf_SS[i][tid].f_n[1] * pInterf_SS[i][tid].f_n[1] +
								  pInterf_SS[i][tid].f_n[2] * pInterf_SS[i][tid].f_n[2]);

					mag_ft = sqrt(pInterf_SS[i][tid].f_t[0] * pInterf_SS[i][tid].f_t[0] +
								  pInterf_SS[i][tid].f_t[1] * pInterf_SS[i][tid].f_t[1] +
								  pInterf_SS[i][tid].f_t[2] * pInterf_SS[i][tid].f_t[2]);

					if (mag_ft > (mag_fn * myu))
					{
						pPair_SS[i].vt[0] = pPair_SS[i].vt[0] * myu * mag_fn / mag_ft;
						pPair_SS[i].vt[1] = pPair_SS[i].vt[1] * myu * mag_fn / mag_ft;
						pPair_SS[i].vt[2] = pPair_SS[i].vt[2] * myu * mag_fn / mag_ft;

						pInterf_SS[i][tid].f_t[0] = pPar[i].mass * pPair_SS[i].vt[0] / dt;
						pInterf_SS[i][tid].f_t[1] = pPar[i].mass * pPair_SS[i].vt[1] / dt;
						pInterf_SS[i][tid].f_t[2] = pPar[i].mass * pPair_SS[i].vt[2] / dt;
					}

					//calculating contacting force on the structure particles
					linear_interpolate_force(p3, sectp, outp, pInterf_SS[i][tid].f_n,
											 force1_n, force3_n);
					linear_interpolate_force(p3, sectp, outp, pInterf_SS[i][tid].f_t,
											 force1_t, force3_t);

					pInterf_SS[id3][tid].f_n[0] += force1_n[0];
					pInterf_SS[id3][tid].f_n[1] += force1_n[1];
					pInterf_SS[id3][tid].f_n[2] += force1_n[2];

					pInterf_SS[id3][tid].f_t[0] += force1_t[0];
					pInterf_SS[id3][tid].f_t[1] += force1_t[1];
					pInterf_SS[id3][tid].f_t[2] += force1_t[2];

					linear_interpolate_force(p1, p2, sectp, force3_n,
											 force1_n, force2_n);
					linear_interpolate_force(p1, p2, sectp, force3_t,
											 force1_t, force2_t);

					pInterf_SS[id1][tid].f_n[0] -= force1_n[0];
					pInterf_SS[id1][tid].f_n[1] -= force1_n[1];
					pInterf_SS[id1][tid].f_n[2] -= force1_n[2];

					pInterf_SS[id2][tid].f_n[0] -= force2_n[0];
					pInterf_SS[id2][tid].f_n[1] -= force2_n[1];
					pInterf_SS[id2][tid].f_n[2] -= force2_n[2];

					pInterf_SS[id1][tid].f_t[0] -= force1_t[0];
					pInterf_SS[id1][tid].f_t[1] -= force1_t[1];
					pInterf_SS[id1][tid].f_t[2] -= force1_t[2];

					pInterf_SS[id2][tid].f_t[0] -= force2_t[0];
					pInterf_SS[id2][tid].f_t[1] -= force2_t[1];
					pInterf_SS[id2][tid].f_t[2] -= force2_t[2];
				}
				else
					continue;
			}
		}
	}

	//summing contacting force of different threads into the particle information
#pragma omp parallel for schedule(static) private(j)
	for (i = 0; i < ntotal; i++)
	{
		if (pPar[i].type == 1 || pPar[i].type == 2 || pPar[i].type == 4)
		{
			pPar[i].interfss[0] = 0.0;
			pPar[i].interfss[1] = 0.0;
			pPar[i].interfss[2] = 0.0;

			for (j = 0; j < cn; j++)
			{
				pPar[i].interfss[0] += zeta * (pInterf_SS[i][j].f_t[0] + pInterf_SS[i][j].f_n[0]);
				pPar[i].interfss[1] += zeta * (pInterf_SS[i][j].f_t[1] + pInterf_SS[i][j].f_n[1]);
				pPar[i].interfss[2] += zeta * (pInterf_SS[i][j].f_t[2] + pInterf_SS[i][j].f_n[2]);
			}

			pPar[i].ax[0] += pPar[i].interfss[0] / pPar[i].mass;
			pPar[i].ax[1] += pPar[i].interfss[1] / pPar[i].mass;
			pPar[i].ax[2] += pPar[i].interfss[2] / pPar[i].mass;
		}
	}

	return 0;
}

//penalty method
int clSoilStruct_Fun::interaction_ss_penalty(Particle *pPar, Par_Cell *pParCell,
											 clPair_SS *pPair_SS, clInterf_SS **pInterf_SS, clParti_Pair *pPartiPair_SS,
											 const Para_SSInter &pSSInter, const Para_Pro &pPPro, int cn)
{

	int i, j, k, nc, nj, err;
	int id1, id2, id3, tid;
	double dst, mag_fn, mag_ft;
	clPoint p1, p2, p3, p4, outp, sectp;
	Vector3 v1, v2, v3, v4;
	clLine line1;
	clLine3D line31, line32;
	clPlane plane1;
	double temp_ms, temp_vs[3];
	double force1_n[3], force1_t[3], force2_n[3], force2_t[3],
		force3_n[3], force3_t[3];
	bool checker;

	//intial setting
	double dt = pPPro.dt;
	double dr = pPPro.dr;
	int ntotal = pPPro.ntotal;
	int ndim = pPPro.ndim;
	double myu = pSSInter.myu;
	double zeta = pSSInter.zeta;

	//calculating relative distance between cumputing particle and structure particle
#pragma omp parallel for schedule(static) private(j, nc, nj, k, dst)
	for (i = 0; i < ntotal; i++)
	{
		//initializing of particle information
		pPartiPair_SS[i].info_reset();
		pPair_SS[i].info_reset();
		for (j = 0; j < cn; j++)
			pInterf_SS[i][j].info_reset();
		k = 0;
		//searching contacting structure particles
		nc = pParCell[i].ninflu;
		for (nj = 1; nj <= nc; nj++)
		{
			j = pParCell[i].influ[nj];
			if ((pPar[i].type == 1 || pPar[i].type == 2) && pPar[j].type == 4)
			{
				pPartiPair_SS[i].parti_id[k] = j;
				dst = (pPar[i].xp[0] - pPar[j].xp[0]) * (pPar[i].xp[0] - pPar[j].xp[0]) + (pPar[i].xp[1] - pPar[j].xp[1]) * (pPar[i].xp[1] - pPar[j].xp[1]) + (pPar[i].xp[2] - pPar[j].xp[2]) * (pPar[i].xp[2] - pPar[j].xp[2]);
				dst = sqrt(dst);
				pPartiPair_SS[i].parti_dst[k] = dst;
				k += 1;
			}
		}
		pPartiPair_SS[i].total = k;
		pPartiPair_SS[i].dist_sorting();
	}

	//2 dimensional: calculating contacting forces between soil and structure
	if (ndim == 2)
	{
#pragma omp parallel for schedule(static) private(id1, id2, tid, p1, p2, p4, outp, \
																v1, v2, v3, line1, mag_fn, mag_ft, force1_n, force1_t, force2_n, force2_t, checker)
		for (i = 0; i < ntotal; i++)
		{
			tid = omp_get_thread_num();
			if ((pPar[i].type == 1 || pPar[i].type == 2) && pPartiPair_SS[i].total >= 2)
			{
				converse_point(pPar[i], &p4);

				//converting particle information to point and vector
				pPair_SS[i].p1id = pPartiPair_SS[i].parti_id[0];
				id1 = pPartiPair_SS[i].parti_id[0];
				converse_point(pPar[id1], &p1);
				converse_vector(pPar[id1], &v1);

				pPair_SS[i].p2id = pPartiPair_SS[i].parti_id[1];
				id2 = pPartiPair_SS[i].parti_id[1];
				converse_point(pPar[id2], &p2);
				converse_vector(pPar[id2], &v2);

				//calculating line between p1 and p2
				cal_line(p1, p2, &line1);

				//calculating distance between point and line
				cal_dpline(p4, line1, &pPair_SS[i].dp);

				//calculating reflecting point on the line
				cal_mappoint_line(p4, line1, &outp);

				//check if the mapping point is between p1 and p2
				v3.x = outp.x;
				v3.y = outp.y;
				v3.z = 0.0;

				checker = check_ifinside_line(v1, v2, v3);

				if (checker && pPair_SS[i].dp < dr)
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
									   pPar[id1].vxp, pPar[id2].vxp, &pPair_SS[i].ms, pPair_SS[i].vs);

					pPair_SS[i].vs[0] -= pPar[i].vxp[0];
					pPair_SS[i].vs[1] -= pPar[i].vxp[1];
					pPair_SS[i].vs[0] *= pPair_SS[i].ms / (pPair_SS[i].ms + pPar[i].mass);
					pPair_SS[i].vs[1] *= pPair_SS[i].ms / (pPair_SS[i].ms + pPar[i].mass);
					decomposition_velocity2D(&pPair_SS[i], line1);

					//calculating the contacting force
					mag_fn = (1.0 - zeta) * 2.0 * pPar[i].mass * (dr - pPair_SS[i].dp) / dt / dt;
					pInterf_SS[i][tid].f_n[0] = mag_fn * line1.n_x;
					pInterf_SS[i][tid].f_n[1] = mag_fn * line1.n_y;

					pInterf_SS[i][tid].f_t[0] = pPar[i].mass * pPair_SS[i].vt[0] / dt;
					pInterf_SS[i][tid].f_t[1] = pPar[i].mass * pPair_SS[i].vt[1] / dt;

					//telling if the frictional force exceeds u*fn and correction
					mag_fn = sqrt(pInterf_SS[i][tid].f_n[0] * pInterf_SS[i][tid].f_n[0] +
								  pInterf_SS[i][tid].f_n[1] * pInterf_SS[i][tid].f_n[1]);

					mag_ft = sqrt(pInterf_SS[i][tid].f_t[0] * pInterf_SS[i][tid].f_t[0] +
								  pInterf_SS[i][tid].f_t[1] * pInterf_SS[i][tid].f_t[1]);

					if (mag_ft > (mag_fn * myu))
					{
						pPair_SS[i].vt[0] = pPair_SS[i].vt[0] * myu * mag_fn / mag_ft;
						pPair_SS[i].vt[1] = pPair_SS[i].vt[1] * myu * mag_fn / mag_ft;

						pInterf_SS[i][tid].f_t[0] = pPar[i].mass * pPair_SS[i].vt[0] / dt;
						pInterf_SS[i][tid].f_t[1] = pPar[i].mass * pPair_SS[i].vt[1] / dt;
					}

					//calculating contacting force on the structure particles
					linear_interpolate_force(p1, p2, outp, pInterf_SS[i][tid].f_n,
											 force1_n, force2_n);
					linear_interpolate_force(p1, p2, outp, pInterf_SS[i][tid].f_t,
											 force1_t, force2_t);

					pInterf_SS[id1][tid].f_n[0] += force1_n[0];
					pInterf_SS[id1][tid].f_n[1] += force1_n[1];
					pInterf_SS[id1][tid].f_n[2] += force1_n[2];

					pInterf_SS[id2][tid].f_n[0] += force2_n[0];
					pInterf_SS[id2][tid].f_n[1] += force2_n[1];
					pInterf_SS[id2][tid].f_n[2] += force2_n[2];

					pInterf_SS[id1][tid].f_t[0] += force1_t[0];
					pInterf_SS[id1][tid].f_t[1] += force1_t[1];
					pInterf_SS[id1][tid].f_t[2] += force1_t[2];

					pInterf_SS[id2][tid].f_t[0] += force2_t[0];
					pInterf_SS[id2][tid].f_t[1] += force2_t[1];
					pInterf_SS[id2][tid].f_t[2] += force2_t[2];
				}
				else
					continue;
			}
		}
	}
	//3 dimensional: calculating contacting forces between soil and structure
	else if (ndim == 3)
	{
#pragma omp parallel for schedule(static) private(id1, id2, id3, tid, p1, p2, p3, p4, outp, \
																sectp, v1, v2, v3, v4, line31, line32, plane1, temp_ms, temp_vs, mag_fn, mag_ft, force1_n, force1_t, force2_n, force2_t, force3_n, force3_t, checker, err)
		for (i = 0; i < ntotal; i++)
		{
			tid = omp_get_thread_num();
			if ((pPar[i].type == 1 || pPar[i].type == 2) && pPartiPair_SS[i].total >= 3)
			{
				converse_point(pPar[i], &p4);
				//converting particle information to point and vector
				pPair_SS[i].p1id = pPartiPair_SS[i].parti_id[0];
				id1 = pPartiPair_SS[i].parti_id[0];
				converse_point(pPar[id1], &p1);
				converse_vector(pPar[id1], &v1);

				pPair_SS[i].p2id = pPartiPair_SS[i].parti_id[1];
				id2 = pPartiPair_SS[i].parti_id[1];
				converse_point(pPar[id2], &p2);
				converse_vector(pPar[id2], &v2);

				pPair_SS[i].p3id = pPartiPair_SS[i].parti_id[2];
				id3 = pPartiPair_SS[i].parti_id[2];
				converse_point(pPar[id3], &p3);
				converse_vector(pPar[id3], &v3);

				//calculating plane between p1, p2, and p3
				cal_plane(p1, p2, p3, &plane1);

				//calculating distance between point and plane
				cal_dpplane(p4, plane1, &pPair_SS[i].dp);

				//calculating reflecting point on the line
				cal_mappoint_plane(p4, plane1, &outp);

				//check if the mapping point is between p1 and p2
				v4.x = outp.x;
				v4.y = outp.y;
				v4.z = outp.z;

				checker = check_ifinside_triangle(v1, v2, v3, v4);

				if (checker && pPair_SS[i].dp < dr)
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
									   pPar[id1].vxp, pPar[id2].vxp, &temp_ms, temp_vs);

					//interpolate velocity and mass on line of p3 and mapping point
					linear_interpolate(p3, sectp, outp, pPar[id3].mass, temp_ms,
									   pPar[id3].vxp, temp_vs, &pPair_SS[i].ms, pPair_SS[i].vs);

					pPair_SS[i].vs[0] -= pPar[i].vxp[0];
					pPair_SS[i].vs[0] *= pPair_SS[i].ms / (pPair_SS[i].ms + pPar[i].mass);
					pPair_SS[i].vs[1] -= pPar[i].vxp[1];
					pPair_SS[i].vs[1] *= pPair_SS[i].ms / (pPair_SS[i].ms + pPar[i].mass);
					pPair_SS[i].vs[2] -= pPar[i].vxp[2];
					pPair_SS[i].vs[2] *= pPair_SS[i].ms / (pPair_SS[i].ms + pPar[i].mass);

					//velocity decomposing
					decomposition_velocity3D(&pPair_SS[i], plane1);

					//calculating the contacting force
					mag_fn = (1.0 - zeta) * 2.0 * pPar[i].mass * (dr - pPair_SS[i].dp) / dt / dt;

					pInterf_SS[i][tid].f_n[0] = mag_fn * plane1.n_x;
					pInterf_SS[i][tid].f_n[1] = mag_fn * plane1.n_y;
					pInterf_SS[i][tid].f_n[2] = mag_fn * plane1.n_z;

					pInterf_SS[i][tid].f_t[0] = pPar[i].mass * pPair_SS[i].vt[0] / dt;
					pInterf_SS[i][tid].f_t[1] = pPar[i].mass * pPair_SS[i].vt[1] / dt;
					pInterf_SS[i][tid].f_t[2] = pPar[i].mass * pPair_SS[i].vt[2] / dt;

					//telling if the frictional force exceeds u*fn and correction
					mag_fn = sqrt(pInterf_SS[i][tid].f_n[0] * pInterf_SS[i][tid].f_n[0] +
								  pInterf_SS[i][tid].f_n[1] * pInterf_SS[i][tid].f_n[1] +
								  pInterf_SS[i][tid].f_n[2] * pInterf_SS[i][tid].f_n[2]);

					mag_ft = sqrt(pInterf_SS[i][tid].f_t[0] * pInterf_SS[i][tid].f_t[0] +
								  pInterf_SS[i][tid].f_t[1] * pInterf_SS[i][tid].f_t[1] +
								  pInterf_SS[i][tid].f_t[2] * pInterf_SS[i][tid].f_t[2]);

					if (mag_ft > (mag_fn * myu))
					{
						pPair_SS[i].vt[0] = pPair_SS[i].vt[0] * myu * mag_fn / mag_ft;
						pPair_SS[i].vt[1] = pPair_SS[i].vt[1] * myu * mag_fn / mag_ft;
						pPair_SS[i].vt[2] = pPair_SS[i].vt[2] * myu * mag_fn / mag_ft;

						pInterf_SS[i][tid].f_t[0] = pPar[i].mass * pPair_SS[i].vt[0] / dt;
						pInterf_SS[i][tid].f_t[1] = pPar[i].mass * pPair_SS[i].vt[1] / dt;
						pInterf_SS[i][tid].f_t[2] = pPar[i].mass * pPair_SS[i].vt[2] / dt;
					}

					//calculating contacting force on the structure particles
					linear_interpolate_force(p3, sectp, outp, pInterf_SS[i][tid].f_n,
											 force1_n, force3_n);
					linear_interpolate_force(p3, sectp, outp, pInterf_SS[i][tid].f_t,
											 force1_t, force3_t);

					pInterf_SS[id3][tid].f_n[0] += force1_n[0];
					pInterf_SS[id3][tid].f_n[1] += force1_n[1];
					pInterf_SS[id3][tid].f_n[2] += force1_n[2];

					pInterf_SS[id3][tid].f_t[0] += force1_t[0];
					pInterf_SS[id3][tid].f_t[1] += force1_t[1];
					pInterf_SS[id3][tid].f_t[2] += force1_t[2];

					linear_interpolate_force(p1, p2, sectp, force3_n,
											 force1_n, force2_n);
					linear_interpolate_force(p1, p2, sectp, force3_t,
											 force1_t, force2_t);

					pInterf_SS[id1][tid].f_n[0] -= force1_n[0];
					pInterf_SS[id1][tid].f_n[1] -= force1_n[1];
					pInterf_SS[id1][tid].f_n[2] -= force1_n[2];

					pInterf_SS[id2][tid].f_n[0] -= force2_n[0];
					pInterf_SS[id2][tid].f_n[1] -= force2_n[1];
					pInterf_SS[id2][tid].f_n[2] -= force2_n[2];

					pInterf_SS[id1][tid].f_t[0] -= force1_t[0];
					pInterf_SS[id1][tid].f_t[1] -= force1_t[1];
					pInterf_SS[id1][tid].f_t[2] -= force1_t[2];

					pInterf_SS[id2][tid].f_t[0] -= force2_t[0];
					pInterf_SS[id2][tid].f_t[1] -= force2_t[1];
					pInterf_SS[id2][tid].f_t[2] -= force2_t[2];
				}
				else
					continue;
			}
		}
	}

	//summing contacting force of different threads into the particle information
#pragma omp parallel for schedule(static) private(j)
	for (i = 0; i < ntotal; i++)
	{
		if (pPar[i].type == 1 || pPar[i].type == 2 || pPar[i].type == 4)
		{
			pPar[i].interfss[0] = 0.0;
			pPar[i].interfss[1] = 0.0;
			pPar[i].interfss[2] = 0.0;

			for (j = 0; j < cn; j++)
			{
				pPar[i].interfss[0] += pInterf_SS[i][j].f_t[0] + pInterf_SS[i][j].f_n[0];
				pPar[i].interfss[1] += pInterf_SS[i][j].f_t[1] + pInterf_SS[i][j].f_n[1];
				pPar[i].interfss[2] += pInterf_SS[i][j].f_t[2] + pInterf_SS[i][j].f_n[2];
			}

			pPar[i].ax[0] += pPar[i].interfss[0] / pPar[i].mass;
			pPar[i].ax[1] += pPar[i].interfss[1] / pPar[i].mass;
			pPar[i].ax[2] += pPar[i].interfss[2] / pPar[i].mass;
		}
	}

	return 0;
}

//converse the particle information to clPoint type
void clSoilStruct_Fun::converse_point(Particle inPar, clPoint *outPoint)
{
	outPoint->x = inPar.xp[0];
	outPoint->y = inPar.xp[1];
	outPoint->z = inPar.xp[2];
}

//converse the particle information to Vector3 type
void clSoilStruct_Fun::converse_vector(Particle inPar, Vector3 *outVect)
{
	outVect->x = inPar.xp[0];
	outVect->y = inPar.xp[1];
	outVect->z = inPar.xp[2];
}

//calculating line parameters Ax+By+C=0
void clSoilStruct_Fun::cal_line(clPoint p1, clPoint p2, clLine *outLine)
{
	double dst = p2.y * p1.x - p1.y * p2.x;
	outLine->A = (p1.y - p2.y) / dst;
	outLine->B = (p2.x - p1.x) / dst;
}

//calculating plane parameters Ax+By+Cz+D=0
void clSoilStruct_Fun::cal_plane(clPoint p1, clPoint p2, clPoint p3, clPlane *outPlane)
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

//calculating a line by two points in 3D space
void clSoilStruct_Fun::cal_line3D(clPoint p1, clPoint p2, clLine3D *outLine)
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

//calculating the mapping point on line
void clSoilStruct_Fun::cal_mappoint_line(clPoint ref_point, clLine Line, clPoint *OutPoint)
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
void clSoilStruct_Fun::cal_mappoint_plane(clPoint ref_point, clPlane Plane, clPoint *OutPoint)
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

//calculating the interacting point for 3D lines
int clSoilStruct_Fun::cal_intersect_lines3D(clLine3D line1, clLine3D line2, clPoint *outPoint)
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

//calcualting the referencing mass and velocity of particle P in literature
void clSoilStruct_Fun::linear_interpolate(clPoint P1, clPoint P2, clPoint ref_p,
										  double m1, double m2, double v1[3], double v2[3], double *ref_m, double *ref_v)
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
void clSoilStruct_Fun::linear_interpolate_force(clPoint P1, clPoint P2, clPoint ref_p,
												double ss_force[3], double *outforce1, double *outforce2)
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

//calculating distance of a point to a line
void clSoilStruct_Fun::cal_dpline(clPoint ref_point, clLine inLine, double *dp)
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
void clSoilStruct_Fun::cal_dpplane(clPoint ref_point, clPlane inPlane, double *dp)
{
	double temp, dst;

	dst = fabs(inPlane.A * ref_point.x + inPlane.B * ref_point.y + inPlane.C * ref_point.z + inPlane.D);
	temp = sqrt(inPlane.A * inPlane.A + inPlane.B * inPlane.B + inPlane.C * inPlane.C);

	if (temp > ep_err)
		*dp = dst / temp;
	else
		*dp = 0.0;
}

//Check if a vector is outward for 2D case
bool clSoilStruct_Fun::check_ifoutward2D(clPoint base_p, clPoint reflect_p, clLine inLine)
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
bool clSoilStruct_Fun::check_ifoutward3D(clPoint base_p, clPoint reflect_p, clPlane inPlane)
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
bool clSoilStruct_Fun::check_ifinside_line(Vector3 A, Vector3 B, Vector3 P)
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
bool clSoilStruct_Fun::check_ifinside_triangle(Vector3 A, Vector3 B, Vector3 C, Vector3 P)
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

//decompose velocity
void clSoilStruct_Fun::decomposition_velocity2D(clPair_SS *pPair_SS, clLine inLine)
{
	double v;
	v = pPair_SS->vs[0] * inLine.n_x + pPair_SS->vs[1] * inLine.n_y;
	pPair_SS->vn[0] = v * inLine.n_x;
	pPair_SS->vn[1] = v * inLine.n_y;

	pPair_SS->vt[0] = pPair_SS->vs[0] - pPair_SS->vn[0];
	pPair_SS->vt[1] = pPair_SS->vs[1] - pPair_SS->vn[1];
}

//decompose velocity
void clSoilStruct_Fun::decomposition_velocity3D(clPair_SS *pPair_SS, clPlane inPlane)
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