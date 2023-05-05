/************************************************************************************
Multi-Phase and Parallelized SPH program
--parallelized and optimized for SMP (Symmetric multiprocessing ) architecture.
--Copyright (c) 2018-2023, Weijie ZHANG, GeoHohai, Hohai University.
************************************************************************************/

#ifndef CLASS_FUNCTIONS_H_INCLUDED
#define CLASS_FUNCTIONS_H_INCLUDED

#include <stdio.h>
#include "Vector3.h"
#include "Class_NNPS.h"
#include "Class_particle.h"
#include "Class_ParticleProperty.h"
#include "Class_Problem.h"
#include "Class_spatial_variability.h"

class clFIni_Fun
{
public:
	int initial_pro(Para_FL *pFlag, Para_Pro *pPPro, Para_Rain *pPRain, Para_Inter *pPInter, Para_SSInter *pSSInter,
					Para_Boundary *pPBou, Para_Vibra *pPVib, Para_GF *pPGf, Para_EC *pPEc, double *threshold_min, double *threshold_max, FILE *flog, FILE *finp);
	int initial_mat(Para_Model *pMatter, Para_Soil *pSoil, Para_Fluid *pWater,
					Para_Fluid *pAir, Para_Structure *pStruct, const Para_FL &pFlag, Flag_SpatialVariability *pFSV,
					Para_SpatialVariability *pPSV, FILE *flog, FILE *finp);
	void input_par(Particle *pPar, Para_Pro *pPPro, Cell_Con *pCellc, Para_Rain *pPRain, FILE *flog, FILE *finp);
	void setting_vis_cons(Particle *pPar, Para_Pro &pPPro, Para_Fluid *pVFluid, Para_Fluid *pWater, Para_Fluid &pAir,
						  Para_Soil *pParti_ConsPara, Para_Soil *pSoil, Para_Structure &pStruct);
	void setting_par(Particle *pPar, Para_Pro *pPPro, Para_Fluid *pVFluid, Para_Soil *pParti_ConsPara,
					 Para_Soil *pSoil, Para_GF pPGf, const Para_Structure &pStruct, Para_Fluid *pWater, Para_Rain *pPRain, FILE *flog);
	int initial_vib(Para_Vibra *pPVib, FILE *flog, FILE *fvib);
};

class clErr_Fun
{
public:
	void Err_Deal(int err_type, FILE *flog);
	void Err_Debug(int step, int pNum, FILE *flog);
};

class clNNPS_Fun
{
public:
	void celllist_ini1(Cell_Con *pCellc, const Para_Pro &pPPro);
	void celllist_ini2(cell_link *pCell_Link, const Cell_Con &pCellc, int dim, int cn);
	int parttocell(Par_Cell *pParCell, cell_info *pCell_Info, const Cell_Con &pCellc,
						   parti_cellid *pParti_Cell_Sorted, parti_cellid *pParti_Cell_Temp, 
						   const Para_Pro &pPPro, Particle *pPar, int dim, int cn, int size_m, int k);
	int partisearching(Par_Cell *pParCell, cell_info *pCell_Info, cell_link *pCell_Link, parti_cellid *pParti_Cell,
							   const Cell_Con &pCellc, const Para_Pro &pPPro, Particle *pPar, int cn, int flknl);
	void kernel_gradient_correction(Particle *pPar, Par_Cell *pParCell, const Para_Pro &pPPro);

private:
	void inverse_mat_2D(double (*A)[3], double (*A_inv)[3]);
	void inverse_mat_3D(double (*A)[3], double (*A_inv)[3]);
	void Merge_CellID(parti_cellid *arr_A, parti_cellid *arr_temp, int start, int mid_index, int end);
	void Binary_Sort(parti_cellid *pParti_Cell_Sorted, parti_cellid *pParti_Cell_Temp,
							 cell_info *pCell_Info, int ntotal, int ctotal, int size_m, int k);
};

class clDensity_Fun
{
public:
	void density(Particle *pPar, Par_Cell *pParCell, Para_Fluid *pVFluid,
				 cl_bndy_pair *pBndy_Pair, clVar_Boundary *pParti_VariBndy,
				 const Para_Pro &pPPro, int bndytype);
};

class clStraStre_Fun
{
public:
	void strain_noreg(Particle *pPar, Par_Cell *pParCell,
					  cl_bndy_pair *pBndy_Pair, clVar_Boundary *pParti_VariBndy,
					  const Para_Pro &pPPro, int bndytype);
	void strain_to_stress(Particle *pPar, const Para_Soil *pParti_ConsPara, StiffMat *pParStiff, const Para_Fluid *pVFluid,
						  const Para_Structure &pStruct, const Para_Pro &pPPro, int cn, int flts, double alpha);
	void Stress_eigValue(Particle *pPar, const Para_Pro &pPPro, int cn);
	void eos(Particle *pPar, Par_Cell *pParCell, Para_Fluid *pVFluid, const Para_Pro &pPPro, int cn);
	void stress_regu(Particle *pPar, Par_Cell *pParCell, const Para_Pro &pPPro);
	void strain_regu(Particle *pPar, Par_Cell *pParCell, const Para_Pro &pPPro);
	void free_surface_check(Particle *pPar, Par_Cell *pParCell, const Para_Pro &pPPro, double max, double min);

private:
	void elastic(double (*dept)[6], const Particle &pStrePar, const Para_Soil &pSoil);
	void elastic_struct(double (*dept)[6], const Para_Structure &pSoil);
	void sub_camclay(double (*dept)[6], const Particle &pPar, const Para_Soil &pSoil, double *adps);
	void dpmodel1(double (*dept)[6], const Particle &pStrePar, const Para_Soil &pSoil, double *adps);
	void dpmodel_bui(double (*dept)[6], const Particle &pStrePar, const Para_Soil &pSoil, double *adps);
	void DPModel_check(Particle *pStrePar, const Para_Soil &pSoil, double dt);
	void sub_camunsatu(double (*dept)[6], double *erfs, const Particle &pPar, const Para_Soil &pSoil, double *adps);
	void sub_camanisotropy(double (*dept)[6], double *erfs, const Particle &pPar, const Para_Soil &pSoil, double *adps);
	void elastic_plastic(double (*dept)[6], const Particle &pStrePar, const Para_Soil &pSoil, double *adps);
	void MCModel_check(Particle *pStrePar, const Para_Soil &pSoil, double dt);
	void viscous_elastic(Particle *pPar, const Para_Soil &pSoil);
	void fractional_mcc(double (*dept)[6], double (*cijkl)[6], Particle &pPar, const Para_Soil &pSoil);
	void fractional_plastic(double *dsig, double *deps, double *adps, double cijkl[6][6]);

	void TSSMP_before(Particle *pPar, double alpha);
	void TSSMP_after(Particle *pPar);

	double detA(double arcs[][6], int n);
	void Display(double arr[][6], int n);
	void getAStart(double arcs[6][6], int n, double ans[6][6]);
};

class clInterFor_Fun
{
public:
	void saturation(Particle *pPar, Par_Cell *pParCell, const Para_Soil *pParti_ConsPara, const Para_Pro &pPPro, int cn);
	void cal_permeability(Particle *pPar, const Para_Inter &pPInter, const Para_Soil *pParti_ConsPara, const Para_Pro &pPPro, int cn);
	void interact1(Particle *pPar, Par_Cell *pParCell, const Para_Inter &pInter,
				   const Para_Fluid *pVFluid, const Para_Pro &pPPro, Para_GF pgf, int cn);
	void interact2(Particle *pPar, Par_Cell *pParCell, const Para_Inter &pInter,
				   const Para_Soil *pParti_ConsPara, const Para_Fluid *pVFluid, const Para_Pro &pPPro,
				   Para_GF pgf, int cn);
	void interact3(Particle *pPar, Par_Cell *pParCell, const Para_Inter &pInter,
				   const Para_Soil *pParti_ConsPara, const Para_Fluid *pVFluid, const Para_Pro &pPPro,
				   Para_GF pgf, int cn);
};

class clAcce_Fun
{
public:
	void fluid_acceleration_water(Particle *pPar, Par_Cell *pParCell, const Para_Fluid *pVFluid,
								  const Para_Pro &pPPro, const Para_GF &pPGf, int cn);
	void fluid_acceleration_air(Particle *pPar, Par_Cell *pParCell, const Para_Fluid *pVFluid,
								const Para_Pro &pPPro, const Para_GF &pPGf, int cn);
	void fluid_acceleration_minorwater(Particle *pPar, Par_Cell *pParCell, const Para_Fluid *pVFluid,
									   const Para_Pro &pPPro, const Para_GF &pPGf, int cn);
	void fluid_acceleration_minorair(Particle *pPar, Par_Cell *pParCell, const Para_Fluid *pVFluid,
									 const Para_Pro &pPPro, const Para_GF &pPGf, int cn);
	void soil_acceleration_noreg(Particle *pPar, Par_Cell *pParCell, const Para_Soil *pParti_ConsPara,
								 const Para_Pro &pPPro, const Para_GF &pPGf, int cn);
	void structure_acceleration_noreg(Particle *pPar, Par_Cell *pParCell,
									  const Para_Pro &pPPro, const Para_GF &pPGf, int cn);
	void boundary_moving_stress_effect(Particle *pPar, Par_Cell *pParCell, const Para_Pro &pPPro, int cn);
	void boundary_effect_free_struct(Particle *pPar, Par_Cell *pParCell, const Para_Pro &pPPro, int cn);
	void artificial_viscosity(Particle *pPar, Par_Cell *pParCell, const Para_Soil *pParti_ConsPara,
							  const Para_Fluid *pVFluid, const Para_Pro &pPPro, int cn);
	void artificial_stress_2d(Particle *pPar, Par_Cell *pParCell, const Para_Pro &pPPro, clRaij *Rij, int cn);
	void artificial_stress_3d(Particle *pPar, Par_Cell *pParCell, const Para_Pro &pPPro, clRaij *Rij, int cn);
	void vibration_acceleration(Particle *pPar, const Para_Vibra &pPVib, const Para_Pro &pPPro, int cn);
	void vibration_velocity(Particle *pPar, const Para_Vibra &pPVib, const Para_Pro &pPPro, int cn);
	void coupled_Water_Soil_Structure(Particle *pPar, Par_Cell *pParCell, const Para_Fluid *pVFluid, const Para_Soil *pParti_ConsPara,
									  const Para_Pro &pPPro, const Para_GF &pPGf, int cn);
	void soil_damping_geostress(Particle *pPar, Par_Cell *pParCell, const Para_Pro &pPPro, int cn);
	void boundary_stress_effect_tran(Particle *pPar, Par_Cell *pParCell, clVar_Boundary *pParti_VariBndy,
									 const Para_Pro &pPPro, int cn);
	void soil_damping_vibration(Particle *pPar, Par_Cell *pParCell, const Para_Soil *pParti_ConsPara,
								StiffMat *pParStiff, const Para_Pro &pPPro);

private:
	void inverse_mat(double (*A)[3], double (*A_inv)[3]);
};

class clVel_Fun
{
public:
	void veocity_update(Particle *pPar, const Para_Pro &pPPro, const Para_Boundary &pPBou, int flv);
	void xsph(Particle *pPar, Par_Cell *pParCell, cl_bndy_pair *pBndy_pair, clVar_Boundary *pParti_VariBndy,
			  const Para_Pro &pPPro, double (*vx_temp)[3], int bndytype);
};

class clUpdate_Fun
{
public:
	void update(Particle *pPar, const Para_Pro &pPPro, const Para_Boundary &pPBou, int cn);
	void timestep(Para_Pro *pPPro);
};

class clBndy_Fun
{
public:
	void tmboundary(Particle *pPar, Par_Cell *pParCell, cl_bndy_pair *pBndy_Pair, const Para_Pro &pPPro);
	void tm_boundary_2d(Particle *pPar, Par_Cell *pParCell, cl_bndy_pair *pBndy_Pair, clVar_Boundary *pParti_VariBndy,
						const Para_Pro &pPPro, int cn);
	void tm_boundary_3d(Particle *pPar, Par_Cell *pParCell, cl_bndy_pair *pBndy_Pair, clVar_Boundary *pParti_VariBndy, const Para_Pro &pPPro, int cn);
	void no_slip_tran(Particle *pPar, Par_Cell *pParCell, clVar_Boundary *pParti_VariBndy,
					  const Para_Pro &pPPro, int cn);
	void free_slip_tran(Particle *pPar, Par_Cell *pParCell, clVar_Boundary *pParti_VariBndy,
						const Para_Pro &pPPro, int cn);
	void check_domain(Particle *pPar, const Cell_Con &pCellc,
					  const Para_Boundary &pBou, const Para_Pro &pPPro);
	void velocity_set_domain(Particle *pPar, Par_Cell *pParCell, cl_bndy_pair *pBndy_Pair, const Para_Pro &pPPro, double coe_bndy);
	void boundary_free_solid(Particle *pPar, Par_Cell *pParCell, const Para_Pro &pPPro, Para_SSInter pSSInter, int cn);

private:
	void converse_point(Particle inPar, clPoint *outPoint);
	void converse_vector(Particle inPar, Vector3 *outVect);
	void cal_plane(clPoint p1, clPoint p2, clPoint p3, clPlane *outPlane);
	void cal_line(clPoint p1, clPoint p2, clLine *outLine);
	void cal_dpline(clPoint ref_point, clLine inLine, double *dp);
	void cal_dpplane(clPoint ref_point, clPlane inPlane, double *dp);
	void plane_two_points(clPoint p1, clPoint p2, clPlane *outPlane);
	void cal_line3D(clPoint p1, clPoint p2, clLine3D *outLine);
	int cal_intersect_lines3D(clLine3D line1, clLine3D line2, clPoint *outPoint);
	bool check_ifoutward2D(clPoint base_p, clPoint reflect_p, clLine inLine);
	bool check_ifoutward3D(clPoint base_p, clPoint reflect_p, clPlane inPlane);
	bool check_ifinside_line(Vector3 A, Vector3 B, Vector3 P);
	bool check_ifinside_triangle(Vector3 A, Vector3 B, Vector3 C, Vector3 P);
	void linear_interpolate(clPoint P1, clPoint P2, clPoint ref_p,
							double m1, double m2, double v1[3], double v2[3], double *ref_m, double *ref_v);
	void linear_interpolate_force(clPoint P1, clPoint P2, clPoint ref_p,
								  double ss_force[3], double *outforce1, double *outforce2);
	void cal_mappoint_line(clPoint ref_point, clLine Line, clPoint *OutPoint);
	void cal_mappoint_plane(clPoint ref_point, clPlane Plane, clPoint *OutPoint);
	void decomposition_velocity2D(clPair_SS *pPair_SS, clLine inLine);
	void decomposition_velocity3D(clPair_SS *pPair_SS, clPlane inPlane);
};

class clOutput_Fun
{
public:
	void mgf(Particle *pVelPar, const Para_Pro &pPPro, const Cell_Con &pCellc, FILE *fp);
	void output_total(Particle *pVelPar, const Para_Pro &pPPro, char *argv);
	void output_impact(FILE *fpo, Particle *pVelPar, const Para_Pro &pPPro, char *argv);
	void output_flowdist(FILE *fpi, Particle *pVelPar, const Para_Pro &pPPro);
	void output_status(int step, int start, int end);
};

class clExcBac_Fun
{
public:
	void exca_or_back(Particle *pPar, Para_Pro *pPPro, const Para_EC &pPEc, int cn, int nstage);
};

class clRain_Fun
{
public:
	void rainini(const Para_Rain &pPRain, const Para_Pro &pPPro, int *lp, double *vrain);
	void rainpro(Particle *pVelPar, Para_Rain *pPRain, Para_Pro *pPPro,
				 const Para_Fluid *pVFluid, Para_GF pgf, double vrain);
	void partidel(Particle *pVelPar, Para_Rain *pPRain,
				  Para_Pro *pPPro, const Cell_Con &pCellc, int cn);
	void cementpro(Particle *pVelPar, Para_Rain *pPRain, Para_Pro *pPPro,
				   const Para_Fluid *pVFluid, Para_GF pgf, double vrain, double cement_base[3]);
	void center_cement(Particle *pVelPar, Para_Pro &pPPro, double *x);
	void waterdel_min(Particle *pVelPar, Para_Rain *pPRain, Para_Pro *pPPro);
	void waterdel_max(Particle *pVelPar, Para_Rain *pPRain, Para_Pro *pPPro);
};

class clSoilStruct_Fun
{
public:
	int interaction_ss_rigid(Particle *pPar, Par_Cell *pParCell,
							 clPair_SS *pPair_SS, clInterf_SS **pInterf_SS, clParti_Pair *pPartiPair_SS,
							 const Para_SSInter &pSSInter, const Para_Pro &pPPro, int cn);
	int interaction_ss_deformable(Particle *pPar, Par_Cell *pParCell,
								  clPair_SS *pPair_SS, clInterf_SS **pInterf_SS, clParti_Pair *pPartiPair_SS,
								  const Para_SSInter &pSSInter, const Para_Pro &pPPro, int cn);
	int interaction_ss_penalty(Particle *pPar, Par_Cell *pParCell,
							   clPair_SS *pPair_SS, clInterf_SS **pInterf_SS, clParti_Pair *pPartiPair_SS,
							   const Para_SSInter &pSSInter, const Para_Pro &pPPro, int cn);

private:
	void converse_point(Particle inPar, clPoint *outPoint);
	void converse_vector(Particle inPar, Vector3 *outVect);
	void cal_line(clPoint p1, clPoint p2, clLine *outLine);
	void cal_plane(clPoint p1, clPoint p2, clPoint p3, clPlane *outPlane);
	void cal_line3D(clPoint p1, clPoint p2, clLine3D *outLine);
	void cal_mappoint_line(clPoint ref_point, clLine Line, clPoint *OutPoint);
	void cal_mappoint_plane(clPoint ref_point, clPlane Plane, clPoint *OutPoint);
	int cal_intersect_lines3D(clLine3D line1, clLine3D line2, clPoint *outPoint);
	void linear_interpolate(clPoint P1, clPoint P2, clPoint ref_p,
							double m1, double m2, double v1[3], double v2[3], double *ref_m, double *ref_v);
	void linear_interpolate_force(clPoint P1, clPoint P2, clPoint ref_p,
								  double ss_force[3], double *outforce1, double *outforce2);
	void cal_dpline(clPoint ref_point, clLine inLine, double *dp);
	void cal_dpplane(clPoint ref_point, clPlane inPlane, double *dp);
	bool check_ifoutward2D(clPoint base_p, clPoint reflect_p, clLine inLine);
	bool check_ifoutward3D(clPoint base_p, clPoint reflect_p, clPlane inPlane);
	bool check_ifinside_line(Vector3 A, Vector3 B, Vector3 P);
	bool check_ifinside_triangle(Vector3 A, Vector3 B, Vector3 C, Vector3 P);
	void decomposition_velocity2D(clPair_SS *pPair_SS, clLine inLine);
	void decomposition_velocity3D(clPair_SS *pPair_SS, clPlane inPlane);
};

class clFreeFieldBndy_Fun
{
public:
	void free_field_solid(Particle *pPar, cl_bndy_pair *pBndy_pair,
						  const Para_Soil pFFBound, Par_Cell *pParCell, const Para_Pro &pPPro, int cn);

private:
	void ff_stress_2d_solid(double *sigxx, double *sigyy, double vx_a[3], double vx_b[3],
							double e, double poi, double dens, int type);
	void ff_stress_3d_solid(double *sigxx, double *sigyy, double *sigzz, double vx_a[3], double vx_b[3],
							double e, double poi, double dens, int type);
};

class clFreeFieldBndy_Fun_soil
{
public:
	void free_field_solid(Particle *pPar, const Para_Soil *pParti_ConsPara, Par_Cell *pParCell, const Para_Pro &pPPro, int cn);

private:
	void ff_stress_2d_solid(double *sigxx, double *sigyy, double vx_a[3], double vx_b[3],
							double e, double poi, double dens, int type);
	void ff_stress_3d_solid(double *sigxx, double *sigyy, double *sigzz, double vx_a[3], double vx_b[3],
							double e, double poi, double dens, int type);
};
#endif
