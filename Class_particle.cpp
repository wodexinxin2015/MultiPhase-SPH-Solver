/************************************************************************************
Multi-Phase and Parallelized SPH program
--parallelized and optimized for SMP (Symmetric multiprocessing ) architecture.
--Copyright (c) 2018-2023, Weijie ZHANG, GeoHohai, Hohai University.
************************************************************************************/

#include <cmath>
#include "Class_particle.h"

/*initializing of mass*/

Particle::Particle()
{
	int i;

	for (i = 0; i < 3; i++)
	{
		this->x[i] = 0.0;
		this->xp[i] = 0.0;
		this->vx[i] = 0.0;
		this->vxp[i] = 0.0;
		this->ax[i] = 0.0;
		this->axp[i] = 0.0;

		this->ux[i] = 0.0; /*2016-7-5*/

		this->strep[i] = 0.0; /*2016-7-7*/

		this->divx[i][0] = 0.0;
		this->divx[i][1] = 0.0;
		this->divx[i][2] = 0.0;

		this->interf[i] = 0.0;
		this->interfss[i] = 0.0;

		this->weps[i][0] = 0.0;
		this->weps[i][1] = 0.0;
		this->weps[i][2] = 0.0;

		this->beta[i][0] = 0.0;
		this->beta[i][1] = 0.0;
		this->beta[i][2] = 0.0;

		this->beta_ts[i][0] = 0.0;
		this->beta_ts[i][1] = 0.0;
		this->beta_ts[i][2] = 0.0;
	}

	tsqc = 0.0;

	rho = 0.0;
	rhop = 0.0;
	mass = 0.0;
	hl = 0.0;

	pre = 0.0;
	prep = 0.0;

	e = 0.0;

	for (i = 0; i < 6; i++)
	{
		this->eps[i] = 0.0;
		this->deps[i] = 0.0;
		this->sig[i] = 0.0;
		this->dsig[i] = 0.0;
		this->tssig[i] = 0.0;
		this->adps[i] = 0.0;
	}
	ff_sig[0] = 0.0;
	ff_sig[1] = 0.0;
	ff_sig[2] = 0.0;

	satu = 0.0;
	dsatu = 0.0;
	suct = 0.0;
	cop = 0.0;

	vsig = 0.0;
	veps = 0.0;
	divq = 0.0;
	divr = 0.0;
	meps = 0.0;

	roue = 0.0;
	rous = 0.0;
	zeta = 0.0;
	rsta = 0.0;
	prew = 0.0;
	prea = 0.0;

	porosity = 0.0;
	cd = 0.0;
	poro_cop = 0.0;

	type = 0;
	matype = 0;
	permtype = 0;
	etype = 0;
}

Particle::~Particle()
{
}

//calculating mass for particles
void Particle::massini(int dim)
{
	if (dim == 2)
	{
		mass = rhop * hl * hl;
	}
	else if (dim == 3)
	{
		mass = rhop * hl * hl * hl;
	}
}

/*calculation of the mean stress, mean strain,
	deviator stress and deviator stress ratio*/
void Particle::deq()
{
	double sxx, syy, szz, sxy, syz, szx;
	double j2;
	/*mean strain*/
	veps = eps[0] + eps[1] + eps[2];
	/*mean stress*/
	vsig = (sig[0] + sig[1] + sig[2]) / 3;
	/*deviatoric stress and deviatoric stress ratio*/
	divq = 0.0;
	sxx = sig[0] - vsig;
	syy = sig[1] - vsig;
	szz = sig[2] - vsig;
	sxy = sig[3];
	syz = sig[4];
	szx = sig[5];
	j2 = 0.5 * (sxx * sxx + syy * syy + szz * szz) + sxy * sxy + syz * syz + szx * szx;
	divq = sqrt(3 * fabs(j2));
	if (vsig != 0.0)
		divr = fabs(divq / vsig);
	else
		divr = 0.0;
}

/*calculation of permeability according to degree of saturation*/
//variable permeability
void Particle::permeablel(double sitas, double sitar, double pers, double perl)
{
	if (satu < sitar)
		cop = pers;
	else if (satu > sitas)
		cop = perl;
	else
		cop = (perl - pers) * satu / (sitas - sitar) + pers;
}

//three layers
void Particle::permeabletl(double pers, double perl, double perm)
{
	//adding permtype
	if (permtype == 1)
	{
		cop = pers; //small permeability
	}
	if (permtype == 2)
	{
		cop = perm; // middle permeability
	}
	if (permtype == 3)
	{
		cop = perl; //large permeability
	}
}

void Particle::shearstrain()
{
	//sqrt(2*I2) I2 is the second invariant of deviatoric strain tensor

	double sxx, syy, szz, sxy, syz, szx;
	double j2;

	sxx = eps[0] - veps / 3.0;
	syy = eps[1] - veps / 3.0;
	szz = eps[2] - veps / 3.0;
	sxy = eps[3];
	syz = eps[4];
	szx = eps[5];

	j2 = 0.5 * (sxx * sxx + syy * syy + szz * szz) + sxy * sxy + syz * syz + szx * szx;

	meps = sqrt(2.0 * fabs(j2));
}

void Particle::poro_cal(double dt)
{
	double de, et;
	de = (1 + e) * (deps[0] + deps[1] + deps[2]) * dt;
	et = e + de;
	if (et < 0.0)
		et = 0.0;
	porosity = et / (1 + et);
	e = et;
}

void Particle::stress_adjust() {
	if (vsig > 0.0) {
		double temp = -vsig;
		sig[0] = sig[0] + temp;
		sig[1] = sig[1] + temp;
		sig[2] = sig[2] + temp;
		sig[3] = 0.0;
		sig[4] = 0.0;
		sig[5] = 0.0;
	}
}

//initializing of StiffMat
StiffMat::StiffMat()
{
	int i, j;

	for (i = 0; i < 6; i++)
	{
		for (j = 0; j < 6; j++)
		{
			depp0[i][j] = 0.0;
		}
	}
}

StiffMat::~StiffMat()
{
}

clPair_SS::clPair_SS()
{
	p1id = -1;
	p2id = -1;
	p3id = -1;

	dp = 0.0;
	ms = 0.0;
	vs[0] = 0.0;
	vs[1] = 0.0;
	vs[2] = 0.0;

	vn[0] = 0.0;
	vn[1] = 0.0;
	vn[2] = 0.0;

	vt[0] = 0.0;
	vt[1] = 0.0;
	vt[2] = 0.0;

	ms_2p = 0.0;
	vs_2p[0] = 0.0;
	vs_2p[1] = 0.0;
	vs_2p[2] = 0.0;
}

clPair_SS::~clPair_SS()
{
}

void clPair_SS::info_reset()
{
	p1id = -1;
	p2id = -1;
	p3id = -1;

	dp = 0.0;
	ms = 0.0;
	vs[0] = 0.0;
	vs[1] = 0.0;
	vs[2] = 0.0;

	vn[0] = 0.0;
	vn[1] = 0.0;
	vn[2] = 0.0;

	vt[0] = 0.0;
	vt[1] = 0.0;
	vt[2] = 0.0;

	ms_2p = 0.0;
	vs_2p[0] = 0.0;
	vs_2p[1] = 0.0;
	vs_2p[2] = 0.0;
}

clInterf_SS::clInterf_SS()
{
	f_n[0] = 0.0;
	f_n[1] = 0.0;
	f_n[2] = 0.0;

	f_t[0] = 0.0;
	f_t[1] = 0.0;
	f_t[2] = 0.0;
}

clInterf_SS::~clInterf_SS()
{
}

void clInterf_SS::info_reset()
{
	f_n[0] = 0.0;
	f_n[1] = 0.0;
	f_n[2] = 0.0;

	f_t[0] = 0.0;
	f_t[1] = 0.0;
	f_t[2] = 0.0;
}

clParti_Pair::clParti_Pair()
{
	int i;
	total = 0;
	for (i = 0; i < 100; i++)
	{
		parti_id[i] = -1;
		parti_dst[i] = 500000.0;
	}
}

clParti_Pair::~clParti_Pair()
{
}

void clParti_Pair::info_reset()
{
	int i;
	total = 0;
	for (i = 0; i < 100; i++)
	{
		parti_id[i] = -1;
		parti_dst[i] = 5000000.0;
	}
}

void clParti_Pair::dist_sorting()
{
	int i, j, tempid;

	for (i = 0; i < 4; i++)
	{
		for (j = i + 1; j < total; j++)
		{
			if (parti_dst[i] > parti_dst[j])
			{
				tempid = parti_id[i];
				parti_id[i] = parti_id[j];
				parti_id[j] = tempid;
			}
		}
	}
}

cl_bndy_pair::cl_bndy_pair()
{
	vel_bndy[0] = 0.0;
	vel_bndy[1] = 0.0;
	vel_bndy[2] = 0.0;
}


cl_bndy_pair::~cl_bndy_pair()
{
}

clRaij::clRaij() {
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++)
			Rij[i][j] = 0.0;
	}
}

clRaij::~clRaij() {

}

clVar_Boundary::clVar_Boundary() {
	for (int i = 0; i < 3; i++) {
		vx_water[i] = 0.0;
		vx_soil[i] = 0.0;
		vx_air[i] = 0.0;
		vx_struct[i] = 0.0;

		sig_water[i] = 0.0;
		sig_soil[i] = 0.0;
		sig_air[i] = 0.0;
		sig_struct[i] = 0.0;

		sig_water[i + 3] = 0.0;
		sig_soil[i + 3] = 0.0;
		sig_air[i + 3] = 0.0;
		sig_struct[i + 3] = 0.0;
	}
	rho_water = 0.0;
	rho_soil = 0.0;
	rho_air = 0.0;
	rho_struct = 0.0;
}

clVar_Boundary::~clVar_Boundary() {

}