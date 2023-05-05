/************************************************************************************
Multi-Phase and Parallelized SPH program
--parallelized and optimized for SMP (Symmetric multiprocessing ) architecture.
--Copyright (c) 2018-2023, Weijie ZHANG, GeoHohai, Hohai University.
************************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "Class_Problem.h"

Para_FL::Para_FL()
{
	fla = 0;
	fls = 0;
	flr = 0;
	fli = 0;
	flv = 0;
	fleob = 0;
	flts = 0;
	flspv = 0;
	flknl = 1;
	flbndy = 2;

	ntd = 1;
}

Para_FL::~Para_FL()
{
}

Para_Pro::Para_Pro()
{
	ndim = 1;

	fop = 0;
	inip = 0;
	loop = 0;
	ntotal = 0;
	step_regu = 1;

	nboun = 0;
	nwater = 0;
	nsoil = 0;
	nair = 0;
	nstruct = 0;

	xsph[0] = 0.0;
	xsph[1] = 0.0;
	xsph[2] = 0.0;

	art_vis[0][0] = 0.0;
	art_vis[0][1] = 0.0;
	art_vis[1][0] = 0.0;
	art_vis[1][1] = 0.0;
	art_vis[2][0] = 0.0;
	art_vis[2][1] = 0.0;

	art_stre_coe = 0.0;

	dttype = 0;
	dr = 0.0;
	dt = 0.0;
	tt = 0.0;

	l = 1;
	t = 0.0;
}

Para_Pro::~Para_Pro()
{
}

Para_Rain::Para_Rain()
{
	rtype = 0;
	nrain = 0;
	nhor = 0;
	press_factor = 0.0;
	cement_flag = 1;
	drop_velocity = 0.0;
	x_max = -10000.0;
	x_min = 10000.0;
}

Para_Rain::~Para_Rain()
{
}

Para_Inter::Para_Inter()
{
	interc = 0.0;

	permtype = 0;
	pers = 0.0;
	perl = 0.0;
	perm = 0.0;
}

Para_Inter::~Para_Inter()
{
}

Para_Boundary::Para_Boundary()
{
	if_cal = 0;
	if_move = 0;

	move_vx[0] = 0.0;
	move_vx[1] = 0.0;
	move_vx[2] = 0.0;

	coe_bndy = 0.10;
}

Para_Boundary::~Para_Boundary()
{
}

//allocate memory for loads
Para_Vibra::Para_Vibra()
{
	loads = NULL;
	ttime = 0.0;
	nsteps = 0;

	sttime = 0.0;
	edtime = 0.0;
}

void Para_Vibra::memloads()
{
	loads = new double[nsteps][3];
}
//delete memory for loads
Para_Vibra::~Para_Vibra()
{
	if (loads != NULL)
		delete[] loads;
	loads = NULL;
}

Para_GF::Para_GF()
{
	gx = 0.0;
	gy = 0.0;
	gz = 0.0;
}

Para_GF::~Para_GF()
{
}

Para_EC::Para_EC()
{
	ecp = NULL;
	ecount = 0;
}

Para_EC::~Para_EC()
{
	if (ecp != NULL)
		delete[] ecp;
	ecp = NULL;
}
void Para_EC::memloads()
{
	ecp = new int[ecount][3];
}

Para_SSInter::Para_SSInter()
{
	type = 0;
	myu = 0.0;
	zeta = 0.5;
	acc_type = 1;
}

Para_SSInter::~Para_SSInter()
{
}

clPoint::clPoint()
{
	x = 0.0;
	y = 0.0;
	z = 0.0;
}

clPoint::~clPoint()
{
}

clLine::clLine()
{
	A = 0.0;
	B = 0.0;
	C = 1.0;

	n_x = 0.0;
	n_y = 0.0;
}

void clLine::set_normal()
{
	double total;
	total = sqrt(A * A + B * B);

	if (total > 1.0e-10)
	{
		n_x = A / total;
		n_y = B / total;
	}
}

clLine::~clLine()
{
}

clLine3D::clLine3D()
{
	x0 = 0.0;
	y0 = 0.0;
	z0 = 0.0;

	a = 0.0;
	b = 0.0;
	c = 0.0;
}

clLine3D::~clLine3D()
{
}

clPlane::clPlane()
{
	A = 0.0;
	B = 0.0;
	C = 0.0;
	D = 1.0;

	n_x = 0.0;
	n_y = 0.0;
	n_z = 0.0;
}

void clPlane::set_normal()
{
	double total;
	total = sqrt(A * A + B * B + C * C);

	if (total > 1.0e-10)
	{
		n_x = A / total;
		n_y = B / total;
		n_z = C / total;
	}
}

clPlane::~clPlane()
{
}