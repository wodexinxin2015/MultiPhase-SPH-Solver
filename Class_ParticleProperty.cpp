/************************************************************************************
Multi-Phase and Parallelized SPH program
--parallelized and optimized for SMP (Symmetric multiprocessing ) architecture.
--Copyright (c) 2018-2023, Weijie ZHANG, GeoHohai, Hohai University.
************************************************************************************/

#include "Class_ParticleProperty.h"

Para_Model::Para_Model()
{
	flagem = 0;
	flagsb = 0;
	flagdp = 0;
	flagum = 0;
	flagam = 0;
	flagbui = 0;
	flagpl = 0;
	flagve = 0;
	flagfocc = 0;
	flagep = 0;
	alphats = 0.0;
}

Para_Model::~Para_Model()
{
}

Para_Soil::Para_Soil()
{
	e = 0.0;
	poi = 0.0;
	eps0 = 0.0;

	dens = 0.0;

	damp1 = 0.0;
	damp2 = 0.0;

	cop = 0.0;

	ramda = 0.0;
	kapa = 0.0;
	zmf = 0.0;

	ann = 0.0;
	bnn = 0.0;
	beta = 0.0;

	sitas = 0.0;
	sitar = 0.0;
	epsr = 0.0;
	epss = 0.0;

	sd = 0.0;
	sw = 0.0;
	kse = 0.0;
	c1 = 0.0;
	c2 = 0.0;
	c3 = 0.0;

	fai = 0.0;
	c = 0.0;

	ocr = 0.0;
	sme0 = 0.0;

	mzm = 0.0;
	azm = 0.0;
	r0 = 0.0;
	bl = 0.0;
	bz = 0.0;
	zet = 0.0;

	vpu = 0.0;
	porosity = 0.0;
	fs = 1.0;

	etao = 0.0;
	g0 = 0.0;

	ds = 1.0;

	damp_coe = 0.0;
}

Para_Soil::~Para_Soil()
{
}

Para_Fluid::Para_Fluid()
{
	vis = 0.0;
	vpu = 0.0;
	dens = 0.0;
	paa0 = 0.0;

	porosity = 0.0;

	fai = 0.0;
	c = 0.0;

	bhtype = 0;
	press_type = 0;
	acc_type = 0;
	pre_coe = 0.0;

	strain_min = 1.0;
}

Para_Fluid::~Para_Fluid()
{
}

Para_Structure::Para_Structure()
{
	dens = 0.0;
	eps0 = 0.0;
	e = 0.0;
	poi = 0.0;

	vpu = 0.0;
	porosity = 0.0;
	cop = 0.0;
}

Para_Structure::~Para_Structure()
{
}