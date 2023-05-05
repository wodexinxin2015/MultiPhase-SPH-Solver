/************************************************************************************
Multi-Phase and Parallelized SPH program
--parallelized and optimized for SMP (Symmetric multiprocessing ) architecture.
--Copyright (c) 2018-2023, Weijie ZHANG, GeoHohai, Hohai University.
************************************************************************************/

#ifndef CLASS_PARTICLEPROPERTY_H_INCLUDE
#define CLASS_PARTICLEPROPERTY_H_INCLUDE

class Para_Model
{
public:
	//flag** must be large than 1 and it also means the number of layers
	//flag** must be less than or equation to 3
	int flagem;		/*elastic*/
	int flagsb;		/*modified cam-clay*/
	int flagdp;		/*drucker-prager*/
	int flagum;		/*unsaturated cam-clay*/
	int flagam;		/*anisotropic cam-clay*/
	int flagbui;	/*drucker prager model proposed by Bui H.H.*/
	int flagpl;		/*pile material*/
	int flagve;		/*flag for viscous elastic*/
	int flagfocc;	/*flag for fractional order plastic model of MCC*/
	int flagep;		/*flag for the M-C failure criteria*/
	double alphats; /*alpha coefficient for the stress transformation of
					YP Yao 2015, Chinese Journal of Geotechnical Engineering*/

	//methods
	Para_Model();
	~Para_Model();
};

class Para_Soil
{
public:
	double e;	 /*elastic module*/
	double poi;	 /*Poisson's ratio*/
				 /*paramter poison for FOPM*/
	double eps0; /*Initial void ratio*/

	double dens; /*density*/

	double damp1; /*damping coeffificent h1--alpha damping in Rayleigh damping*/
				  /*parameter h1 for FOPM*/
	double damp2; /*damping coeffificent h2--beta damping in Rayleigh damping*/
				  /*parameter h2 for FOPM*/

	double cop; /*coeffificent of permeability*/

	double ramda; /*cambridge model: lamda*/
				  /*parameter ramda for FOPM*/
	double kapa;  /*cambridge model: kapa*/
	double zmf;	  /*critical stress ratio*/
				  /*paramter mc for FOPM*/

	double ann;	 /*consititutive model controlling paramters*/
				 /*parameter alpha for FOPM*/
	double bnn;	 /*consititutive model controlling paramters*/
				 /*parameter beta for FOPM*/
	double beta; /*consititutive model controlling paramters*/

	double sitas; /*Saturation: maximum*/
	double sitar; /*Saturation: minimum*/
	double epsr;  /*void ratio for dry state*/
	double epss;  /*void ratio for saturated state*/

	double sd;	/*moisture characteristic curve: sd*/
	double sw;	/*moisture characteristic curve: sw*/
	double kse; /*moisture characteristic curve: kse*/
	double c1;	/*moisture characteristic curve: c1*/
	double c2;	/*moisture characteristic curve: c2*/
	double c3;	/*moisture characteristic curve: c3*/

	double fai; /*internal frictional angle: degree*/
	double c;	/*cohesion*/
	/*the ratio between Me and Mc*/

	double ocr;	 /*overconsolidation ratio*/
	double sme0; /*previous consolidation pressure*/

	double mzm; /*parameter for overconsolidation*/
	double azm; /*parameter for structure*/
	double r0;	/*initial structure ratio*/
	double bl;	/*parameter for anisotropy, default=0.95*/
	double bz;	/*parameter for anisotropy*/
	double zet; /*initial anisotropy*/
				/*parameter zeta for FOPM*/

	double vpu;		 /*sound velocity*/
	double porosity; /*porosity*/
	double fs;		 /*factor of safety*/

	//the parameters of Fractional Order Plastic Model
	double etao; /*parameter etao for FOPM*/
	double g0;	 /*parameter G0 for FOPM*/

	//parameter for the soil particle size
	double ds; /*average diameter of soil particles*/

	double damp_coe; /*damping coefficient for the damping effect in Bui's research*/

	//methods
	Para_Soil();
	~Para_Soil();
};

class Para_Fluid
{
public:
	double vis; /*viscosity*/
	double vpu; /*sound velocity*/

	double dens; /*density*/
	double paa0; /*initial pressure*/

	double porosity; /*porosity*/

	double fai; /*frictional angle in Mohr��Coulomb model*/
	double c;	/*internal cohesion in Mohr��Coulomb model*/

	int bhtype; /*Bingham flowing type*/

	int press_type; /*type of pressure used for the calculation of acceleration*/
	/*0--incremental pressure; 1-- total pressure*/
	int acc_type; /*type of the calculation of acceleration*/
	/*0--plus stress; 1--subtract stress*/
	double pre_coe; /*coefficient of pressure for initializing*/
	double strain_min; /*minimum shear strain for the Bingham fluid calculation*/

	//methods
	Para_Fluid();
	~Para_Fluid();
};

class Para_Structure
{
	//constitutive model parameters for structure
public:
	double dens; //density
	double eps0; //initial void ratio
	double e;	 //elastic module
	double poi;	 //poison ratio

	double vpu;		 /*sound velocity*/
	double porosity; /*porosity*/
	double cop;		 /*coeffificent of permeability*/

	//methods
	Para_Structure();
	~Para_Structure();
};

#endif // CLASS_PARTICLEPROPERTY_H_INCLUDE
