/************************************************************************************
Multi-Phase and Parallelized SPH program
--parallelized and optimized for SMP (Symmetric multiprocessing ) architecture.
--Copyright (c) 2018-2023, Weijie ZHANG, GeoHohai, Hohai University.
************************************************************************************/

#ifndef CLASS_PARTICLE_H_INCLUDE
#define CLASS_PARTICLE_H_INCLUDE

#include "Class_Problem.h"

//particles
class Particle
{
public:
	//variables
	double x[3];   /*position vector */
	double xp[3];  /*previous position vector */
	double vx[3];  /*velocity vector or boundary virtual velocity for water*/
	double vxp[3]; /*previous velocity vector or boundary virtual velocity for soil*/
	double ax[3];  /*acceleration vector or boundary virtual velocity for air*/
	double axp[3]; /*previous acceleration vector or boundary virtual velocity for structure*/

	double ux[3]; /*displacement*/

	double strep[3]; /*eig value for stress also principle stress*/
					 /*2016-7-7*/

	double divx[3][3];	/*velocity derivatives*/
	double interf[3];	/*interaction forces 
					  or boundary wall velocity for seismic or moving acceleration*/
	double interfss[3]; /*interaction forces of soil-structure considering the frictional force 
						or boundary wall velocity for seismic or moving velocity*/

	double rho;	 /*density */
	double rhop; /*previous  density*/
	double mass; /*mass */
	double hl;	 /*smooth length */

	double pre;	 /*dynamic pressure*/
	double prep; /*previous pressure*/

	double e; /*void ratio*/

	double eps[6];	   /*strain */
	double deps[6];	   /*strain rate*/
	double sig[6];	   /*stress*/
	double dsig[6];	   /*stress rate*/
	double ff_sig[3];  /*stress component for free field particles*/
	double weps[3][3]; /*derivate strain rate */
	double adps[6];	   /*accumulated plastic strain*/
	/*V2.7 2016-7-25*/

	double tssig[6]; /*stress after TS transformation*/
	double tsqc;	 /*deviatoric stress q after TS transformation*/

	double satu;  /*degree of saturation */
	double dsatu; /*degree of saturation */
	double suct;  /*suction */
	double cop;	  /*permeability coefficient */

	double vsig; /*volumetric stress*/
	double veps; /*volumetric strain*/
	double divq; /*deviatoric stress q */
	double divr; /*deviatoric stress ratio ita*/
	double meps; /*maximium of shear strain*/

	double roue;	 /*state variable for overconsolidation*/
	double rous;	 /*state variable for degree of saturation*/
	double zeta;	 /*state variable for anisotropy*/
	double rsta;	 /*state variable for structure*/
	double prew;	 /*water pressure for soil*/
	double prea;	 /*air pressure for soil*/
	double porosity; /*porosity*/
	double cd;		 /*damping coefficient*/
	double poro_cop; /*void ratio change by the soil deformation*/

	double beta[3][3];	  /*transformation matrix of stress tensor*/
	double beta_ts[3][3]; /*transposed transformation matrix of stress tensor*/

	int type; /*type of particle */
	//0, 1, 2, 3, 4, and 7
	int matype; /*sub type of material */
	/*0, 10, 20; 1, 11, 21; 2, 12, 22; 3, 13, 23, 33;4, 14, 24;*/
	/*5, 15, 25; 6, 16, 26; 7, 17, 27; 8, 18, 28; 9, 19, 29*/
	//>10 indicates different layers
	int permtype; /*permeability type*/

	int etype; /*existing type: 0--not exist; 1--exist*/

	//methods
	void massini(int dim);												   /*calculating mass for particles*/
	void deq();															   /*the function of derivation shear stress*/
	void permeablel(double sitas, double sitar, double pers, double perl); /*two layers*/

	void permeabletl(double pers, double perl, double perm); /*three layers*/
	void shearstrain();										 /*calculate the sqrt(2*I2) : I2 is the second invariant of deviatoric strain tensor*/
	void poro_cal(double);									 /*calculate the porosity*/
	void stress_adjust();

	Particle();
	~Particle();
};

class StiffMat
{
public:
	//variables
	double depp0[6][6]; /*stiffness matrix at t0*/

	//methods
	StiffMat();
	~StiffMat();
};

//defining the pair information of soil-structure interaction
class clPair_SS
{
public:
	int p1id; //id of neighboring structure particle 1 for 2D and 3D
	int p2id; //id of neighboring structure particle 2 for 2D and 3D
	int p3id; //id of neighboring structure particle 3 for 3D

	clPoint p1; //coordinate of neighboring structure particle 1 for 2D and 3D
	clPoint p2; //coordinate of neighboring structure particle 2 for 2D and 3D
	clPoint p3; //coordinate of neighboring structure particle 3 for 3D

	double dp;	  //distance between moving particle and structure particle for 2D and 3D
	double ms;	  //mass of mappping point for 2D and 3D
	double vs[3]; //velocity of mapping point for 2D and 3D

	double vn[3]; //normal velocity vector
	double vt[3]; //tangential velocity vector

	double ms_2p;	 //mass of mappping point for 3D
	double vs_2p[3]; //velocity of mapping point for 3D

	//methods
	clPair_SS();
	~clPair_SS();
	void info_reset();
};

class clInterf_SS
{
public:
	double f_n[3]; //normal interaction force between soil particle and structure particle
	double f_t[3]; //tangential interaction force between soil particle and structure particle
	//methods
	clInterf_SS();
	~clInterf_SS();
	void info_reset();
};

class clParti_Pair
{
	//soil-structure interaction pair information
public:
	int total;			   //total contacting number of particles
	int parti_id[100];	   //contacting structure particle id
	double parti_dst[100]; //contacting structure partcle distance

	//method
	clParti_Pair();
	~clParti_Pair();
	void info_reset();
	//distance sorting
	void dist_sorting();
};

class cl_bndy_pair
{
public:
	double vel_bndy[3];

	//method
	cl_bndy_pair();
	~cl_bndy_pair();
};

class clRaij {
public:
	double Rij[3][3];

	clRaij();
	~clRaij();
};


class clVar_Boundary {
public:
	double vx_water[3]; //boundary particle velocity for water phase;
	double vx_soil[3]; //boundary particle velocity for soil phase;
	double vx_air[3]; //boundary particle velocity for air phase;
	double vx_struct[3]; //boundary particle velocity for structure phase;
	double sig_water[6]; //boundary particle stress for water phase;
	double sig_soil[6]; //boundary particle stress for soil phase;
	double sig_air[6]; //boundary particle stress for air phase;
	double sig_struct[6]; //boundary particle stress for structure phase;
	double rho_water; //density of water
	double rho_soil; //density of soil
	double rho_air; //density of air
	double rho_struct; //density of structure

	clVar_Boundary();
	~clVar_Boundary();
};
#endif // CLASS_PARTICLE_H_INCLUDE