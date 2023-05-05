/************************************************************************************
Multi-Phase and Parallelized SPH program
--parallelized and optimized for SMP (Symmetric multiprocessing ) architecture.
--Copyright (c) 2018-2023, Weijie ZHANG, GeoHohai, Hohai University.
************************************************************************************/

#ifndef CLASS_PROBLEM_H_INCLUDE
#define CLASS_PROBLEM_H_INCLUDE

class Para_FL
{
public:
	/*flags controling air phase, soil phase (number of materials), rianfall,
	nteraction, vibration, excavation or backfill, regularization, 
	stress correction, transformed stress for 3D Space, spatial variability,
	flag for kernel function*/
	int fla, fls, flr, fli, flv, fleob, flts, flspv, flknl;
	//fli: 1--interaction force of Darcy's law;
	//2--interaction including laminar flow and turbulent flow by Shi Zm et al. 2018
	//3--interaction force simulation by Peng C. et al. 2017
	//flknl: 1--B splines function; 2--kernel function by Yang and Liu 2016
	//3--B splines function with normalised kernel gradient; 
	//4--kernel function by Yang and Liu 2016 with normalized kernel gradient
	//flspv: flag for spatial variability generation
	//-1--the random field input from Parameters.txt
	//1--using the K-L expansion method for the non-stationary random field; 
	//2--using the K-L expansion method for the non-stationary random field with correlated c and fai; 
	//3--using the correlation matrix decomposing method for the non-stationary random field; 
	//4--using the correlation matrix decomposing method for the non-stationary random field with correlated c and fai; 
	//5--using the correlated random variable for the homogeneous field
	int ntd; /*number of threads used*/
	int flbndy; //boundary selection: 1--velocity correction method; 2--virtual velocity method
	//3--no-slip boundary by Tran et al. 2019 Computers and Geotechnics without velocity setting of free particles; 
	//4--free-slip boundary by Tran et al. 2019 Computers and Geotechnics without velocity setting of free particles; 
	//5--using the free particle-solid wall interaction for the boundary treatment;
	//13--no-slip boundary by Tran et al. 2019 Computers and Geotechnics with velocity setting of free particles; 
	//14--free-slip boundary by Tran et al. 2019 Computers and Geotechnics with velocity setting of free particles.

	//methods
	Para_FL();
	~Para_FL();
};

class Para_Pro
{
	//Problem parameters
public:
	int ndim;	/*dimensions*/

	int fop;	/*steps for output to file*/
	int inip;	/*steps for initial stress of soil*/
	int loop;	/*total loops*/
	int ntotal; /*total number of particles*/
	int step_regu; /*steps for strain and stress regularization*/

	int nboun;	 /*number of boundary particles*/
	int nwater;	 /*number of water particles*/
	int nsoil;	 /*number of soil particles*/
	int nair;	 /*number of air particles*/
	int nstruct; /*number of structure particles*/

	int dttype; /*time step type 1-constant 2-variable*/
	double dr;	/*particle spacing*/
	double dt;	/*time increment*/
	double tt;  /*total simulation time*/

	double xsph[3]; /*xsph coefficient for the subroutine of xsph*/
	//0--water; 1--soil; 2--air
	double art_vis[3][2]; /*artificial viscosity coefficient*/
	//0--water; 1--soil; 2--air
	double art_stre_coe; /*coefficient in the calculation of artificial stress*/
	/*if art_stre_coe < 0, no artifical stress; 
	if art_stre_coe > 0, using artificial stress and art_stre_coe is the coefficient.*/

	int l;	  /*current step*/
	double t; /*time*/

	//methods
	Para_Pro();
	~Para_Pro();
};

class Para_Rain
{
	//rainfall parameters-
public:
	int rtype;			  /*rainfall type*/
	int nrain;			  /*total number of rainfall particles*/
	int nhor;			  /*rain particles at every line*/
	double press_factor;  /*factor of initial pressure for rain or cement*/
	double drop_velocity; /*drop velocity of rainfall*/

	int cement_flag; /*controlling the order of cement particles*/
	double x_max;	 /*controlling the deletion of water particles--left coordinate*/
	double x_min;	 /*controlling the deletion of water particles--left coordinate*/

	//methods
	Para_Rain();
	~Para_Rain();
};

class Para_Inter
{
	//interaction force parameters
public:
	double interc; /*interaction force controlling*/

	int permtype; /*permeability calculation type*/
	double pers;  /*small permeability*/
	double perl;  /*large permeability*/
	double perm;  /*middle permeability*/

	//methods
	Para_Inter();
	~Para_Inter();
};

class Para_Boundary
{
	//boudary parameters
public:
	double coe_bndy; //boundary coefficient for the velocity_set_domain of type 2, 3 and 4
	int if_cal;	 //controlling the calculation of impact forces on boundary
	int if_move; //if there are moveable boundary particles: 0--no; 1--yes

	double move_vx[3]; //velocity for the moveable particles

	//methods
	Para_Boundary();
	~Para_Boundary();
};

class Para_Vibra
{
	//Vibration loads
public:
	double ttime; /*one period time*/
	int nsteps;	  /*steps of one period time*/

	double sttime; /*starting time*/
	double edtime; /*end time*/

	double (*loads)[3]; /*loading value at every steps: m/s/s*/

	//methods
	void memloads();
	Para_Vibra();
	~Para_Vibra();
};

class Para_GF
{
	//gravity field
public:
	double gx; //gravity acceleration at x
	double gy; //gravity acceleration at y
	double gz; //gravity acceleration at z

	//methods
	Para_GF();
	~Para_GF();
};

class Para_EC
{
	//excavation or backfill parameters
public:
	int ecount;	   //total number of excavation or backfill stages
	int (*ecp)[3]; //controlling parameters: 0-estart; 1-eend; 2-ecpnp;

	//methods
	void memloads();
	Para_EC();
	~Para_EC();
};

class Para_SSInter
{
	//parameters for soil and structure interaction
public:
	int type; //type of soil and structure interaction:
	//1--rigid structure; 2--deformable structure; 3--penalty method
	double myu;	 //frictional coefficient
	double zeta; //coefficient of soil-structure interaction force
	//zeta is for the rigid method, deformable structure method and penalty method
	int acc_type; //soil or water and structure accelaration calculation type
	//1--seperate calculation of accelaration; 2--coupled calculation of accelaration

	//methods
	Para_SSInter();
	~Para_SSInter();
};

//defining a point class
class clPoint
{
public:
	double x; //x coordinate
	double y; //y coordinate
	double z; //z coordinate

	//methods
	clPoint();
	~clPoint();
};

//defining a line class
class clLine
{
public:
	double A; //coefficient A
	double B; //coefficient B
	double C; //coefficient C

	double n_x; //unit normal vector x
	double n_y; //unit normal vector y

	//methods
	void set_normal(void);
	clLine();
	~clLine();
};

//defining a line in 3D space
class clLine3D
{
public:
	double x0; //start point x
	double y0; //start point y
	double z0; //start point z

	double a; //direction vector x
	double b; //direction vector y
	double c; //direction vector z

	//methods
	clLine3D();
	~clLine3D();
};

//defining a plane class
class clPlane
{
public:
	double A; //coefficient A
	double B; //coefficient B
	double C; //coefficient C
	double D; //coefficient D

	double n_x; //unit normal vector x
	double n_y; //unit normal vector y
	double n_z; //unit normal vector z

	//methods
	void set_normal(void);
	clPlane();
	~clPlane();
};
#endif // CLASS_PROBLEM_H_INCLUDE