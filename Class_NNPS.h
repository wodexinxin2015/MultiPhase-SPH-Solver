/************************************************************************************
Multi-Phase and Parallelized SPH program
--parallelized and optimized for SMP (Symmetric multiprocessing ) architecture.
--Copyright (c) 2018-2023, Weijie ZHANG, GeoHohai, Hohai University.
************************************************************************************/

#ifndef CLASS_NNPS_H_INCLUDE
#define CLASS_NNPS_H_INCLUDE

// define the class for NNPS nearest neighboring particle searching

class Cell_Con
{ // controlling information of cells
public:
	int cxm;	/*total number of cells at x*/
	int cym;	/*total number of cells at y*/
	int czm;	/*total number of cells at z*/
	int ctotal; /*total number of cells*/

	double xmin, xmax; /*problem domain*/
	double ymin, ymax; /*problem domain*/
	double zmin, zmax; /*problem domain*/

	double xmin_nb, xmax_nb; /*problem domain of non boundary*/
	double ymin_nb, ymax_nb; /*problem domain of non boundary*/
	double zmin_nb, zmax_nb; /*problem domain of non boundary*/

	double water_min[3], water_max[3]; /*problem domain of water phase*/
	double soil_min[3], soil_max[3];   /*problem domain of soil phase*/
	double air_min[3], air_max[3];	   /*problem domain of air phase*/
	double stru_min[3], stru_max[3];   /*problem domain of structure phase*/

	// methods
	Cell_Con();
	~Cell_Con();
};

class parti_cellid
{
public:
	int cell_id;
	int parti_id;

	parti_cellid();
	~parti_cellid();
};

class cell_link
{
public:
	int nncell[27];
	cell_link();
	~cell_link();
};

class cell_info
{
public:
	int start;
	int end;
	cell_info();
	~cell_info();
};

class Par_Cell
{ // storing NNP
public:
	int cell_id;   /*cell id for this particle*/
	int ninflu;			/*total number of supporting particles*/
	int influ[100];		/*NNP ID*/
	double wij[100][4]; /*smoothing function values*/

	// methods
	Par_Cell();
	~Par_Cell();
	void initial();
};

#endif // CLASS_NNPS_H_INCLUDE