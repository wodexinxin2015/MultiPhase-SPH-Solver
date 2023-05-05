/************************************************************************************
Multi-Phase and Parallelized SPH program
--parallelized and optimized for SMP (Symmetric multiprocessing ) architecture.
--Copyright (c) 2018-2023, Weijie ZHANG, GeoHohai, Hohai University.
************************************************************************************/

#include <cmath>
#include <algorithm>
#include "Class_Functions.h"
#include "Header_Option.h"

//excavation or backfill
void clExcBac_Fun::exca_or_back(Particle *pPar, Para_Pro *pPPro, const Para_EC &pPEc, int cn, int nstage)
{

	int i;
	int type, start, end;

	start = pPEc.ecp[nstage][0];
	end = pPEc.ecp[nstage][1];

#pragma omp parallel for schedule(static) private(type)
	for (i = start - 1; i < end; i++)
	{
		type = pPar[i].type;
		pPar[i].type = pPar[i].etype;
		pPar[i].etype = type;

		if (pPar[i].type == 1)
		{
#pragma omp critical
			{
				pPPro->nwater += 1;
			}
		}
		else if (pPar[i].type == 2)
		{
#pragma omp critical
			{
				pPPro->nsoil += 1;
			}
		}
		else if (pPar[i].type == 4)
		{
#pragma omp critical
			{
				pPPro->nstruct += 1;
			}
		}
		else if (pPar[i].type == 7)
		{
			if (type == 1)
			{
#pragma omp critical
				{
					pPPro->nwater -= 1;
				}
			}
			else if (type == 2)
			{
#pragma omp critical
				{
					pPPro->nsoil -= 1;
				}
			}
			else if (type == 4)
			{
#pragma omp critical
				{
					pPPro->nstruct -= 1;
				}
			}
		}
	}
}