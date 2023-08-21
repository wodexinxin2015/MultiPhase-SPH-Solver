/************************************************************************************
Multi-Phase and Parallelized SPH program
--parallelized and optimized for SMP (Symmetric multiprocessing ) architecture.
--Copyright (c) 2018-2023, Weijie ZHANG, GeoHohai, Hohai University.
************************************************************************************/

#ifndef HEADER_OPTION_H_INCLUDED
#define HEADER_OPTION_H_INCLUDED

#ifdef __ANDROID__
#define win32 0
#elif __linux__
#define win32 0
#elif _WIN32
#define win32 1
#endif

#define ep_err 1.0e-3

#endif // HEADER_OPTION_H_INCLUDED
