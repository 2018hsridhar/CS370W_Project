#ifndef GLOBS_DEF_H
#define GLOBS_DEF_H

#include "tutorial_shared_path.h"

// #TODO :: read up how globals.h is set up as 
// note :: most of these should be a bunch of consts
// and you should also describe what they are

// Avoid the use of magic numbers and multiple defines.
// Put all the globally-used numbers in here.

// remember :: constants are ALWAYS capitalized!

// #TODO :: in glob_defs, have variables you can set ( for say, edge lengths in CGAl code, or files you want to test ... wait, for files you want to tet, that's a good @(#)$u idea ) . Your files should or should not be different, based on what you do and do not want to test. A.K.A. make it configurable!

namespace GLOBAL
{
	extern const char* remeshInputFile;
    extern const char* remeshOutputFile;

	extern const char* offsetGenScan1File;
	extern const char* offsetGenScan2File;

	extern const char* rigidIcpInputScanFile;
	extern const char* rigidIcpOutputScanFile;

	extern const char* meanCurvFlowInputScanFile;
	extern const char* meanCurvFlowOutputScanFile;

	extern const char* viewMesh;

	// inspired from Randall's code for rendering
	// to adopt for later
	/*
    extern const char* output_directory;
    extern const char* model_path;
    extern unsigned int window_width;
    extern unsigned int window_height;
    extern unsigned int num_samples;
	extern enum sample_type; // unsure if useful
    extern enum algo_comp; // unsure if useful
	*/

	extern const unsigned int MAXITERS;	
	extern const unsigned int MCF_STEPS;
	extern const double TOLERANCE; 
	
}
#endif




