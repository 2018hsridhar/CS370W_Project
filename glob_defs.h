#ifndef GLOBS_DEF_H
#define GLOBS_DEF_H

//#include "path.h"
#include "path.h"

// #TODO :: read up how globals.h is set up as 
// #TODO :: describe what variables actually dop
// #TODO Put all the globally-used numbers in here.
// #TODO :: capitalize all consts [ generally a good idea ]

// Avoid the use of following
// [1] magic numbers 
// [2] multiple defines.

namespace GLOBAL
{
	// most of these variables - for debugging of pipeline steps mostly. Not much for other tasks.
	extern const char* stitching_output;
	extern const char* remeshInputFile;
	extern const char* sleOutputFile;
	extern const char* isoOutputFile;
    extern const char* remeshOutputFile;
	extern const char* remesh_after_succ_file;
	extern const char* mcfDebugText;
	extern const char* surfaceBlowUpCase;

	extern const char* stitchedSurfGenScan1File;
	extern const char* stitchedSurfGenScan2File;

	extern const char* rigidIcpInputScanFile;
	extern const char* rigidIcpOutputScanFile;

	extern const char* meanCurvFlowInputScanFile;
	extern const char* meanCurvFlowOutputScanFile;

	extern const char* viewMesh;

	// INSPIRED from Randall's code for rendering
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

	extern const unsigned int MCF_STEPS; // for MCF ... maybe ... later
	extern const unsigned int MAXITERS;	 // for RigidICP
	extern const double TOLERANCE;  // for RigidICP

	extern const char* pipelineScan1File;
	extern const char* pipelineScan2File;
	extern const char* pipelineOutputFile;

	extern const char* pipelineDebugFile;
	extern const char* mcfDebugFile;

	extern const char* boolOpersMeshOne;
	extern const char* boolOpersMeshTwo;
	extern const char* boolOpersFragOne;
	extern const char* boolOpersFragTwo;

	// triangle-test
	extern const char* simpleTriangleFile;
	extern const char* twoTrianglesFile;
}
#endif




