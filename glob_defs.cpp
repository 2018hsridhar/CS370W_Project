#ifndef GLOBS_DEF_H
#define GLOBS_DEF_H

#include "path.h"
#include "glob_defs.h"

namespace GLOBAL
{

	// #TODO :: usea  different directory aside from <shared>. 
	// These two directories should've been seperated in the beginning.
	// this is getting incredibly frustrating.
	// NOTE :: only for the two inputs should the shared paths be used. Not for 
	const char* stitchedSurfGenScan1File = SHARED_PATH "/planexy.off";
    const char* stitchedSurfGenScan2File = SHARED_PATH "/mirrorPlane.off";

	const char* stitching_output = PROC_PATH "/stitching.off";
	const char* remeshInputFile = PROC_PATH "/offset.off";
	const char* sleOutputFile = PROC_PATH "/sle_post.off"; 
	const char* isoOutputFile = PROC_PATH "/iso_post.off";
    const char* remeshOutputFile = PROC_PATH "/offsetRemesh.off";
	const char* remesh_after_succ_file = PROC_PATH "/remesh_after_succ.off";
	const char* mcfDebugText = PROC_PATH "/mcf_debug.txt";
	const char* surfaceBlowUpCase = PROC_PATH "/surfaceBlowUpCase.off";

	//const char* stitchedSurfGenScan1File = SHARED_PATH "/planexy.off";
    //const char* stitchedSurfGenScan2File = SHARED_PATH "/mirrorPlane.off";
//	const char* remeshInputFile = SHARED_PATH "/offset.off";
 //   const char* remeshOutputFile = SHARED_PATH "/offsetRemesh.off";

	const char* rigidIcpInputScanFile = SHARED_PATH "/bunny.off";
	const char* rigidIcpOutputScanFile = SHARED_PATH "/bunny90.off";

	const char* meanCurvFlowOutputScanFile = SHARED_PATH "/mcf.off";

	// a nice user setting, for my configuration!
	const char* viewMesh = SHARED_PATH "/badstitchedSurf.off";

//	const char* pipelineScan1File = SHARED_PATH "/2triangles.off";
	// #TODO :: test non-animated pipeline here with boundary mesh
	//const char* pipelineScan1File = SHARED_PATH "/planexy.off";
	//const char* pipelineScan1File = SHARED_PATH "/testBoundary.off";
	//const char* pipelineScan2File = SHARED_PATH "/mirrorPlane.off";
	const char* pipelineScan1File = SHARED_PATH "/camelhead.off";
	const char* pipelineScan2File = SHARED_PATH "/camelhead2.off";
	const char* pipelineOutputFile = SHARED_PATH "/result.off";

	// DEBUG VARS
	const char* pipelineDebugFile = SHARED_PATH "/debug.txt";
	const char* mcfDebugFile = SHARED_PATH "/mcfLog.txt";

	// BOOL OPERS FILES
	// NOTE :: <.OBJ> files are used since <blender> lacks <.OFF> file support.
	const char* boolOpersMeshOne = SHARED_PATH "/sphere.obj";
	const char* boolOpersMeshTwo= SHARED_PATH "/sphere.obj";
	//const char* boolOpersMeshTwo= SHARED_PATH "/armadillo.obj";
	//const char* boolOpersMeshOne = SHARED_PATH "/armadillo.obj";
	//const char* boolOpersMeshOne = SHARED_PATH "/cube.obj";
	//const char* boolOpersMeshTwo = SHARED_PATH "/cube.obj";
	const char* boolOpersFragOne = "fragment1.obj";
	const char* boolOpersFragTwo = "fragment2.obj";
}
#endif




