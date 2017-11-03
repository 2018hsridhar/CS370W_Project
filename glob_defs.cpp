#ifndef GLOBS_DEF_H
#define GLOBS_DEF_H

#include "tutorial_shared_path.h"
#include "glob_defs.h"

namespace GLOBAL
{

	// #TODO :: usea  different directory aside from <shared>. 
	// These two directories should've been seperated in the beginning.
	// this is getting incredibly frustrating.
	const char* interpSurfGenScan1File = TUTORIAL_SHARED_PATH "/planexy.off";
    const char* interpSurfGenScan2File = TUTORIAL_SHARED_PATH "/mirrorPlane.off";
	const char* remeshInputFile = TUTORIAL_SHARED_PATH "/offset.off";
	const char* sleOutputFile = "sle_post.off";
	const char* isoOutputFile = "iso_post.off";
    const char* remeshOutputFile = "offsetRemesh.off";

	//const char* interpSurfGenScan1File = TUTORIAL_SHARED_PATH "/planexy.off";
    //const char* interpSurfGenScan2File = TUTORIAL_SHARED_PATH "/mirrorPlane.off";
//	const char* remeshInputFile = TUTORIAL_SHARED_PATH "/offset.off";
 //   const char* remeshOutputFile = TUTORIAL_SHARED_PATH "/offsetRemesh.off";

	const char* rigidIcpInputScanFile = TUTORIAL_SHARED_PATH "/bunny.off";
	const char* rigidIcpOutputScanFile = TUTORIAL_SHARED_PATH "/bunny90.off";

	const char* meanCurvFlowOutputScanFile = TUTORIAL_SHARED_PATH "/mcf.off";

	// a nice user setting, for my configuration!
	const char* viewMesh = TUTORIAL_SHARED_PATH "/badInterpSurf.off";

//	const char* pipelineScan1File = TUTORIAL_SHARED_PATH "/2triangles.off";
	// #TODO :: test non-animated pipeline here with boundary mesh
	//const char* pipelineScan1File = TUTORIAL_SHARED_PATH "/planexy.off";
	//const char* pipelineScan1File = TUTORIAL_SHARED_PATH "/testBoundary.off";
	//const char* pipelineScan2File = TUTORIAL_SHARED_PATH "/mirrorPlane.off";
	const char* pipelineScan1File = TUTORIAL_SHARED_PATH "/camelhead.off";
	const char* pipelineScan2File = TUTORIAL_SHARED_PATH "/camelhead2.off";
	const char* pipelineOutputFile = TUTORIAL_SHARED_PATH "/result.off";

	// DEBUG VARS
	const char* pipelineDebugFile = TUTORIAL_SHARED_PATH "/debug.txt";
	const char* mcfDebugFile = TUTORIAL_SHARED_PATH "/mcfLog.txt";

	// BOOL OPERS FILES
	// NOTE :: <.OBJ> files are used since <blender> lacks <.OFF> file support.
	const char* boolOpersMeshOne = TUTORIAL_SHARED_PATH "/sphere.obj";
	const char* boolOpersMeshTwo= TUTORIAL_SHARED_PATH "/sphere.obj";
	//const char* boolOpersMeshTwo= TUTORIAL_SHARED_PATH "/armadillo.obj";
	//const char* boolOpersMeshOne = TUTORIAL_SHARED_PATH "/armadillo.obj";
	//const char* boolOpersMeshOne = TUTORIAL_SHARED_PATH "/cube.obj";
	//const char* boolOpersMeshTwo = TUTORIAL_SHARED_PATH "/cube.obj";
	const char* boolOpersFragOne = "fragment1.obj";
	const char* boolOpersFragTwo = "fragment2.obj";
}
#endif




