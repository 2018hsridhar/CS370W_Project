#ifndef GLOBS_DEF_H
#define GLOBS_DEF_H

#include "tutorial_shared_path.h"
#include "glob_defs.h"

namespace GLOBAL
{
	const char* remeshInputFile = TUTORIAL_SHARED_PATH "/offset.off";
    const char* remeshOutputFile = TUTORIAL_SHARED_PATH "/offsetRemesh.off";

//	const char* interpSurfGenScan1File = TUTORIAL_SHARED_PATH "/planexy.off";
 //   const char* interpSurfGenScan2File = TUTORIAL_SHARED_PATH "/mirrorPlane.off";
	const char* interpSurfGenScan1File= TUTORIAL_SHARED_PATH "/camelhead.off";
    const char* interpSurfGenScan2File= TUTORIAL_SHARED_PATH "/camelhead2.off";

	const char* rigidIcpInputScanFile = TUTORIAL_SHARED_PATH "/bunny.off";
	const char* rigidIcpOutputScanFile = TUTORIAL_SHARED_PATH "/bunny90.off";

	const char* meanCurvFlowOutputScanFile = TUTORIAL_SHARED_PATH "/result.off";

	// a nice user setting, for my configuration!
	const char* viewMesh = TUTORIAL_SHARED_PATH "/planexy.off";

	const char* pipelineScan1File = TUTORIAL_SHARED_PATH "/planexy.off";
	const char* pipelineScan2File = TUTORIAL_SHARED_PATH "/mirrorPlane.off";

}
#endif




