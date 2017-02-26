#ifndef GLOBS_DEF_H
#define GLOBS_DEF_H

#include "tutorial_shared_path.h"

// note :: most of these should be a bunch of consts
// and you should also describe what they are

// Avoid the use of magic numbers and multiple defines.
// Put all the globally-used numbers in here.
namespace GLOBAL
{
//	const char* remeshInputFile = TUTORIAL_SHARED_PATH "/pig.off";
//	const char* remeshOutputFile = TUTORIAL_SHARED_PATH "/pigRemesh.off";
	extern const char* remeshInputFile;
    extern const char* remeshOutputFile;

	extern const char* offsetGenScan1File;
	extern const char* offsetGenScan2File;
}
#endif




