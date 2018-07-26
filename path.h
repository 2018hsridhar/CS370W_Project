// path.h :: represents routes to commonly access directories. 
// ---- Akin to <router> logic in MVC apps
// All paths include <..> since we perform external builds, NOT, internal builds
// KEEP tests consolidated into one directory.

#ifndef path_h_included
#define path_h_included

// path to <mesh> directory. Represents commonly accessed meshes. 
#ifndef SHARED_PATH
//#define SHARED_PATH "../MESHES"
//#define SHARED_PATH "../shared"
#endif

#ifndef MESHES_PATH 
#define MESHES_PATH "../MESHES"
#endif


// path to <processing> dir. Represents meshes in the stage of being processed.
#ifndef PROC_PATH
#define PROC_PATH "../PROCESSING"
#endif

/*
// path to <debug> dir. Represents debug output
#ifndef DEBUG_PATH
#define DEBUG_PATH "../DEBUG"
#endif

// path to <testing> dir. Represents test cases being run
#ifndef TEST_PATH
#define TEST_PATH "../TEST"
#endif

// path to <libs> dir. Represents code base
#ifndef TEST_PATH
#define TEST_PATH "../LIB"
#endif

// path to CGAL library
// #TODO :: update this library later. For now, it doesn't need to be updated
#ifndef CGAL_PATH
#define CGAL_Path "/p/CGAL-4.9"
#endif

// path to LibIgl library
#ifndef IGL_PATH 
#define IGL_Path "/p/libigl"
#endif
*/

#endif
