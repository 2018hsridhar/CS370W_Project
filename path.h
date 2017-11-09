// path.h --- represents routes to commonly access directories. 
// it's akin to <router> logic in MVC apps
// all of these paths include <..> since we perform external, vs internal, builds
// also, you need to keep your tests consolidated. This will prove even more useful in the future
// #TODO :: slowly shift your files from <tutroial_shared_path.h> to <path.h>. 

#ifndef path_h_included
#define path_h_included

// path to <mesh> directory. Represents commonly accessed meshes. 
#ifndef SHARED_PATH
//#define SHARED_PATH "../MESHES"
#define SHARED_PATH "../shared"
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
