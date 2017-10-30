// *** NOTE :: the # bndry vertices differs in both scans. You must get the cardinality here, for proper index starts.
// my lib includes
#include "interpSurface.h"
#include "remesh.h"
#include "glob_defs.h"

// libIGL includes [ note :: can put in .h file instead ]
#include <igl/writeOFF.h>
#include <igl/writeOBJ.h>
#include <igl/readOFF.h>
#include <igl/boundary_loop.h>
#include <igl/slice.h>

// namespace includes
using namespace Eigen;  
using namespace igl;
using namespace std;

#include <algorithm> // needed for std-lib [containers + iterators]

int main(int argc, char *argv[])
{
	
	// LOAD IN MESHES, COMMENCE TESTS
	// std::cout << "Executing (TEST) for offset surface generation." << std::endl;

	struct Mesh
	{
		Eigen::MatrixXd V; 
		Eigen::MatrixXi F;
	} scan1, scan2, interp, remeshed;

	Eigen::MatrixXd boundaryOne;
	Eigen::MatrixXd boundaryTwo;
	int bvOne_num;
	int bvTwo_num;

	if(!readOFF(GLOBAL::interpSurfGenScan1File,scan1.V,scan1.F)) {
		std::cout << "Failed to load partial scan one." << std::endl;
	} 

	if(!readOFF(GLOBAL::interpSurfGenScan2File,scan2.V,scan2.F)) {
		std::cout << "Failed to load partial scan two." << std::endl;
	}

	// generate interp surf, and extract boundry verts from both scans.
	INTERP_SURF::generateInterpolatingSurface(scan1.V,scan1.F,scan2.V,scan2.F, interp.V,interp.F);
	// generate remeshed surface
	double rel = REMESH::avgEdgeLenInputMeshes(scan1.V, scan1.F, scan2.V, scan2.F);
	rel /= 2.0;
	std::cout << "Remesh avg edge len = [" << rel << "]\n";
	REMESH::remeshSurface(interp.V, interp.F,remeshed.V,remeshed.F, rel);

	// INCORPORATE THIS CODE BLOCK INTO YOUR ALGO PIPELINE NOW!
	INTERP_SURF::generateBoundaryVertices(scan1.V, scan1.F,boundaryOne); 
	INTERP_SURF::generateBoundaryVertices(scan2.V, scan2.F,boundaryTwo); 
	bvOne_num = boundaryOne.rows();
	bvTwo_num = boundaryTwo.rows();
	
	// get REMESHED mesh boundary vertices, for both scan parts
	Eigen::MatrixXd remeshBoundaryOne;
	Eigen::MatrixXd remeshBoundaryTwo;
	REMESH::getRemeshedBoundaryVerts(bvOne_num, bvTwo_num, remeshed.V, remeshed.F, remeshBoundaryOne, remeshBoundaryTwo);


	// oh wait a sec, I can test normals here too! just print out the normal vector field, and see if it matches what you expect :-)

}
