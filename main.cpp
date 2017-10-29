// my lib includes
#include "interpSurface.h"
#include "remesh.h"
#include "glob_defs.h"

// libIGL includes [ note :: can put in .h file instead ]
#include <igl/writeOFF.h>
#include <igl/writeOBJ.h>
#include <igl/readOFF.h>
#include <igl/boundary_loop.h>

// namespace includes
using namespace Eigen;  
using namespace igl;
using namespace std;

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

	if(!readOFF(GLOBAL::interpSurfGenScan1File,scan1.V,scan1.F)) {
		std::cout << "Failed to load partial scan one." << std::endl;
	} 

	if(!readOFF(GLOBAL::interpSurfGenScan2File,scan2.V,scan2.F)) {
		std::cout << "Failed to load partial scan two." << std::endl;
	}

	// generate interp surf, and extract boundry verts from both scans.
	INTERP_SURF::generateInterpolatingSurface(scan1.V,scan1.F,scan2.V,scan2.F, interp.V,interp.F);
	INTERP_SURF::generateBoundaryVertices(scan1.V, scan1.F,boundaryOne); 
	INTERP_SURF::generateBoundaryVertices(scan2.V, scan2.F,boundaryTwo); 

	// generate remeshed surface
	double rel = REMESH::avgEdgeLenInputMeshes(scan1.V, scan1.F, scan2.V, scan2.F);
	rel /= 2.0;
	std::cout << "Remesh avg edge len = [" << rel << "]\n";
	REMESH::remeshSurface(interp.V, interp.F,remeshed.V,remeshed.F, rel);

	// write OFF [interp, remeshed] surfaces
	igl::writeOFF("interpSurf.off", interp.V, interp.F);
	igl::writeOFF("remeshedSurf.off", remeshed.V, remeshed.F);

	// ITERATE over resmeshed.V vertices. 
	// KNOW that the first [boundaryOne,boundaryTwo] vertices are already set.
	// LET us construct boundary loops for the remeshsed surface and see what we get as the results
	Eigen::MatrixXd remeshBndry;
	std::vector<std::vector<int>> bndryLoop; // #NOTE :: thsi is the appropriate call
	// Refer to <hessian_energy.cpp> for your example.

	// GET boundary_loops for the remeshed vertex
	// wait a sec ... as long as I know one of these indices ( in remeshed loop ) is in the ( interp ) one too ... I'm good :-) SHIT!!
	// ... this makes my life so easy! OMG I don't have to go through Vouga's solution now! AAHHHH!
	igl::boundary_loop(remeshed.F, bndryLoop);
	for( const std::vector<int> &v : bndryLoop )
	{
		for ( int x : v )
		{
			std::cout << x << ' ';
		}
		cout << endl;
		cout << endl;
		cout << endl;
	}

	//igl::slice(remeshed.V, remeshBndIndexes, 1, remeshBndry);
	//cout << remeshBndry << endl;
	
	

	int startIdx = boundaryOne.rows() + boundaryTwo.rows();
	for(int i = startIdx; i < remeshed.V.rows(); ++i)
	{
	}
	

}

