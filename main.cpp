// PURPOSE :: a greedy surface reconstruction algorithm, from 2 range images, based of closest distances between boundary vertices

#include "interpSurface.h"
#include "remesh.h"
#include "glob_defs.h"
using namespace Eigen;  
using namespace igl;

// #TODO :: if time available, convert main() to a testing method
int main(int argc, char *argv[])
{
	// LOAD IN MESHES 
	std::cout << "Executing (TEST) for offset surface generation." << std::endl;

	struct Mesh
	{
		Eigen::MatrixXd V; 
		Eigen::MatrixXi F;
	} scan1,scan2,scene, remeshed;

	if(!readOFF(GLOBAL::offsetGenScan1File,scan1.V,scan1.F)) {
		std::cout << "Failed to load partial scan one." << std::endl;
	} 

	if(!readOFF(GLOBAL::offsetGenScan2File,scan2.V,scan2.F)) {
		std::cout << "Failed to load partial scan two." << std::endl;
	}

	INTERP_SURF::generateOffsetSurface(scan1.V,scan1.F,scan2.V,scan2.F, scene.V,scene.F);
	//REMESH::remeshSurface(scene.V,scene.F,remeshed.V,remeshed.F);

	//  SETUP LibIgl Viewer 
	igl::viewer::Viewer viewer;
    viewer.data.set_mesh(scene.V, scene.F); 
	//viewer.data.set_mesh(remeshed.V, remeshed.F); 
	viewer.launch();
}

