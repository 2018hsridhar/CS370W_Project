﻿// MY FILES
#include "glob_defs.h"
#include "meanCurvatureFlow.h"
#include "interpSurface.h"
#include "remesh.h"

using namespace Eigen;  
using namespace igl;

int main(int argc, char *argv[])
{
	// code to test offset surface gen + remeshing
	// LOAD IN MESHES 
	std::cout << "Executing (TEST) for offset surface generation." << std::endl;
	struct Mesh
	{
		Eigen::MatrixXd V; 
		Eigen::MatrixXi F;
	} scan1,scan2,scene, remeshed;

	if(!readOFF(GLOBAL::interpSurfGenScan1File,scan1.V,scan1.F)) {
		std::cout << "Failed to load partial scan one." << std::endl;
	} 

	if(!readOFF(GLOBAL::interpSurfGenScan2File,scan2.V,scan2.F)) {
		std::cout << "Failed to load partial scan two." << std::endl;
	}

	INTERP_SURF::generateOffsetSurface(scan1.V,scan1.F,scan2.V,scan2.F, scene.V,scene.F);
	REMESH::remeshSurface(scene.V,scene.F,remeshed.V,remeshed.F);

	//  SETUP LibIgl Viewer 
	igl::viewer::Viewer viewer;
	//viewer.data.set_mesh(scene.V, scene.F); 
	viewer.data.set_mesh(remeshed.V, remeshed.F); 
	viewer.launch();

/*
	// code to test in data generation
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    // LOAD mesh data ( OFF format )
    //igl::readOFF(TUTORIAL_SHARED_PATH "/camelhead.off",V,F);
    igl::readOFF(TUTORIAL_SHARED_PATH "/cow.off",V,F);
    bool hasBndry = MCF::meshHasBoundary(V,F);
    if(hasBndry) std::cout << "INPUT Mesh does have a boundary.\n";
	Eigen::MatrixXd Vc;
	MCF::computeMeanCurvatureFlow(V,F,0.001, Vc);
    igl::viewer::Viewer viewer;
	viewer.data.set_mesh(Vc,F);
    viewer.launch();
*/

    return 0;
}

    
