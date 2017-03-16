// MY LIBRARIES  
#include "glob_defs.h"
#include "meanCurvatureFlow.h"
#include "interpSurface.h"
#include "remesh.h"
#include "sgd.h"

// LibIgl includes
#include <igl/writeOFF.h>

using namespace Eigen;  
using namespace std;
using namespace igl;

// of course you will have an isuse here, since glob_defs.h does not get linked to this program, ONLY to the main program!
int main(int argc, char *argv[])
{
	// code to test offset surface gen + remeshing
	// LOAD IN MESHES 
	std::cout << "Executing (TEST) for offset surface generation." << std::endl;
	struct Mesh
	{
		Eigen::MatrixXd V; 
		Eigen::MatrixXi F;
	} scan1,scan2,interp,remeshed;

	if(!readOFF(GLOBAL::interpSurfGenScan1File,scan1.V,scan1.F)) {
		std::cout << "Failed to load partial scan one." << std::endl;
	} 

	if(!readOFF(GLOBAL::interpSurfGenScan2File,scan2.V,scan2.F)) {
		std::cout << "Failed to load partial scan two." << std::endl;
	}

	INTERP_SURF::generateOffsetSurface(scan1.V,scan1.F,scan2.V,scan2.F, interp.V,interp.F);
	double remeshEdgeLen = REMESH::avgEdgeLenInputMeshes(scan1.V,scan1.F,scan2.V,scan2.F);
	REMESH::remeshSurface(interp.V,interp.F,remeshed.V,remeshed.F, remeshEdgeLen);

	//  SETUP LibIgl Viewer 
/*
	igl::viewer::Viewer viewer;
	Eigen::MatrixXd V = remeshed.V;
    Eigen::MatrixXi F = remeshed.F;
    bool hasBndry = MCF::meshHasBoundary(V,F);
    if(hasBndry) std::cout << "INPUT Mesh does have a boundary.\n";
	Eigen::MatrixXd Vc;
	MCF::computeMeanCurvatureFlow(V,F,0.001, Vc);
	double energy = SGD::calculateSurfaceEnergy(Vc,F);
*/
//	igl::writeOFF(GLOBAL::meanCurvFlowOutputScanFile,Vc,F);

/*
    igl::viewer::Viewer viewer;
	viewer.data.set_mesh(Vc,F);
    viewer.launch();
*/
		// Develop an initial transformation matrix, filled with random data 
		Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
		std::vector<Eigen::Matrix4d> transMats;	
		SGD::generateTransMats(T,transMats);
/*
		for(int i = 0; i < transMats.size(); ++i)
		{
			cout << "i = [" << i << "]" << endl;
			cout << transMats[i] << endl;
			cout << "--------------------------" << endl;
		}
*/
	


    return 0;
}

    
