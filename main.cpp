/*** FILE SHOULD BE CALLED :: test_J_Align.cpp ***/
// C++ includes
#include <iostream>
#include <fstream>

// MY LIBRARIES  
#include "glob_defs.h"
#include "meanCurvatureFlow.h"
#include "interpSurface.h"
#include "remesh.h"
#include "j_align.h"
#include "helpers.h"

// LibIgl includes
#include <igl/writeOFF.h>
#include <igl/viewer/Viewer.h>

// Namespace includes 
using namespace Eigen;  
using namespace std;
using namespace igl;

// PREALLOCATION
struct Mesh
{
	Eigen::MatrixXd V; Eigen::MatrixXi F;
} scan1,scan2,interp,remeshed, result;

ofstream debugFile;
Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
double prev_energy = std::numeric_limits<double>::max();
int k = 0;

// METHOD HEADERS
int runPipeline();
bool key_down(igl::viewer::Viewer& viewer, unsigned char key, int mod);

// METHOD BODY
int main(int argc, char *argv[])
{
	runPipeline();
}

int runPipeline()
{
	igl::viewer::Viewer viewer;
	std::cout << "Executing Zero-Overlap, Impuslse-Driven, rigid Alignment Pipeline" << std::endl;

	if(!readOFF(GLOBAL::pipelineScan1File,scan1.V,scan1.F)) {
		std::cout << "Failed to load partial scan one." << std::endl;
	} 

	if(!readOFF(GLOBAL::pipelineScan2File,scan2.V,scan2.F)) {
		std::cout << "Failed to load partial scan two." << std::endl;
	}

	std::cout << "Executing pipeline for the following meshes" << std::endl;
	std::cout << "Mesh one [" << GLOBAL::pipelineScan1File << "]" << std::endl;
	std::cout << "Mesh two [" << GLOBAL::pipelineScan2File << "]" << std::endl;

	std::cout << "Opening debug file" << std::endl;
	debugFile.open(GLOBAL::pipelineDebugFile);

	igl::cat(1,scan1.V,scan2.V,result.V);
	igl::cat(1,scan1.F, MatrixXi(scan2.F.array() + scan1.V.rows()), result.F);
	viewer.data.clear();
	viewer.data.set_mesh(result.V,result.F);
	viewer.callback_key_down = key_down;
	std::cout<<"Press [space] to smooth."<<std::endl;
	std::cout<<"Press [r] to reset."<<std::endl;
	return viewer.launch();
}

bool key_down(igl::viewer::Viewer& viewer, unsigned char key, int mod)
{
	switch(key)
	{
		case 'r':
		case 'R':
			k = 0;
			T = Eigen::Matrix4d::Identity();
			break;
		case ' ':
		{
			std::cout << "Running J_Align for [" << k << "] th iteration" << std::endl;;
			debugFile << "Running J_ALIGN for [" << k << "] th iteration.\n";

			cout << "Applying pipeline, to transition matrix" << endl;
			cout << T << endl;

			//cout << "Apply Rigid Transformation to Scan 1 Mesh" << endl;
			Eigen::MatrixXd transScan1;	
			HELPER::applyRigidTransformation(scan1.V,curT,transScan1); 
			//cout << "Applied Transformation" << endl;

			//cout << "Generating interpolating surface" << endl;
			INTERP_SURF::generateOffsetSurface(transScan1,scan1.F,scan2.V,scan2.F, interp.V,interp.F);
			double rEL = REMESH::avgEdgeLenInputMeshes(transScan1,scan1.F,scan2.V,scan2.F);
			bool remSucc = REMESH::remeshSurface(interp.V,interp.F,remeshed.V,remeshed.F, rEL);
			if(!remSucc)
			{
				std::cout << "BAD REMESHING!" << std::endl;
				std::cout << "Printing out interpolating surface" << std::endl;
				std::string output_bad_mesh = "badInterpSurf[" + std::to_string(k) + "].off";
				std::string scan1_bad = "badScan1.off";
				std::string scan2_bad = "badScan2.off";
				igl::writeOFF(output_bad_mesh, interp.V, interp.F);
				igl::writeOFF(scan1_bad, transScan1, scan1.F);
				igl::writeOFF(scan2_bad, scan2.V, scan2.F);
				exit(0);
			}

//			cout << "Flowing the interpolating surface" << endl;
			Eigen::MatrixXd Vc;
			assert(MCF::meshHasBoundary(remeshed.V,remeshed.F));
			MCF::computeMeanCurvatureFlow(remeshed.V,remeshed.F,0.001, Vc);

			// CHECK IF any transformation matrices do better than the existing one
			// COMPARE their corresponding energy values
			double local_energy = SGD::calculateSurfaceEnergy(Vc,remeshed.F);

			// this energy calculation seems awfully FISHY!
			std::cout << "Prev energy = [" << prev_energy << "]" << std::endl;
			std::cout << "Local energy = [" << local_energy << "]" << std::endl;
			debugFile << "Prev energy = [" << prev_energy << "].\n"; 
			debugFile << "Local energy = [" << local_energy << "].\n"; 

//////////////////////////////////////////////////////////////////////////////////////////////	

			if(local_energy < 0.1) // #TODO :: find a better approach here??
			{
				std::cout << "Impulse based alignment Converged to a solution" << std::endl;
				std::cout << "Printing optimal transformation matrix" << std::endl;	
				std::cout << T << std::endl;
				debugFile << "SGD Converged to a solution.\n";
				debugFile << "Printing optimal transformation matrix.\n";
				debugFile << T << std::endl;
				std::cout << "Closing debug file, exiting program execution.\n";
				// close file
				debugFile.close();	
				// exit program
				//exit(0);
			}
			else
			{
				T = T_local;
				prev_energy = local_energy;
			}
			
			// CREATE one global mesh, contain the two inputs
			Eigen::MatrixXd transScan1;	
			HELPER::applyRigidTransformation(scan1.V,T,transScan1);
			igl::cat(1,transScan1,scan2.V,result.V);
			igl::cat(1,scan1.F, MatrixXi(scan2.F.array() + scan1.V.rows()), result.F);
			viewer.data.clear();
			viewer.data.set_mesh(result.V,result.F);

			++k;
			break;
		}
		default:
			return false;	
	}
	return true;
}


