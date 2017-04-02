// MY LIBRARIES  
#include "glob_defs.h"
#include "meanCurvatureFlow.h"
#include "interpSurface.h"
#include "remesh.h"
#include "sgd.h"
#include "helpers.h"

// LibIgl includes
#include <igl/writeOFF.h>
#include <igl/viewer/Viewer.h>

using namespace Eigen;  
using namespace std;
using namespace igl;

struct Mesh
{
	Eigen::MatrixXd V; Eigen::MatrixXi F;
} scan1,scan2,interp,remeshed, result;

int runPipeline();
bool key_down(igl::viewer::Viewer& viewer, unsigned char key, int mod);

// Develop an initial transformation matrix, filled with random data 
std::vector<Eigen::Matrix4d> transMats;	


/*
 MatrixXd m(2,2);
  m(0,0) = 3;
  m(1,0) = 2.5;
  m(0,1) = -1;
  m(1,1) = m(1,0) + m(0,1);
*/


Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
//Eigen::Matrix4d T;

int k = 0;
double prev_energy = std::numeric_limits<double>::max();

int main(int argc, char *argv[])
{
// dumb reaosn ... must be in a func!!
// testing out some existing "optimal" transformation matrices ... 
/*
T << 0.905716, 0.0236395, 0.423226,  0.100115,
	-0.055592, 0.996444, 0.0633116, -0.0254154,
	-0.420225, -0.0808702, 0.903809, -0.0150992,
	0		,0			,0			,1;
*/

	runPipeline();
}
/*
{
  //igl::readOFF(TUTORIAL_SHARED_PATH "/bunny.off", V, F);
  //igl::readOFF(GLOBAL::viewMesh, result.V, result.F);

  // Plot the mesh
  //igl::viewer::Viewer viewer;
  //viewer.data.set_mesh(result.V, result.F);
  //viewer.launch();
}
*/

int runPipeline()
{
	igl::viewer::Viewer viewer;
	std::cout << "Executing Zero-Overlap Rigid Alignment Pipeline" << std::endl;

	if(!readOFF(GLOBAL::pipelineScan1File,scan1.V,scan1.F)) {
		std::cout << "Failed to load partial scan one." << std::endl;
	} 

	if(!readOFF(GLOBAL::pipelineScan2File,scan2.V,scan2.F)) {
		std::cout << "Failed to load partial scan two." << std::endl;
	}

	std::cout << "Executing pipeline for the following meshes" << endl;
	std::cout << "Mesh one [" << GLOBAL::pipelineScan1File << "]" << std::endl;
	std::cout << "Mesh two [" << GLOBAL::pipelineScan2File << "]" << std::endl;

	// Develop an initial transformation matrix, filled with random data 
	transMats.push_back(T);

	// initial gradient descent!
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
			transMats.clear();
			transMats.push_back(T);
			k = 0;
			break;
		case ' ':
		{
			cout << "Running SGD for [" << k << "] th iteration" << endl;;
			std::vector<double> energies;
			for(int i = 0; i < transMats.size(); ++i)
			{
//				cout << "Applying pipeline, to [" << i << "th] transition matrix" << endl;
//				cout << transMats[i] << endl;

				//cout << "Apply Rigid Transformation to Scan 1 Mesh" << endl;
				Eigen::MatrixXd transScan1;	
				Eigen::Matrix4d curT = transMats[i];
				HELPER::applyRigidTransformation(scan1.V,curT,transScan1);

				//cout << "Generating interpolating surface" << endl;
				INTERP_SURF::generateOffsetSurface(transScan1,scan1.F,scan2.V,scan2.F, interp.V,interp.F);
				double remeshEdgeLen = REMESH::avgEdgeLenInputMeshes(transScan1,scan1.F,scan2.V,scan2.F);
				bool remSucc = REMESH::remeshSurface(interp.V,interp.F,remeshed.V,remeshed.F, remeshEdgeLen);
				cout << "Remeshed interpolating surface" << endl;
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
				Eigen::MatrixXd V = remeshed.V;
				Eigen::MatrixXi F = remeshed.F;
				bool hasBndry = MCF::meshHasBoundary(V,F);
				assert(hasBndry);
				Eigen::MatrixXd Vc;
				MCF::computeMeanCurvatureFlow(V,F,0.001, Vc);
				double energy = SGD::calculateSurfaceEnergy(Vc,F);
				//cout << "Energy value is [" << energy << "]" << endl;
				energies.push_back(energy);
			}

			// CHECK IF any transformation matrices do better than the existing one
			// COMPARE their corresponding energy values
			Eigen::Matrix4d T_local;

			//for ( auto i : transMats )
			//	std::cout << i << ' ';
			//std::cout << std::endl;


			/*for ( auto i : energies)
				std::cout << i << ' ';
			std::cout << std::endl;
*/

			double local_energy = SGD::findOptimalTransMat(transMats,energies,T_local);

			// this energy calculation seems awfully FISHY!
			std::cout << "Prev energy = [" << prev_energy << "]" << std::endl;
			std::cout << "Local energy = [" << local_energy << "]" << std::endl;
			if(local_energy < prev_energy)
			{
				T = T_local;
				prev_energy = local_energy;
			}
			else
			{
				std::cout << "SGD Converged to a solution" << std::endl;
				std::cout << "Printing optimal transformation matrix" << std::endl;	
				std::cout << T << std::endl;
				exit(0);
			}

			transMats.clear();
			SGD::generateTransMats(T,transMats);
			cout << "OPTIMAL transition matrix is : " << endl;
			cout << T << endl;
			cout << "-----------------------------------------------" << endl;
			energies.clear();

			// CREATE one global mesh, contain the two inputs
			Eigen::MatrixXd transScan1;	
			HELPER::applyRigidTransformation(scan1.V,T,transScan1);
			igl::cat(1,transScan1,scan2.V,result.V);
			igl::cat(1,scan1.F, MatrixXi(scan2.F.array() + scan1.V.rows()), result.F);
			viewer.data.clear();
			viewer.data.set_mesh(result.V,result.F);

			// PRINT out the interpolating surface ... something seems off here!
			std::string output_bad_mesh = "interpSurf[" + std::to_string(k) + "].off";
			igl::writeOFF(output_bad_mesh, interp.V, interp.F);
	
			if(k % 10 == 0)
			{	
				std::string output_of_mesh = "pipelineOutput[" + std::to_string(k) + "].off";
				std::cout << "WRITING PIPELINE data for iteration [" << std::to_string(k) << " ].\n"; 
				igl::writeOFF(output_of_mesh, result.V, result.F);	
			}

			++k;
			break;
		}
		default:
			return false;	
	}
	return true;
}


