// MY LIBRARIES  
// #TODO :: cut down. This seriously affects compilation time.
#include "glob_defs.h"
#include "meanCurvatureFlow.h"
#include "interpSurface.h"
#include "remesh.h"
#include "sgd.h"  // #TODO :: eliminate later
#include "j_align.h"
#include "rigidbodytemplate.h"
#include "rigidbodyinstance.h"
#include "helpers.h"
#include "vectormath.h"

// LibIgl includes
#include <igl/readOFF.h>
#include <igl/writeOFF.h>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/per_edge_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/print_vector.h>
#include <igl/signed_distance.h> 
#include <igl/boundary_loop.h>
#include <igl/slice.h>
#include <igl/random_dir.h>
#include <igl/boundary_loop.h>
#include <igl/per_face_normals.h>

// C++ includes
#include <iostream>
#include <fstream>
#include <algorithm> 			// needed for std-lib [containers + iterators]

// viewer logic [ For later debugging ... crucial only for main method ]
#include <igl/viewer/Viewer.h>
#include <igl/cat.h>

// Namespace includes 
using namespace Eigen;  
using namespace std;
using namespace igl;

// Forward declaration of classes
class RigidBodyTemplate;
class RigidBodyInstance;

// PREALLOCATION
struct Mesh
{
	Eigen::MatrixXd V; 
	Eigen::MatrixXi F;

	Eigen::MatrixXd bV;  // boundary vertices
	Eigen::VectorXi bI; // boundary indices; just indxs to specific faces ( not a 1,0 boolean thing ) 

	Eigen::MatrixXd N;

	int numBV; 			// number of boundary vertices

} scan1,scan2,scan1_rest,scan2_rest,interp,remeshed, flowed, debug_result, result, normals, normals_info;
// #TODO :: note that these <normals> are for <scan1> primarily, and <normals_info> is used to compute said <normals> 3D triangle-mesh structure.
// #TODO :: please improve naming conventions here.
// #TODO :: organize file better. I'm getting lost in these variable names and what not.
// #TODO :: note the stage in which speciifc data is written ( i.e. interp, remesh, flow ) 

ofstream debugFile;
Eigen::Matrix4d T1 = Eigen::Matrix4d::Identity();
Eigen::Matrix4d T2 = Eigen::Matrix4d::Identity();
// I expect initialy COMS and orientations to be equal to 0. 
const Vector3d scanOneCenter(0,0,0);
const Vector3d scanOneOrient(0,0,0);
const Vector3d scanTwoCenter(0,0,0);
const Vector3d scanTwoOrient(0,0,0);
double prev_energy = std::numeric_limits<double>::max();
int iters = 0;
const double h = 0.01;

// FOR impulse based alignment scheme
RigidBodyInstance *ScanOneBody;
RigidBodyInstance *ScanTwoBody;
RigidBodyTemplate *ScanOneTemplate;
RigidBodyTemplate *ScanTwoTemplate;

// METHOD HEADERS
int runPipeline();
bool key_down(igl::viewer::Viewer& viewer, unsigned char key, int mod);

// METHOD BODY
int main(int argc, char *argv[])
{
	// EXECUTE J_ALIGN PIPELINE
	runPipeline();
}

// note :: only 2 rigid bodies to consider ... no point in putting them in an iterator
// FOCUS on getting one of them to work atm
int runPipeline() 
{
	// #TODO :: work on setting up the visualization correctly!
	igl::viewer::Viewer viewer;
	std::cout << "Executing Zero-Overlap, Impuslse-Driven, rigid Alignment Pipeline" << std::endl;

	if(!readOFF(GLOBAL::pipelineScan1File,scan1.V,scan1.F)) {
		std::cout << "Failed to load partial scan one." << std::endl;
		return -1;
	} 
	readOFF(GLOBAL::pipelineScan1File,scan1_rest.V,scan1_rest.F);

	if(!readOFF(GLOBAL::pipelineScan2File,scan2.V,scan2.F)) {
		std::cout << "Failed to load partial scan two." << std::endl;
		return -1;
	}
	readOFF(GLOBAL::pipelineScan2File,scan2_rest.V,scan2_rest.F);

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

	// SETUP templates and body instances
	// Initial Configuration vector values also set up here too!
	std::cout << "Setting up template and body instances.\n"; 
	ScanOneTemplate = new RigidBodyTemplate(scan1.V,scan1.F);
	ScanTwoTemplate = new RigidBodyTemplate(scan2.V,scan2.F);

	ScanOneBody = new RigidBodyInstance(*ScanOneTemplate,scanOneCenter,scanOneOrient);
	ScanTwoBody = new RigidBodyInstance(*ScanTwoTemplate,scanTwoCenter,scanTwoOrient);

	// any debug output here btw?
	std::cout << "Finished setting up template and body instances.\n"; 

	viewer.launch();
	return 0;
}

bool key_down(igl::viewer::Viewer& viewer, unsigned char key, int mod)
{
	switch(key)
	{
		case 'r':
		case 'R':
			iters = 0;
			// restore mesh data
			scan1.V = scan1_rest.V;
			scan1.F = scan1_rest.F;
			scan2.V = scan2_rest.V;
			scan2.F = scan2_rest.F;
			// restore templates 
			ScanOneTemplate = new RigidBodyTemplate(scan1.V,scan1.F);
			ScanTwoTemplate = new RigidBodyTemplate(scan2.V,scan2.F);
			ScanOneBody = new RigidBodyInstance(*ScanOneTemplate,scanOneCenter,scanOneOrient);
			ScanTwoBody = new RigidBodyInstance(*ScanTwoTemplate,scanTwoCenter,scanTwoOrient);
			break;
		case ' ':
		{
			std::cout << "Running J_Align for [" << iters << "] th iteration" << std::endl;;
			debugFile << "Running J_ALIGN for [" << iters << "] th iteration.\n";

			// GENERATE INTERPOLATING SURFACE
			//cout << "Generating interpolating surface" << endl;
			INTERP_SURF::generateInterpolatingSurface(scan1.V,scan1.F,scan2.V,scan2.F, interp.V,interp.F);
			double rEL = REMESH::avgEdgeLenInputMeshes(scan1.V,scan1.F,scan2.V,scan2.F);
			rEL = rEL / 3.0; // #NOTE :: can actually visualize pinching occuring here
			bool remSucc = REMESH::remeshSurface(interp.V,interp.F,remeshed.V,remeshed.F, rEL);

			if(!remSucc)
			{
				std::cout << "BAD REMESHING!" << std::endl;
				std::cout << "Printing out interpolating surface" << std::endl;
				std::string output_bad_mesh = "badInterpSurf[" + std::to_string(iters) + "].off";
				std::string scan1_bad = "badScan1.off";
				std::string scan2_bad = "badScan2.off";
				igl::writeOFF(output_bad_mesh, interp.V, interp.F);
				exit(0);
			}

			// APPLY MCF TO SURFACE [ denoted as 'flowed' here ]
			assert(MCF::meshHasBoundary(remeshed.V,remeshed.F));
			MCF::computeMeanCurvatureFlow(remeshed.V,remeshed.F,0.001, flowed.V);
			flowed.F = remeshed.F;

			// CHECK IF impulses improved alignment
			double local_energy = SGD::calculateSurfaceEnergy(flowed.V,flowed.F);

			// this energy calculation seems awfully FISHY!
			std::cout << "Prev energy = [" << prev_energy << "]" << std::endl;
			std::cout << "Local energy = [" << local_energy << "]" << std::endl;
			debugFile << "Prev energy = [" << prev_energy << "].\n"; 
			debugFile << "Local energy = [" << local_energy << "].\n"; 

			///////////////////////////////////////////////////////////////////////////	
			// something tells me that I need to be careufl with this local energy value
			// and need to double check magnitudes of F_ext being used here!
			// [T1,T2] need to be properly updated here!
			if(local_energy < 0.1) // #TODO :: find a better approach here??
			{
				std::cout << "Impulse based alignment Converged to a solution" << std::endl;
				std::cout << "Printing optimal transformation matrices" << std::endl;	
				std::cout << "Transformation matrix, for scan one [t1] is:\n";
				std::cout << T1 << "\n";
				std::cout << "Transformation matrix, for scan one [t2] is:\n";
				std::cout << T2 << "\n";
				debugFile << "SGD Converged to a solution.\n";
				debugFile << "Printing optimal transformation matrix.\n";
				debugFile << "Transformation matrix, for scan one [t1] is:\n";
				debugFile << T1 << "\n";
				debugFile << "Transformation matrix, for scan one [t2] is:\n";
				debugFile << T2 << "\n";
				std::cout << "Closing debug file, exiting program execution.\n";
				// close file
				debugFile.close();	
				// exit program
			}
			else
			{
				prev_energy = local_energy;
			}

			/*******************************
			 *** IMPULSE BASED ALIGNMENT ***
			 *******************************/
			// [1] CALCULATE boundary vertices for both scans
			// generate boundary vertices from interpolating surface
			Eigen::MatrixXd boundaryOne;
			Eigen::MatrixXd boundaryTwo;

			INTERP_SURF::generateBoundaryVertices(scan1.V, scan1.F,boundaryOne); 
			INTERP_SURF::generateBoundaryVertices(scan2.V, scan2.F,boundaryTwo); 
			int bvOne_num = boundaryOne.rows();
			int bvTwo_num = boundaryTwo.rows();
			
			// get FLOWED mesh ( after remeshing too, btw ) boundary vertices, for both scan parts
			std::vector<int> remeshBoundaryOne;
			std::vector<int> remeshBoundaryTwo;
			REMESH::getRemeshedBoundaryVerts(bvOne_num, bvTwo_num, flowed.V, flowed.F, remeshBoundaryOne, remeshBoundaryTwo);
			cout << "Done collecting remeshed boundary vertices" << endl;

			// [2] SOLVE for delta_cVel, delta_omega, for both scans
			// NOTE :: Can assume $\rho$ Vol = 1 [ that is, our setting is of const unit mass ]
			Eigen::Vector3d delta_cVel_1 = Eigen::Vector3d(0,0,0);
			Eigen::Vector3d delta_omega_1 = Eigen::Vector3d(0,0,0);

			Eigen::Vector3d delta_cVel_2 = Eigen::Vector3d(0,0,0);
			Eigen::Vector3d delta_omega_2 = Eigen::Vector3d(0,0,0);

			// [3] SOLVE for external forces assoc w/boundary 1
			Eigen::MatrixXd boundaryV_1;
			Eigen::MatrixXd extBndForceMat_1;
			J_ALIGN::solveExternalForcesForBoundary(remeshBoundaryOne, flowed.V, flowed.F, boundaryV_1, extBndForceMat_1);  

			// [4] SOLVE for external forces assoc w/boundary 2 
			Eigen::MatrixXd boundaryV_2;
			Eigen::MatrixXd extBndForceMat_2;
			J_ALIGN::solveExternalForcesForBoundary(remeshBoundaryTwo, flowed.V, flowed.F, boundaryV_2, extBndForceMat_2);  

			// [5] solve for config_deltas
			J_ALIGN::solve_config_deltas(boundaryV_1, extBndForceMat_1, ScanOneBody, delta_cVel_1,delta_omega_1);
			J_ALIGN::solve_config_deltas(boundaryV_2, extBndForceMat_2, ScanTwoBody, delta_cVel_2,delta_omega_2);

			// [6] update rigid body instance configurations
			J_ALIGN::updateConfiguration(ScanOneBody,delta_cVel_1,delta_omega_1);
			J_ALIGN::updateConfiguration(ScanTwoBody,delta_cVel_2,delta_omega_2);

			// [7] GENERATE transformation [ rot + trans ] components
			// ##TODO :: fix this transformation. It isn't correct. Something about it is off apparantely.
			Eigen::Matrix4d T1 = Eigen::Matrix4d::Identity();
			Eigen::Matrix4d T2 = Eigen::Matrix4d::Identity();
			J_ALIGN::solveTransformation(ScanOneBody, T1);
			J_ALIGN::solveTransformation(ScanTwoBody, T2);

			// [7] APPLY rigid transformations to both boundaries.
			HELPER::applyRigidTransformation(scan1_rest.V,T1,scan1.V);
			HELPER::applyRigidTransformation(scan2_rest.V,T2,scan2.V);
			++iters;

			
			// [8] CREATE global mesh, contain the two inputs aligned based on impulses
			igl::cat(1,scan1.V,scan2.V,result.V);
			igl::cat(1,scan1.F, MatrixXi(scan2.F.array() + scan1.V.rows()), result.F);
			viewer.data.clear();
			viewer.data.set_mesh(result.V,result.F);
			break;
		}
		default:
			return false;	
	}
	return true;
}
