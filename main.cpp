/*** FILE SHOULD BE CALLED :: test_J_Align.cpp ***/
// C++ includes
#include <iostream>
#include <fstream>

// MY LIBRARIES  
#include "glob_defs.h"
#include "meanCurvatureFlow.h"
#include "interpSurface.h"
#include "remesh.h"
#include "sgd.h" // get rid of later
#include "j_align.h"
#include "rigidbodytemplate.h"
#include "rigidbodyinstance.h"
#include "helpers.h"

// LibIgl includes
#include <igl/writeOFF.h>
#include <igl/viewer/Viewer.h>
// I have a good feeling that it is the vertex normals, versus the edge normals, that are needed in this situation! 
#include <igl/per_edge_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/print_vector.h>
#include <igl/signed_distance.h> 
#include <igl/boundary_loop.h>
#include <igl/slice.h>

// remember to set each rigid body instance's {c,theta} and {cvel,w} to 0 
	// ... nvm, {c,theta} zeroed out @ constructor time!
// #TODO :: for debug output, use print_vector
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
	Eigen::MatrixXd bV;  // boundary vertices
	Eigen::MatrixXi F;
	Eigen::VectorXi bI; // boundary indices; just indxs to specific faces ( not a 1,0 boolean thing ) 
	Eigen::MatrixXd N;
	int numBV; 			// number of boundary vertices
} scan1,scan2,scan1_orig,scan2_orig,interp,remeshed, result;

ofstream debugFile;
Eigen::Matrix4d T1 = Eigen::Matrix4d::Identity();
Eigen::Matrix4d T2 = Eigen::Matrix4d::Identity();
double prev_energy = std::numeric_limits<double>::max();
int iters = 0;
const Vector3d scanOneCenter(0,0,0);
const Vector3d scanOneOrient(0,0,0);
const Vector3d scanTwoCenter(0,0,0);
const Vector3d scanTwoOrient(0,0,0);
const double rho = 1;
const double h = 0.1;

// FOR impulse based alignment scheme
RigidBodyInstance *ScanOneBody;
RigidBodyInstance *ScanTwoBody;
RigidBodyTemplate *ScanOneTemplate;
RigidBodyTemplate *ScanTwoTemplate;

// TBH, your configurations ... lie in the RigidBodyInstances themselves. You update here!
// How do I get the points, on the mesh, that corresponds to the face normals? Do I need to perform barycentric interp here, to find them! ... this does not seem right at all!
// also ... for the code, I need to account for weighting ( based on edge len )! BE CAREFUL!
// ... perhaps the usefulness of the signed distance field?? 

// METHOD HEADERS
int runPipeline();
bool key_down(igl::viewer::Viewer& viewer, unsigned char key, int mod);
void solve_config_deltas(Eigen::Vector3d& delta_cvel, Eigen::Vector3d& delta_omega ,struct Mesh m);
void updateConfiguration(RigidBodyInstance* ScanBody, const Eigen::Vector3d& delta_cVel,const Eigen::Vector3d& delta_omega);

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
	readOFF(GLOBAL::pipelineScan1File,scan1_orig.V,scan1_orig.F);

	if(!readOFF(GLOBAL::pipelineScan2File,scan2.V,scan2.F)) {
		std::cout << "Failed to load partial scan two." << std::endl;
		return -1;
	}
	readOFF(GLOBAL::pipelineScan2File,scan2_orig.V,scan2_orig.F);


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


	// SHIT ... I still need the initial configuratino ( to capture changes; gaaah! ) 
	// SETUP templates and body instances
	// Initial Configuration vector values also set up here too!
	std::cout << "Setting up template and body instances.\n"; 
	ScanOneTemplate = new RigidBodyTemplate(scan1.V,scan1.F);
	ScanTwoTemplate = new RigidBodyTemplate(scan2.V,scan2.F);

	ScanOneBody = new RigidBodyInstance(*ScanOneTemplate,scanOneCenter,scanOneOrient);
	ScanTwoBody = new RigidBodyInstance(*ScanTwoTemplate,scanTwoCenter,scanTwoOrient);

	// any debug output here btw?
	std::cout << "Finished setting up template and body instances.\n"; 

	return viewer.launch();
}

bool key_down(igl::viewer::Viewer& viewer, unsigned char key, int mod)
{
	switch(key)
	{
		case 'r':
		case 'R':
			iters = 0;
			// restore mesh data
			scan1.V = scan1_orig.V;
			scan1.F = scan1_orig.F;
			scan2.V = scan2_orig.V;
			scan2.F = scan2_orig.F;
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
			INTERP_SURF::generateOffsetSurface(scan1.V,scan1.F,scan2.V,scan2.F, interp.V,interp.F);
			double rEL = REMESH::avgEdgeLenInputMeshes(scan1.V,scan1.F,scan2.V,scan2.F);
			bool remSucc = REMESH::remeshSurface(interp.V,interp.F,remeshed.V,remeshed.F, rEL);
			if(!remSucc)
			{
				std::cout << "BAD REMESHING!" << std::endl;
				std::cout << "Printing out interpolating surface" << std::endl;
				std::string output_bad_mesh = "badInterpSurf[" + std::to_string(iters) + "].off";
				std::string scan1_bad = "badScan1.off";
				std::string scan2_bad = "badScan2.off";
				igl::writeOFF(output_bad_mesh, interp.V, interp.F);
				//igl::writeOFF(scan1_bad, transScan1, scan1.F);
				//igl::writeOFF(scan2_bad, scan2.V, scan2.F);
				exit(0);
			}

			// APPLY MCF TO SURFACE
//			cout << "Flowing the interpolating surface" << endl;
			Eigen::MatrixXd Vc;
			assert(MCF::meshHasBoundary(remeshed.V,remeshed.F));
			MCF::computeMeanCurvatureFlow(remeshed.V,remeshed.F,0.001, Vc);

			// CHECK IF impulses improved alignment
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
				//exit(0);
			}
			else
			{
				prev_energy = local_energy;
			}

			/*******************************
			 *** IMPULSE BASED ALIGNMENT ***
			 *** #TODO :: modularize     ***
			 *******************************/
			// UPDATE the configuration
			// REWRITE information to scan1.V,scan2.V

			// compute vertex normals for both scans
			igl::PerVertexNormalsWeightingType weighting = PER_VERTEX_NORMALS_WEIGHTING_TYPE_UNIFORM; 
			igl::per_vertex_normals(scan1.V,scan1.F,weighting,scan1.N);
			igl::per_vertex_normals(scan2.V,scan2.F,weighting,scan2.N);

			// CALCULATE boundary vertices
			igl::boundary_loop(scan1.F,scan1.bI);
			igl::slice(scan1.V,scan1.bI,1,scan1.bV);
			scan1.numBV = scan1.bV.rows();
				
			igl::boundary_loop(scan2.F,scan2.bI);
			igl::slice(scan2.V,scan2.bI,1,scan2.bV);
			scan2.numBV = scan2.bV.rows();

			// we can just assume \rho Vol = 1 ( unit mass, really ) ... or modify this to our wish

			std::cout << "size of boundary loop, scan 1 = [" << scan1.numBV << "]\n";
			std::cout << "size of boundary loop, scan 2 = [" << scan2.numBV << "]\n";

			// ITERATE over boundary verts; calculate interpolating bnd {norm, vert}
			// you also already possess edge lengths ... difference in vertex position data ( take the norm here ) 
			// you also have the bndry indices ... bemember, this covers all vertices, right? 

			// SOLVE for delta_cVel, delta_omega, for both scans
			Eigen::Vector3d delta_cVel_1 = Eigen::Vector3d(0,0,0);
			Eigen::Vector3d delta_omega_1 = Eigen::Vector3d(0,0,0);
			Eigen::Vector3d delta_cVel_2 = Eigen::Vector3d(0,0,0);
			Eigen::Vector3d delta_omega_2 = Eigen::Vector3d(0,0,0);

			solve_config_deltas(delta_cVel_1,delta_omega_1,scan1);
			solve_config_deltas(delta_cVel_2,delta_omega_2,scan2);

			updateConfiguration(ScanOneBody,delta_cVel_1,delta_omega_1);
			updateConfiguration(ScanTwoBody,delta_cVel_2,delta_omega_2);

			// GENERATE transformation ( rotation + translation ) components
			// for both partial scan boundaries
			Eigen:Vector3d t1_comp = ScanOneBody->c - ScanOneBody->c_0;
			Eigen::Matrix3d r1_comp = Eigen::Matrix3d::Identity(); // #TODO :: get actual rotational part! 
			T1.block<3,1>(0,3) = t1_comp; // yes, thsi block operation works as expected
			T1.block<3,3>(0,0) = r1_comp;
			std::cout << "Transformation matrix, for scan one [t1] is:\n";
			std::cout << T1 << "\n";

			Eigen::Vector3d t2_comp = Eigen::Vector3d(0,0,0);
			Eigen::Matrix3d r2_comp = Eigen::Matrix3d::Identity();
			T2.block<3,1>(0,3) = t2_comp;
			T2.block<3,3>(0,0) = r2_comp;

			// note :: since normals update with scan1.v data ... must update scan1.v @ each timestep, with respect to scan1_orig.v, too!

			/*******************************************************
			 *** HOW TO I GET THIS TRANSFORMATION MATRIX THOUGH? ***
			 ******************************************************/
			// well, it's easy to determine the translation  :: c - c_0
			// For the rotational piece :: set [theta - theta_0] to axis-angle formula
			// 		and generate your 3d rotation matrix ( refer to Vouga's code ) 
			// CREATE one global mesh, containing the two inputs, perturbed by (c - c_0)
			// wait ... isn't it the previous configuratino, to store here? 
			// well, it depends! Are you creating an animated, or non-animated version? For now, let us focus on the animated version! ... actually, as long as you use the initial mesh ONLY ( the template coordinates ) ... this shouldn't be that bad!
	
			// CREATE global mesh, contain the two inputs aligned based on impulses
			HELPER::applyRigidTransformation(scan1_orig.V,T1,scan1.V);
			HELPER::applyRigidTransformation(scan2_orig.V,T2,scan2.V);
			igl::cat(1,scan1.V,scan2.V,result.V);
			igl::cat(1,scan1.F, MatrixXi(scan2.F.array() + scan1.V.rows()), result.F);
			viewer.data.clear();
			viewer.data.set_mesh(result.V,result.F);

			++iters;
			break;
		}
		default:
			return false;	
	}
	return true;
}

void solve_config_deltas(Eigen::Vector3d& delta_cvel, Eigen::Vector3d& delta_omega ,struct Mesh m)
{
	for(int k = 0; k < m.numBV; ++k) // iterate through SIZE of bV
	{
		// get idx of 2 bndry verts of interest	
		int i = k % m.numBV;
		int j = (k + 1) % m.numBV;

		Eigen::Vector3d v_i = m.bV.row(i);
		Eigen::Vector3d v_j = m.bV.row(j);

		Eigen::Vector3d n_i = m.N.row(i);
		Eigen::Vector3d n_j = m.N.row(j);

		// average two normals
		Eigen::Vector3d n_avg = 0.5 * (n_i + n_j);
	
		// Rescale F_k_ext, by edge length, and solve for J_k_ext
		n_avg.normalize();
		double edge_len = (v_j - v_i).norm();
		Eigen::Vector3d n_avg_rescaled = edge_len * n_avg;

		Eigen::Vector3d F_k_ext = n_avg_rescaled; 
		Eigen::Vector3d J_k_ext = h * F_k_ext;

		// add up contributions to delta_cVel
		Eigen::Vector3d delta_cVel_k = J_k_ext;
		delta_cvel += delta_cVel_k;

		// add up contributions to delta_omega
		Eigen::Vector3d delta_omega_k = Eigen::Vector3d(0,0,0);
		delta_omega += delta_omega_k;
	}
}

void updateConfiguration(RigidBodyInstance* ScanBody, const Eigen::Vector3d& delta_cVel,const Eigen::Vector3d& delta_omega)
{
	ScanBody->cvel += delta_cVel;
	ScanBody->w += delta_omega; // #TODO :: fix this calculation!
	ScanBody->c += (ScanBody->cvel) * h; 
	ScanBody->theta += (ScanBody->w) * h;
	// vanish out the velocities ( 'molasses' set up ) 
	ScanBody->cvel.setZero();
	ScanBody->w.setZero();
}


// well, the output in the impulse based alignment scheme is two independent transformation matrices
// #TODO :: seperate out the impulse based alignment methods + code, from the pipeline code, itself! 
