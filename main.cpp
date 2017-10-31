// #TODO :: refactor file; some of these methods truly do belong elsewhere, and are contributing to clutter.
/*** FILE SHOULD REALLY BE CALLED :: test_J_Align.cpp ***/
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
#include "vectormath.h"

// LibIgl includes
#include <igl/readOFF.h>
#include <igl/writeOFF.h>
#include <igl/viewer/Viewer.h>
#include <igl/per_edge_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/print_vector.h>
#include <igl/signed_distance.h> 
#include <igl/boundary_loop.h>
#include <igl/slice.h>
#include <igl/random_dir.h>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>

// libIGL includes [ note :: can put in .h file instead ] - needed for addition boundary logic
#include <igl/boundary_loop.h>
#include <igl/per_face_normals.h>
#include <algorithm> // needed for std-lib [containers + iterators]
// viewer logic [ For later debugging ... crucial only for main method ]
#include <igl/viewer/Viewer.h>
#include <igl/cat.h>

// remember to set each rigid body instance's {c,theta} and {cvel,w} to 0 
	// ... nvm, {c,theta} zeroed out @ constructor time!
// #TODO :: for debug output, use print_vector? perhaps. Check this out later though!

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
const double rho = 1;
const double h = 0.01;

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
void solve_config_deltas(const Eigen::MatrixXd& vBoundary, const Eigen::MatrixXd bndExtForceMat, const RigidBodyInstance* body, Eigen::Vector3d& delta_cvel, Eigen::Vector3d& delta_omega);
void updateConfiguration(RigidBodyInstance* ScanBody, const Eigen::Vector3d& delta_cVel,const Eigen::Vector3d& delta_omega);
void solveTransformation(const RigidBodyInstance* ScanBody, Eigen::Matrix4d& T);

// added from normal boundary solving
int mod(int a, int b);
int detectFace(const Eigen::Vector2i& e_i, const Eigen::MatrixXi& F);
void solveExternalForcesForBoundary(const std::vector<int>& boundaryLoop, const Eigen::MatrixXd& remeshV, const Eigen::MatrixXi& remeshF, Eigen::MatrixXd& boundaryV, Eigen::MatrixXd& extBndForceMat);


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
			 *******************************/
			// [1] CALCULATE boundary vertices for both scans
			// alrlight, here are some of our first issues. We need to be working with the remshed surface only here! Let's rip code from later, up to here!

			// generate boundary vertices from interpolating surface
			Eigen::MatrixXd boundaryOne;
			Eigen::MatrixXd boundaryTwo;
			int bvOne_num;
			int bvTwo_num;

			INTERP_SURF::generateBoundaryVertices(scan1.V, scan1.F,boundaryOne); 
			INTERP_SURF::generateBoundaryVertices(scan2.V, scan2.F,boundaryTwo); 
			bvOne_num = boundaryOne.rows();
			bvTwo_num = boundaryTwo.rows();
			
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
			solveExternalForcesForBoundary(remeshBoundaryOne, flowed.V, flowed.F, boundaryV_1, extBndForceMat_1);  

			// [4] SOLVE for external forces assoc w/boundary 2 
			Eigen::MatrixXd boundaryV_2;
			Eigen::MatrixXd extBndForceMat_2;
			solveExternalForcesForBoundary(remeshBoundaryTwo, flowed.V, flowed.F, boundaryV_2, extBndForceMat_2);  

			// [5] solve for config_deltas
			solve_config_deltas(boundaryV_1, extBndForceMat_1, ScanOneBody, delta_cVel_1,delta_omega_1);
			solve_config_deltas(boundaryV_2, extBndForceMat_2, ScanTwoBody, delta_cVel_2,delta_omega_2);

			// [6] update rigid body instance configurations
			updateConfiguration(ScanOneBody,delta_cVel_1,delta_omega_1);
			updateConfiguration(ScanTwoBody,delta_cVel_2,delta_omega_2);

			// [7] GENERATE transformation [ rot + trans ] components
			Eigen::Matrix4d T1 = Eigen::Matrix4d::Identity();
			Eigen::Matrix4d T2 = Eigen::Matrix4d::Identity();

			solveTransformation(ScanOneBody, T1);
			solveTransformation(ScanTwoBody, T2);

			// [8] CREATE global mesh, contain the two inputs aligned based on impulses
			HELPER::applyRigidTransformation(scan1_rest.V,T1,scan1.V);
			HELPER::applyRigidTransformation(scan2_rest.V,T2,scan2.V);
			igl::cat(1,scan1.V,scan2.V,result.V);
			igl::cat(1,scan1.F, MatrixXi(scan2.F.array() + scan1.V.rows()), result.F);

			// concat the interpolating surface too [ debug purposes ]
			//igl::cat(1,result.V, flowed.V, debug_result.V);
			//igl::cat(1,result.F, MatrixXi(flowed.F.array() + result.V.rows()), debug_result.F);

			// concat normals of scan 1 [ debug purposes ]
			//igl::cat(1,result.V, normals.V, debug_result.V);
			//igl::cat(1,result.F, MatrixXi(normals.F.array() + result.V.rows()), debug_result.F);

			//igl::writeOFF(output_bad_mesh, interp.V, interp.F);
			//std::string output_mesh = "iterMesh[" + std::to_string(iters) + "].off";
			//igl::writeOFF(output_mesh, result.V,result.F);
			//igl::writeOFF(output_mesh, debug_result.V,debug_result.F);
			//viewer.data.set_mesh(debug_result.V,debug_result.F);
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


/*
	ASSUME such a call is executed :: 
	std::vector<int> remeshBoundaryOne;
	Eigen::MatrixXd boundaryV_1;
	Eigen::MatrixXd extBndForceMat_1;
	solveExternalForcesForBoundary(remeshBoundaryOne, remeshed.V, remeshed.F, boundaryV_1, extBndForceMat_1);  
*/	

// NOTE that the <body> is assoc with a <boundary> in thsi set up.
void solve_config_deltas(const Eigen::MatrixXd& vBoundary, const Eigen::MatrixXd bndExtForceMat, const RigidBodyInstance* body, Eigen::Vector3d& delta_cvel, Eigen::Vector3d& delta_omega)
{
	Eigen::MatrixXd bndExtImpMat = h * bndExtForceMat;
	int numBV = bndExtImpMat.rows();
	
	for(int k = 0; k < numBV; ++k)
	{
		const Eigen::Vector3d J_k_ext = bndExtImpMat.row(k);

		// ADD UP CONTRIBUTIONS TO DELTA_CVEL
		Eigen::Vector3d delta_cVel_k = J_k_ext;
		delta_cvel += delta_cVel_k;
		
		// ADD UP CONTRIBUTIONS TO DELTA_OMEGA
		const Eigen::Vector3d c = vBoundary.row(k); 
		const Eigen::Vector3d C = body->c; 
		const Eigen::Vector3d theta = body->theta;

		Eigen::Vector3d delta_omega_k = VectorMath::rotationMatrix(theta) * ((VectorMath::rotationMatrix(-1.0 * theta) * (C-c)).cross(VectorMath::rotationMatrix(-1.0 * theta) * J_k_ext));
		delta_omega += delta_omega_k;
		// keeping track of normals data and boundary vertices
		//normals_info.V.row(k) = c;
		//normals_info.N.row(k) = delta_omega_k * 100; // for visualization purposes ( * 100 )
	}
	//HELPER::constructNormalField(normals_info.V, normals_info.N, normals.V, normals.F);
}

// #TODO :: something weird is happenign here.
void updateConfiguration(RigidBodyInstance* ScanBody, const Eigen::Vector3d& delta_cVel,const Eigen::Vector3d& delta_omega)
{
	ScanBody->cvel += delta_cVel;

	ScanBody->w += delta_omega; 

	// UPDATE CONFIGURATION 
	// VANISH OUT ( LINEAR & ANGULAR ) VELOCITIES - akin to 'molasses'
	ScanBody->c += (ScanBody->cvel) * h; 
	const Eigen::Matrix3d rotMat = VectorMath::rotationMatrix(h * ScanBody->w) * VectorMath::rotationMatrix(ScanBody->theta);
    ScanBody->theta = VectorMath::axisAngle(rotMat);
	ScanBody->cvel.setZero();
	ScanBody->w.setZero();
}

// NEED axis-angle recovery formulas here!
// #TODO :: doing something wrong in the transformation here for rotatinoal component. Needs to be fixed!
// bring up with Dr. Vouga --- is this axis angle represnetation correct!
void solveTransformation(const RigidBodyInstance* ScanBody, Eigen::Matrix4d& T)
{
	const Eigen::Vector3d t_comp = ScanBody->c - ScanBody->c_0;
	//#TODO ::  wait a sec, should this vector be getting normalized? Could be an issue
	//const Eigen::Vector3d theta_diff = (ScanBody->theta - ScanBody->theta_0).normalized();
	const Eigen::Vector3d theta_diff = (ScanBody->theta - ScanBody->theta_0);
    // #TODO :: fix this. theta_diff = NaN, is definetely an error!
	const Eigen::Vector3d temp_diff = Eigen::Vector3d(0,0,0);
	// #TODO :: we expect <theta_diff> to be a more non-zero value. What is happening to ScanBody->theta here?
	const Eigen::Matrix3d r_comp = VectorMath::rotationMatrix(temp_diff);
	// I expect r_comp = I_3. if not, something else is wrong in the computational of these rotation matrices
	// #NOTE :: this is definetely getting correctly written :-)
	T.block<3,1>(0,3) = t_comp; 
	T.block<3,3>(0,0) = r_comp;
}

/*
 *  Just a basic modulus operation
 */
int mod(int a, int b)
{
	return (a%b+b)%b;
}

/*
 * GIVEN a specific edge, <e_i> and set of faces <F>
 * RETURN the unique face INDEX <f_i> corresponding to those two vertices indexed by the edge <e_i>
 * NOTE that <e_i> is a boundary edge, and thus always indexes into boundary vertices.
 * USE face_index due to <per_face_normals.h> library.
 */
int detectFace(const Eigen::Vector2i& e_i, const Eigen::MatrixXi& F)
{
	// SOLVE for the intersection of a 3d and a 2d vector.
	// JUST use a HM frequency map honestly. Iterate over <e_i> as if those were your keys instead.
	// COMPLEXITY :: [O(F*E) time, O(E) space]
	bool hasFace = false;
	int f_i = -1; 
	for(int i = 0; i < F.rows(); ++i)
	{
		// CONSTRUCT HM frequency of [face,edge] vertex indices.
		Eigen::Vector3i myFace = F.row(i);
		std::map<int,int> vFreq;
		vFreq[myFace(0)]++;
		vFreq[myFace(1)]++;
		vFreq[myFace(2)]++;
		vFreq[e_i(0)]++;
		vFreq[e_i(1)]++;

		// IF we hit two vertices, we have a matching face
		if(vFreq[e_i(0)] == 2 && vFreq[e_i(1)] == 2) 
		{
			f_i = i;
			break;
		}
	}
	return f_i;
}


/*
 * RETURNS a set of boundary vertices and assoc external force normals for each vertex
 * INPUT = [remeshBoundaryOne, remeshed.V, remeshed.F]
 * OUTPUT = [boundryV, extBndForceMat]
 */
void solveExternalForcesForBoundary(const std::vector<int>& boundaryLoop, const Eigen::MatrixXd& remeshV, const Eigen::MatrixXi& remeshF, Eigen::MatrixXd& boundaryV, Eigen::MatrixXd& extBndForceMat)
{
	// [1] Calculate all per face normals, and normalize them
	// #TODO :: <Z> := should I be concerned about cases of faces with degenerate normals?
	// #TODO :: why in the example is it that degenerate faces are assigned (1/3. 1/3, 1/3)^0.5? 
	Eigen::Vector3d Z = Eigen::Vector3d(1,1,1).normalized();
	Eigen::MatrixXd N_face;
	igl::per_face_normals(remeshV, remeshF, Z, N_face);

	// [2] construct set of edges from <boundary_loop_vertices>
	// Eigen::Matrix3i <vFIndices> is of form (v_i, v_j, f_i). It'll serve as a useful index map.
	// #TODO :: ASSERT all faces are unique [ do this later ]
	int numBV = boundaryLoop.size();
	Eigen::MatrixXi vfIndices;	 // size = |bndryLoop| 
	vfIndices.resize(numBV,3);	
	for ( int i = 0; i < numBV; ++i)
	{
		int v_i = boundaryLoop[mod(i,numBV)];
		int v_i_plus_1 = boundaryLoop[mod((i+1), numBV)];
		Eigen::Vector2i e_i = Eigen::Vector2i(v_i,v_i_plus_1);
		int f_i = detectFace(e_i, remeshF);
		if(f_i == -1)
			cout << "ERROR :: INVALID FACE INDEX\n";
		Eigen::Vector3i myData = Eigen::Vector3i(v_i,v_i_plus_1,f_i);
		vfIndices.row(i) = myData;
	}

	// [3] Construct normalized edge vectors and face normals.
	// TAKE their cross product :: create new matrix of f_i_ext. 
	// #TODO :: ordering of X-prod vectors is tied to how the boundary is iterated over too.
	Eigen::MatrixXd f_i_ext_Mat;
	Eigen::MatrixXd V_norm;

	f_i_ext_Mat.resize(numBV, 3);
	V_norm.resize(numBV, 3);

	for(int i = 0; i < numBV; ++i)
	{
		Eigen::Vector3i myData = vfIndices.row(i);
		int v_i = myData(0);
		int v_i_plus_1 = myData(1);
		int f_i = myData(2);

		// Solve for face normal
		Eigen::Vector3d f_norm;
		if(f_i == -1) 
		{ 
			// NOTE :: I do not expect this err to arise at all.
			std::cout << "f_i = [-1] ERROR" << std::endl;
			f_norm = Eigen::Vector3d(1,0,0);
			exit(0);
		}
		f_norm = N_face.row(f_i);	

		// solve for edge normal
		Eigen::Vector3d e_i = (remeshV.row(v_i_plus_1) - remeshV.row(v_i));
		double e_i_len = e_i.norm();
		e_i.normalize();

		// Solve for the external force for bndry vertex
		Eigen::Vector3d f_i_ext = (f_norm.cross(e_i)).normalized() * e_i_len; 
		f_i_ext_Mat.row(i) = f_i_ext;

		// Construct boundary vertex assoc w/specific external force <f_i_ext>  
		Eigen::Vector3d c = 0.5 * (remeshV.row(v_i) + remeshV.row(v_i_plus_1));
		V_norm.row(i) =  c;
	}

	// ASSERT :: cardinality(N_edges) matches cardinality(boundary_one_vertices)
	if(numBV != f_i_ext_Mat.rows())
	{
		std::cout << "SIZE MISMATCH ERROR" << std::endl;
		exit(0);
	}

	boundaryV = V_norm;
	extBndForceMat = f_i_ext_Mat;

}
