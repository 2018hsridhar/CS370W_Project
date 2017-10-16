// one question to raise about this alignment scheme, is this at least
// is alignment happening to the interpolating surface, or the flowed surface? It doesn't seem like it's actually happening to the flowed surface, but rather the interpolating surface. This raises a couple of concerns at the moment.
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
#include "vectormath.h"

// LibIgl includes
#include <igl/writeOFF.h>
#include <igl/viewer/Viewer.h>
#include <igl/per_edge_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/print_vector.h>
#include <igl/signed_distance.h> 
#include <igl/boundary_loop.h>
#include <igl/slice.h>
#include <igl/random_dir.h>

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
} scan1,scan2,scan1_rest,scan2_rest,interp,remeshed, flowed, debug_result, result, normals, normals_info;
// #TODO :: note that these <normals> are for <scan1> primarily, and <normals_info> is used to compute said <normals> 3D triangle-mesh structure.

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
void solve_config_deltas(Eigen::Vector3d& delta_cvel, Eigen::Vector3d& delta_omega ,struct Mesh m, RigidBodyInstance* body);
void updateConfiguration(RigidBodyInstance* ScanBody, const Eigen::Vector3d& delta_cVel,const Eigen::Vector3d& delta_omega);
void solveTransformation(const RigidBodyInstance* ScanBody, Eigen::Matrix4d& T);
void constructNormalField(Eigen::MatrixXd& V, Eigen::MatrixXd& N,
						Eigen::MatrixXd& V_res, Eigen::MatrixXi& N_edges);
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
			INTERP_SURF::generateOffsetSurface(scan1.V,scan1.F,scan2.V,scan2.F, interp.V,interp.F);
			double rEL = REMESH::avgEdgeLenInputMeshes(scan1.V,scan1.F,scan2.V,scan2.F);
			//rEL = rEL / 10.0; // #NOTE :: can actually visualize pinching occuring here
			rEL = rEL / 5.0; // #NOTE :: can actually visualize pinching occuring here
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

			// APPLY MCF TO SURFACE [ denoted as 'flowed' here ]
//			cout << "Flowing the interpolating surface" << endl;
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
			 *** MAKE SURE SURFACE IS PINCHED IN ***
			 *** J_ext as a vec field, ... in pinched surface *** 
			 *******************************/
			// [1] Compute vertex normals for both scans
			// #TODO :: are these normals updating as expected, and do they correspond to scans under the transformation? I believe so, due to the set of conditions stipulated later here.
			igl::PerVertexNormalsWeightingType weighting = PER_VERTEX_NORMALS_WEIGHTING_TYPE_UNIFORM; 
			igl::per_vertex_normals(scan1.V,scan1.F,weighting,scan1.N);
			igl::per_vertex_normals(scan2.V,scan2.F,weighting,scan2.N);

			// [2] Calculate boundary vertices for both scans
			igl::boundary_loop(scan1.F,scan1.bI);
			igl::slice(scan1.V,scan1.bI,1,scan1.bV);
			scan1.numBV = scan1.bV.rows();
				
			igl::boundary_loop(scan2.F,scan2.bI);
			igl::slice(scan2.V,scan2.bI,1,scan2.bV);
			scan2.numBV = scan2.bV.rows();

			// we can just assume \rho Vol = 1 ( unit mass, really )
			// [3] Solve for delta_cVel, delta_omega, for both scans
			Eigen::Vector3d delta_cVel_1 = Eigen::Vector3d(0,0,0);
			Eigen::Vector3d delta_omega_1 = Eigen::Vector3d(0,0,0);

			Eigen::Vector3d delta_cVel_2 = Eigen::Vector3d(0,0,0);
			Eigen::Vector3d delta_omega_2 = Eigen::Vector3d(0,0,0);

			solve_config_deltas(delta_cVel_1,delta_omega_1,scan1, ScanOneBody);
			solve_config_deltas(delta_cVel_2,delta_omega_2,scan2, ScanTwoBody);

			updateConfiguration(ScanOneBody,delta_cVel_1,delta_omega_1);
			updateConfiguration(ScanTwoBody,delta_cVel_2,delta_omega_2);

			// [4] GENERATE transformation ( rotation + translation ) components
			// for both partial scan boundaries
			Eigen::Matrix4d T1 = Eigen::Matrix4d::Identity();
			Eigen::Matrix4d T2 = Eigen::Matrix4d::Identity();

			solveTransformation(ScanOneBody, T1);
			solveTransformation(ScanTwoBody, T2);

            //std::cout << T1 << std::endl;
            //std::cout << T2 << std::endl;

			// [5] CREATE global mesh, contain the two inputs aligned based on impulses
			HELPER::applyRigidTransformation(scan1_rest.V,T1,scan1.V);
			HELPER::applyRigidTransformation(scan2_rest.V,T2,scan2.V);
			igl::cat(1,scan1.V,scan2.V,result.V);
			igl::cat(1,scan1.F, MatrixXi(scan2.F.array() + scan1.V.rows()), result.F);

			// concat the interpolating surface too, in this example [ going to visualize this too ]
			//igl::cat(1,result.V, flowed.V, debug_result.V);
			//igl::cat(1,result.F, MatrixXi(flowed.F.array() + result.V.rows()), debug_result.F);

			// concat the normals of scan 1 too [ too visualize ] 
			igl::cat(1,result.V, normals.V, debug_result.V);
			igl::cat(1,result.F, MatrixXi(normals.F.array() + result.V.rows()), debug_result.F);

			//igl::writeOFF(output_bad_mesh, interp.V, interp.F);
			//std::string output_mesh = "iterMesh[" + std::to_string(iters) + "].off";
			//igl::writeOFF(output_mesh, result.V,result.F);
			//igl::writeOFF(output_mesh, debug_result.V,debug_result.F);
			viewer.data.clear();
			//viewer.data.set_mesh(result.V,result.F);
			viewer.data.set_mesh(debug_result.V,debug_result.F);

			++iters;
			break;
		}
		default:
			return false;	
	}
	return true;
}

void solve_config_deltas(Eigen::Vector3d& delta_cvel, Eigen::Vector3d& delta_omega ,struct Mesh m, RigidBodyInstance* body)
{
	// ITERATE over boundary verts; calculate interpolating bnd {norm, vert}


	/*
     **************************************************************************
	 * You'll have numBV normals and numBV new vertices. Construct from here!
	 * use your global vars ( def above ) to help debug this
     **************************************************************************
     */
	normals_info.V = Eigen::MatrixXd(m.numBV,3);
	normals_info.N = Eigen::MatrixXd(m.numBV,3);	
	for(int k = 0; k < m.numBV; ++k) 
	{
		// get idx of 2 bndry verts of interest	
		int i = k % m.numBV;
		int j = (k + 1) % m.numBV;

		Eigen::Vector3d v_i = m.bV.row(i);
		Eigen::Vector3d v_j = m.bV.row(j);

		Eigen::Vector3d n_i = m.N.row(i);
		Eigen::Vector3d n_j = m.N.row(j);

		// eval NEGATED average two normals
		Eigen::Vector3d n_avg = -0.5 * (n_i + n_j);
	
		// Rescale F_k_ext, by edge length, and solve for J_k_ext
		n_avg.normalize();
		double edge_len = (v_j - v_i).norm();
		//Eigen::Vector3d n_avg_rescaled = edge_len * n_avg;
		// #TODO :: is this too slow or not?? unsure
		Eigen::Vector3d n_avg_rescaled = edge_len * n_avg;

		Eigen::Vector3d F_k_ext = n_avg_rescaled; 
		Eigen::Vector3d J_k_ext = h * F_k_ext;

		// ADD UP CONTRIBUTIONS TO DELTA_CVEL
		Eigen::Vector3d delta_cVel_k = J_k_ext;
		delta_cvel += delta_cVel_k;

		// ADD UP CONTRIBUTIONS TO DELTA_OMEGA
/*
		const Eigen::Vector3d c = 0.5 * (v_i + v_j); 
		const Eigen::Vector3d C = body->c; 			// ... isnt this from template?
		const Eigen::Vector3d theta = body->theta;
		Eigen::Vector3d delta_omega_k = VectorMath::rotationMatrix(theta) * ((VectorMath::rotationMatrix(-theta) * (C-c)).cross(VectorMath::rotationMatrix(-theta) * J_k_ext));
		delta_omega += delta_omega_k;
*/
		const Eigen::Vector3d c = 0.5 * (v_i + v_j); 
		const Eigen::Vector3d C = body->c; 			// ... isnt this from template?
		const Eigen::Vector3d theta = body->theta;
		// #TODO :: ask if this calculation of <delta_omega_k> is correct or NOT!
		Eigen::Vector3d delta_omega_k = VectorMath::rotationMatrix(theta) * ((VectorMath::rotationMatrix(-1.0 * theta) * (C-c)).cross(VectorMath::rotationMatrix(-1.0 * theta) * J_k_ext));
		delta_omega += delta_omega_k;

		// keeping track of normals data and boundary vertices
		normals_info.V.row(k) = c;
		//normals_info.N.row(k) = n_avg_rescaled; // according to edge len ( also F_k_ext too)
		// so F_k_ext looks right. now for delta_omega k's
		normals_info.N.row(k) = delta_omega_k * 100; // for visualization purposes ( * 100 )
	}
	constructNormalField(normals_info.V, normals_info.N, normals.V, normals.F);
	// we expect delta_omega = 0 , in the 2 planes case. if not, there is a bug somewhere!
	// this is not high priority anymore!
	std::cout << "delta_omega vector = \n";
	std::cout << delta_omega << std::endl;
	std::cout << "delta_omega magnitude = \n";
	std::cout << delta_omega.norm() << std::endl;
}

// #TODO :: something weird is happenign here.
void updateConfiguration(RigidBodyInstance* ScanBody, const Eigen::Vector3d& delta_cVel,const Eigen::Vector3d& delta_omega)
{
	ScanBody->cvel += delta_cVel;

	std::cout << "Before, value of w = " << ScanBody->w << std::endl;
	ScanBody->w += delta_omega; // #TODO :: fix this calculation!
	std::cout << "After, value of w = " << ScanBody->w << std::endl;

	// UPDATE CONFIGURATION 
	// VANISH OUT ( LINEAR & ANGULAR ) VELOCITIES - akin to 'molasses'
	ScanBody->c += (ScanBody->cvel) * h; 
	const Eigen::Matrix3d rotMat = VectorMath::rotationMatrix(h * ScanBody->w) * VectorMath::rotationMatrix(ScanBody->theta);
	// #TODO :: check that this update logic makes sense! huge difference between "=" and "+=" here!
    std::cout << "BEFORE, theta = " << ScanBody->theta << std::endl;
    ScanBody->theta = VectorMath::axisAngle(rotMat);
    std::cout << "AFTEr, theta = " << ScanBody->theta << std::endl;
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
	std::cout << "Comparing two axis-angle vectors" << std::endl;
	// #TODO :: we expect <theta_diff> to be a more non-zero value. What is happening to ScanBody->theta here?
    std::cout << theta_diff << std::endl;
    std::cout << temp_diff << std::endl;
	const Eigen::Matrix3d r_comp = VectorMath::rotationMatrix(temp_diff);
	// I expect r_comp = I_3. if not, something else is wrong in the computational of these rotation matrices
	// #NOTE :: this is definetely getting correctly written :-)
	T.block<3,1>(0,3) = t_comp; 
	T.block<3,3>(0,0) = r_comp;
}

void constructNormalField(Eigen::MatrixXd& V, Eigen::MatrixXd& N,
						Eigen::MatrixXd& V_res, Eigen::MatrixXi& N_edges)
{
	// #TODO :: convert this to it's own method body


	// #NOTE :: you'll convert this into a method later that takes in (V,N) and returns (V',E) to represent these visualized normals.
	// apply the normal scaling operation here
	Eigen::MatrixXd offsetedV = V + N;
	//scan2.V = (scan1.V).rowwise() + n.transpose(); 

	// concat the two vertex meshes too ( need this for indexing to properly work )
	igl::cat(1,V,offsetedV,V_res);

	// set up edges between [Scan1.V,scan2.V]
	std::cout << "COMMENCING set up of edges corresponding to normals\n";
	int offset = V.rows();
	N_edges = Eigen::MatrixXi(offset,3);
	for(int i = 0; i < offset; ++i)
	{
		Eigen::VectorXi e_i = Eigen::Vector3i(i,i + offset, i); // yeah, this is a problem #RIP
		// hmm ... cheat, set one of these the same ????
		N_edges.row(i) = e_i;	
	}	
	std::cout << "DONE setting up edges corresponding to normals\n";
}
