// CHECK if stitch surface has Nans, or [V,F] Nana,s [COMs, orientations], and if Forces are NaN
// which step is causing the blow up
// @ 17th iteration of CamelHeads ---> something is off. Fix this!
// so we get a stiffness matrix with NaN entries @ iteration 17 --- didn't Kazhdan mention something like this in his paper. Should corroborate his findings from earlier BTW!

// Vouga's questions
// [1]  are the triangels really flat or not
// [2]  print out the initial and final states of your test images. If it's close to zero, this is expected
// [3] Run on more tests cases. Will be easy to do at the moment :-)
// [4] gen stitched and remesh only every couple of stages. Must flow @ each stage though. You can't help that. This is a later extension after other issues get resolved 
// [5] what sort of entries are in the [MASS, STIFFNESS] matrices. We should check Nan's or invalid rows here, and ask why
// [6] do visualize your 3 surfaces of interest [ 2 boundaries, stitching surface ] too
// so for the Laplacian, it's given a decently valid remeshed mesh. aspect ratio doesn't seem to be an issue. And the call is basic:: it's <igl::cotMatrix(V,F,L)>, and computed once before we apply the MCF algorithm iteratively.
// Current theory :: your CGAL remeshing might be creating two accidental copies, instead of using same vertex. basically, one of the underlying operations here could be introducing FPE [ Floating-Point Exceptions ]. Thus, what you expect to be right, is not actually right.

// #TODO :: reduce time for compilation.
// #TODO :: crud, I have [a] refactoring and [b] automation work for today - Friday. Well, this is just fantastic! ... ehhh, what to do ...
// #TODO :: remove dependencies for some targets - {glfw, CoMISo, glew }. Just need dependencies for CS370W_PROJ_bin really.

// MY LIBRARIES  
#include "glob_defs.h"
#include "helpers.h"
#include "stitchedSurface.h"
#include "remesh.h"
#include "meanCurvatureFlow.h"
#include "sgd.h"  
#include "vectormath.h"
#include "rigidbodytemplate.h"
#include "rigidbodyinstance.h"
#include "j_align.h"

// LibIgl includes
#include <igl/writeOFF.h>
// viewer logic ( viewing mostly )
#include <igl/viewer/Viewer.h>
#include <igl/cat.h>
#include <igl/boundary_loop.h>

// DBEUGGING ... cotmatrix
#include <igl/cotmatrix.h>
#include <igl/doublearea.h>

// C++ includes
#include <iostream>
#include <fstream>
#include <algorithm> 

// C++ Namespace includes 
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
	Eigen::VectorXi bI; // boundary indices 
	Eigen::MatrixXd N;
	int numBV; 			// number of boundary vertices
} scan1, scan2, scan1_rest, scan2_rest, stitched, remeshed, flowed, result;

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
ofstream debugFile;

// METHOD BODY
// also, output the mesh into the remeshed mesh
// INPUT mesh can have some problems { your stitching algorithm )
	// ... crud, is this wrong again! OMG!
// GOOD ADVICE - write out each tests as seperate methods; Lest you get incredibly lost here.


int main(int argc, char *argv[])
{
	// EXECUTE IMPUSLE-BASED, J_ALIGN PIPELINE
	runPipeline();

/* igl::readOFF(TUTORIAL_SHARED_PATH "/stitchedSurfCase.off", scan1.V, scan1.F);
	Eigen::VectorXd dbla;
	igl::doublearea(scan1.V,scan1.F,dbla);
	cout << dbla << endl;
*/
	// note :: just put in a better script for <smallSurfaceArea>
	// output your boundary verts too!

	/* #TODO :: write a method to generate data. Put <runPipeline()> in its own file. It intrudes on what I desire in main.cpp.	
	 * or better yet ... seperate <data_generation> from <pipeline> code.
	 * ... hmmm ... a desired structure could be the following?
 	 *
	 * 
	 */

	// read the number of boundary-vertices for the camel
/*
	igl::readOFF(GLOBAL::pipelineScan1File,scan1.V,scan1.F);
	// NOTE :: get boundary loops from <stitched.F> instead
	// u need to output the sttiched mitch. Code is deterministic. Print to file here.
	std::vector<int> bndLoop
	igl::boundary_loop(scan1.F, bndLoop);
	cout << "# bndry verts = " << bndLoop.size() << endl;

			/
	// HERE, load in caes where the stiffness matrix blows up. Analyze why
	// also #TODO :: get that bool_opers code fixed later. This is still an issue.
*/

/*
 *** COTANGENT MATRIX CALCULATIONS 
	Eigen::MatrixXd V; 
	Eigen::MatrixXi F;
  	igl::readOBJ(TUTORIAL_SHARED_PATH "/surfaceBlowUpCase.obj", V, F);

	Eigen::SparseMatrix<double> L;
	igl::cotmatrix(V,F,L);

	cout << L.norm() << endl;
	double L_norm = L.norm();
	cout << L_norm << endl;
*/
	// uhhh, I don't think <L> at first is the issue, but rather, something later in the algorithm itself, that is probably provoking this issue. Need to take a more careful look at this now!

	// eval sum of diagonal coeffs too --- uhh, can't eval traces(matrices) for sparse(matrix)
	// take rowwise sum
	// wait a sec ... does it even make sense to take these norms? or not? HUH?
	// Eigen::VectorXd row_sum = V.rowwise().sum();
	// cout << row_sum << endl;
	// NOTE:: the sum of these row elements seems completely valid
	// NOTE :: we operate specifically on the Frobenius norm here.
	// #TODO :: omit rows containing [a] Nans or [b] Infs?
	// #TODO :: does L work off of facse, or vertices? I believe faces. 

	// nan check
	//cout << "Num vertices = " << V.rows() << endl;
	//cout << "Num faces = " << F.rows() << endl;

	// check if your vertices are on the boundary or NOT
	// we know vertex zero is on the boundary for sure ( based on indexing scheme )


	// laplacian is a [V,V] matrix. 
	// so we need to get the set of adj faces, and analyze those specific faces in more depth!
	// hmm ... these are cotangents ... when do they approach [-inf,inf]? 
	// cot(x) -> \inf as x -> 0
	// cot(x) -> -\inf as x -> \pi
	// so the aspect ratio needs to be checked, as well sa triangle flatness too, I guess, in the MCF algorithm.


	// let us also run a test over each area of the triangle, and see what areas we get for said triagnles. If we get something close to zero, that'll be an issue for us
	// hmm, some of them are small, but not small enough to make me suspect that they would blow up
/*
	cout << "Will print out triangle areas" << endl;
	Eigen::VectorXd dbla;
	igl::doublearea(V,F,dbla);
	// note :: just seeing how small triangle areas can get
	// there's only one really tiny triangle. 
	// but it does correspond to vertex indices that are failing [ e.g. 0, 100, 200 ]
	for(int i = 0; i < dbla.size(); ++i)
	{
		if(dbla(i) < 0.000001)
		{
			cout << dbla(i) << endl;
			cout << i << endl;
			Eigen::Vector3i myFace = F.row(i);
			cout << "[" << myFace(0) << ", " << myFace(1) << ", " << myFace(2) << "]\n";
		}
	}
	cout << "Done printing out triangle areas" << endl;
	exit(0);
*/
/*
	cout << "will exec nan check" << endl;
	// wait a sec... this isn't the right set of dims
	cout << "num rows in L " << L.rows() << endl;
	cout << "num cols in L " << L.cols() << endl; // 2409
	for(int i = 0; i < L.rows(); ++i)
	{
		for(int j = 0; j < L.cols(); ++j)
		{
			// #NOTE :: there's a faster way. [coeffRef] is quite slow.
			// wait a sec ... the first entry itself is off? what??
			if(isinf(L.coeffRef(i,j))) 
			{
				cout << "[i,j] = [" << i << ", " << j << "]\n";
			}
		}
	}
*/

	// indices corresponding to inf :: [0,100], [0,200], [100,0], [100,200], [200,0], [200,100], [200,200]
	// indices correspoding to nan :: [0,0], [100,100]. Two of these diagonal entries
	// let's print out those rows. let's try the first row

	// okay, there is a [-inf] and [inf] here too.
	// I preusme [-inf] + [inf] = [nan], right? 
	// where does this stem from though?
	// let us check those too!
/*
	cout << "------------------------------------------------------" << endl;
	for(int j = 0; j < L.cols(); ++j)
	{
		// wait a sec ... the first entry itself is off? what??
		cout << L.coeffRef(0,j) << endl;
	}
	cout << "------------------------------------------------------" << endl;
	cout << endl;

	// in lapaclian, equals -(sum) of all other row entries

	cout << "done with nan check" << endl;
*/

/*
	cout << L << endl;
	cout << "V_vertics" << endl;
	cout << V.norm() << endl;
	cout << "------------------------------------------------------" << endl;
	Eigen::MatrixXd flowed_V;
	MCF::computeMeanCurvatureFlow(V,F,0.001, flowed_V);
	cout << "------------------------------------------------------" << endl;
	cout << "flowed_V_vertices" << endl;
	cout << flowed_V.norm() << endl;
	cout << "------------------------------------------------------" << endl;
*/

	



	// so this mesh seems completely valid
	// are any of these triangles flat?
	// any thing to note about computation of stiffness matrices too?
	// okay, this code itself can't be changed.
	// which of these rows is Nan? let's take a sum of each row
	// eval sum, see triangle

/*
	igl::viewer::Viewer viewer;
	viewer.data.clear();
	viewer.data.set_mesh(V,F);
	return viewer.launch();
*/

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

			// GENERATE stitchedOLATING SURFACE
			//cout << "Generating stitchedolating surface" << endl;
			// let's output STITCHED_SURF for camel heads case, to same file, all the time
			// #TODO :: change <stitched> to <STITCH>, as its a better naming convention.
			STITCHED_SURF::generateStitchedSurface(scan1.V,scan1.F,scan2.V,scan2.F, stitched.V,stitched.F);
			igl::writeOFF(GLOBAL::stitching_output, stitched.V, stitched.F);

			double rEL = REMESH::avgEdgeLenInputMeshes(scan1.V,scan1.F,scan2.V,scan2.F);
			bool remSucc = REMESH::remeshSurface(stitched.V,stitched.F,remeshed.V,remeshed.F, rEL);
			igl::writeOFF(GLOBAL::remesh_after_succ_file, remeshed.V, remeshed.F);

			if(!remSucc)
			{
				std::cout << "BAD REMESHING!" << std::endl;
				std::cout << "Printing out stitchedolating surface" << std::endl;
				std::string output_bad_mesh = "badstitchedSurf[" + std::to_string(iters) + "].off";
				std::string scan1_bad = "badScan1.off";
				std::string scan2_bad = "badScan2.off";
				igl::writeOFF(output_bad_mesh, stitched.V, stitched.F);
				exit(0);
			}


			// APPLY MCF TO SURFACE [ denoted as 'flowed' here ]
			assert(MCF::meshHasBoundary(remeshed.V,remeshed.F));
			// #TODO :: fix MCF algorithm. failure is here. Nan's are getting introduced here. 
			// Need we re-read Kazhdan's ( surface pinch cases? Perhaps. Or what val is set too also ) 
			ofstream myFile (GLOBAL::mcfDebugText, std::ios_base::app); // need append mode
			if (myFile.is_open())
			{
				myFile << "Iteration [" << std::to_string(iters) << "]\n";
			}
			myFile.close();
			MCF::computeMeanCurvatureFlow(remeshed.V,remeshed.F,0.001, flowed.V);
			flowed.F = remeshed.F;

			// write mcf file ... see if Nan's
/*
			std::string remesh_mesh = "remeshMesh[" + std::to_string(iters) + "].off";
			igl::writeOFF(remesh_mesh, remeshed.V, remeshed.F);
			std::string mcf_mesh = "mcfMesh[" + std::to_string(iters) + "].off";
			igl::writeOFF(mcf_mesh, flowed.V, flowed.F);
*/
			// CHECK IF impulses improved alignment
			// #TODO :: why do we get a <NaN> value here on 87th iteration of camel heads? That's somewhat off.
			// IF we get a area to have NaN ... one of these vertices/triangles posses a NaN ... where is the problem in the pipeline. 
			// PUT some ASSERTIONS here!
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
			// #TODO :: change this threshold.  
			// simplest thing ... large # of iterations ( i.e. 500 ), don't stop early!
			// use derivative f of energy ( close to 0 ). Magnitude of grad(E)
			// compare both [trans,rot] of rigid bodies

			//if(local_energy < 0.1) // #TODO :: find a better approach here??
			if(local_energy < 0.025) // #TODO :: find a better approach here??
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
			// generate boundary vertices from stitchedolating surface
			// #TODO :: doesn't <generateBoundaryVertices> return a std::vector<int> in actuality? use that instead!
			// #TODO :: not the insanely limited scope of <boundaryOne|boundaryTwo> here.

			/*
			Eigen::MatrixXd boundaryOne;
			Eigen::MatrixXd boundaryTwo;

			STITCHED_SURF::generateBoundaryVertices(scan1.V, scan1.F,boundaryOne); 
			STITCHED_SURF::generateBoundaryVertices(scan2.V, scan2.F,boundaryTwo); 
			int bvOne_num = boundaryOne.rows();
			int bvTwo_num = boundaryTwo.rows();
			*/

			// #TODO :: write a better test. We want to compare boundary loops, not pointers. That test might've been bad.
			// I note that meshes passed here are [scan1.F, scan2.F], versus [stitched.F]
			// perhaps this is causing the difference for us? Let us check later!
/*
			std::vector<int> stitchBoundaryOne;
			std::vector<int> stitchBoundaryTwo;
			igl::boundary_loop(scan1.F,stitchBoundaryOne);
			igl::boundary_loop(scan2.F,stitchBoundaryTwo);
			int bvOne_num = stitchBoundaryOne.size();
			int bvTwo_num = stitchBoundaryTwo.size();
*/


			// NOTE :: get boundary loops from <stitched.F> instead
			std::vector<std::vector<int>> stitch_bndryLoop;
			igl::boundary_loop(stitched.F, stitch_bndryLoop);

			// GET both boundary loop vectors
			std::vector<int> stitch_firstBL = stitch_bndryLoop[0];
			std::vector<int> stitch_secondBL = stitch_bndryLoop[1];
			int bvOne_num = stitch_firstBL.size();
			int bvTwo_num = stitch_secondBL.size();

			// get FLOWED mesh ( after remeshing too, btw ) boundary vertices, for both scan parts
			std::vector<int> remeshBoundaryOne;
			std::vector<int> remeshBoundaryTwo;

			REMESH::getRemeshedBoundaryVerts(bvOne_num, bvTwo_num, flowed.V, flowed.F, remeshBoundaryOne, remeshBoundaryTwo);
			cout << "Done collecting remeshed boundary vertices" << endl;

			// #TODO :: eyeball both sets here
			//std::sort(stitchBoundaryOne.begin(), stitchBoundaryOne.end());
			//std::sort(stitchBoundaryTwo.begin(), stitchBoundaryTwo.end());
			
			// now sort and print them out ( let's focus on 1st boundary ) 
/*
			std::sort(stitch_firstBL.begin(), stitch_firstBL.end());
			std::sort(stitch_secondBL.begin(), stitch_secondBL.end());
			std::sort(remeshBoundaryOne.begin(), remeshBoundaryOne.end());
			std::sort(remeshBoundaryTwo.begin(), remeshBoundaryTwo.end());


			for ( int i : stitch_firstBL )
				std::cout << i << " ";
			std::cout << std::endl;
			std::cout << "----------------------------------------------------------------------" << std::endl;
			for ( int i : remeshBoundaryOne )
				std::cout << i << " ";
			std::cout << std::endl;
			std::cout << "----------------------------------------------------------------------" << std::endl;

			std::cout << "----------------------------------------------------------------------" << std::endl;

			for ( int i : stitch_secondBL)
				std::cout << i << " ";
			std::cout << std::endl;
			std::cout << "----------------------------------------------------------------------" << std::endl;
			for ( int i : remeshBoundaryTwo )
				std::cout << i << " ";
			std::cout << std::endl;
			std::cout << "----------------------------------------------------------------------" << std::endl;

			exit(0);
*/
			// #TODO :: move method below to another file? Possibly an issue?

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
