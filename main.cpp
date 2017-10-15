// idea - to visualize the field ; take each vertex v_i, it's normal n_i, and append the magnitude of the normal to another vertex v_j = v_i + n_i. generate your mesh from that! Assoc each edge properly too
// in essense, we generate this interesting pseudo-triangle mesh
// we can try out a basic one with 6 vertices ... why not?
// we can then place this back atop and see if it matches up with our intuition!
// unfortunately, LibIgl just does not provide this capability. I wish it did!
// RENAME THIS FUNCTIONALITY ::: impulse-based normals visualization 

// all this includes might not be necessary BTW!

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
// I have a good feeling that it is the vertex normals, versus the edge normals, that are needed in this situation! 
#include <igl/per_edge_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/print_vector.h>
#include <igl/signed_distance.h> 
#include <igl/boundary_loop.h>
#include <igl/slice.h>

// needed for randomized of "(N)ormals" here
#include <igl/random_dir.h>

// remember to set each rigid body instance's {c,theta} and {cvel,w} to 0 
	// ... nvm, {c,theta} zeroed out @ constructor time!
// #TODO :: for debug output, use print_vector
// Namespace includes 
using namespace Eigen;  
using namespace std;
using namespace igl;

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

// let us now test this with a random perturbation of vectors instead!
// #TODO :: see if functionality can be added to increase edge thickness or edit their color.
//const Vector3d n(0,0,1); // unit normal, in z-dir ( in other dirs, we can't see anything )
// better idea ... let's just read the planeXy image. we can tell from that itself, and construct edges from there!!

void constructNormalField(Eigen::MatrixXd& V, Eigen::MatrixXd& N, Eigen::MatrixXd& V_res, Eigen::MatrixXi& N_edges);

int main(int argc, char *argv[])
{
	// #TODO :: work on setting up the visualization correctly!
	igl::viewer::Viewer viewer;
	std::cout << "Executing creation of normals mesh here\n";

	// this is defined to be "planeXy.off" btw
	if(!readOFF(GLOBAL::pipelineScan1File,scan1.V,scan1.F)) {
		std::cout << "Failed to load partial scan one." << std::endl;
		return -1;
	} 
	readOFF(GLOBAL::pipelineScan1File,scan1.V,scan1.F);

	// we'll construct a normal field here, for testing purposes
	// make a set of randomized vectors corresponding to scan1.V.rows() in place of N
	Eigen::MatrixXd N(scan1.V.rows(), 3);
	for(int i = 0; i < N.rows(); ++i)
		N.row(i) = igl::random_dir();

	constructNormalField(scan1.V,N, result.V,result.F);

	viewer.data.clear();
	viewer.data.set_mesh(result.V, result.F);
	return viewer.launch();


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
