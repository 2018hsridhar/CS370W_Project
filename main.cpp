// *** NOTE :: the # bndry vertices differs in both scans. You must get the cardinality here, for proper index starts.
// my lib includes
#include "interpSurface.h"
#include "remesh.h"
#include "glob_defs.h"

// libIGL includes [ note :: can put in .h file instead ]
#include <igl/writeOFF.h>
#include <igl/writeOBJ.h>
#include <igl/readOFF.h>
#include <igl/boundary_loop.h>
#include <igl/slice.h>
// #TODO :: you might not use all of these <normal> types BTW. 
#include <igl/per_corner_normals.h>
#include <igl/per_edge_normals.h>
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>

// viewer logic ( for debugging later )
#include <igl/viewer/Viewer.h>

// namespace includes
using namespace Eigen;  
using namespace igl;
using namespace std;

#include <algorithm> // needed for std-lib [containers + iterators]

int mod(int a, int b);
std::vector<int> instersection(std::vector<int> &v1, std::vector<int> &v2);
int detectFace(const Eigen::Vector2i& e_i, const Eigen::MatrixXi& F);

int mod(int a, int b)
{
	return (a%b+b)%b;
}

/*
 * vector intersection code - ripped from <StackOverFlow>
 * https://stackoverflow.com/questions/19483663/vector-intersection-in-c
 * IDK why one must sort vectors for set intersections
 */
// oh crud ... those are eigen vectors. Could get intersection, BUT would take too much time! 
std::vector<int> instersection(std::vector<int> &v1, std::vector<int> &v2)
{
    std::vector<int> v3;

    std::sort(v1.begin(), v1.end());
    std::sort(v2.begin(), v2.end());

    std::set_intersection(v1.begin(),v1.end(),v2.begin(),v2.end(),back_inserter(v3));

    return v3;
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
	// O(F*E) time, O(E) space
	bool hasFace = false;
	int f_i = -1; 
	for(int i = 0; i < F.rows(); ++i)
	{
		Eigen::Vector3i myFace = F.row(i);

		// construct HM frequency of [face,edge] vertex indices.
		// #TODO :: refactor code to something simpler BTW
		std::map<int,int> vFreq;
		vFreq[myFace(0)]++;
		vFreq[myFace(1)]++;
		vFreq[myFace(2)]++;
		vFreq[e_i(0)]++;
		vFreq[e_i(1)]++;

		int numIncV = 0;
		if(vFreq[e_i(0)] == 2)
			numIncV++;
		if(vFreq[e_i(1)] == 2)
			numIncV++;
		if(numIncV == 2)
		{
			f_i = i;
			break;
		}
	}
	return f_i;
}

int main(int argc, char *argv[])
{
	
	// LOAD IN MESHES, COMMENCE TESTS
	// std::cout << "Executing (TEST) for offset surface generation." << std::endl;

	struct Mesh
	{
		Eigen::MatrixXd V; 
		Eigen::MatrixXi F;
	} scan1, scan2, interp, remeshed;

	Eigen::MatrixXd boundaryOne;
	Eigen::MatrixXd boundaryTwo;
	int bvOne_num;
	int bvTwo_num;

	if(!readOFF(GLOBAL::interpSurfGenScan1File,scan1.V,scan1.F)) {
		std::cout << "Failed to load partial scan one." << std::endl;
	} 

	if(!readOFF(GLOBAL::interpSurfGenScan2File,scan2.V,scan2.F)) {
		std::cout << "Failed to load partial scan two." << std::endl;
	}

	// generate interp surf, and extract boundry verts from both scans.
	INTERP_SURF::generateInterpolatingSurface(scan1.V,scan1.F,scan2.V,scan2.F, interp.V,interp.F);
	// generate remeshed surface
	double rel = REMESH::avgEdgeLenInputMeshes(scan1.V, scan1.F, scan2.V, scan2.F);
	rel /= 2.0;
	REMESH::remeshSurface(interp.V, interp.F,remeshed.V,remeshed.F, rel);

	// INCORPORATE THIS CODE BLOCK INTO YOUR ALGO PIPELINE NOW!
	INTERP_SURF::generateBoundaryVertices(scan1.V, scan1.F,boundaryOne); 
	INTERP_SURF::generateBoundaryVertices(scan2.V, scan2.F,boundaryTwo); 
	bvOne_num = boundaryOne.rows();
	bvTwo_num = boundaryTwo.rows();
	
	// get REMESHED mesh boundary vertices, for both scan parts
	/*
	Eigen::MatrixXd remeshBoundaryOne;
	Eigen::MatrixXd remeshBoundaryTwo;
	REMESH::getRemeshedBoundaryVerts(bvOne_num, bvTwo_num, remeshed.V, remeshed.F, remeshBoundaryOne, remeshBoundaryTwo);
	*/
	std::vector<int> remeshBoundaryOne;
	std::vector<int> remeshBoundaryTwo;
	REMESH::getRemeshedBoundaryVerts(bvOne_num, bvTwo_num, remeshed.V, remeshed.F, remeshBoundaryOne, remeshBoundaryTwo);
	cout << "Done collecting remeshed boundary vertices" << endl;

// GOAL :: just apply these algorithm for one of these scans. We'll get to the other boundry ... later!
// remember ... as long as F_ext is right and scaled, working on the remeshed vs interpolated boundary WILL NOT MATTER! scaling helps account for these issues.



	// [1] Calculate all per face normals, and normalize them
	// #TODO :: <Z> := should I be concerned about cases of faces with degenerate normals?
	// why is it that in the example, degenerate faces are assigned (1/3. 1/3, 1/3)^0.5 ? Weird
	Eigen::Vector3d Z = Eigen::Vector3d(1,1,1).normalized();
	Eigen::MatrixXd N_face;
	igl::per_face_normals(remeshed.V, remeshed.F, Z, N_face);

	// current issue ... desire a mapping :: edge -> face. Would be nice to have
	// can't I just collect the edges based on the boundary instead? huh? 

	// [2] construct set of edges from <boundary_loop_vertices>
	// BETTER IDEA - perform a test to discover a faceIdx, given two vertices. That'd be nice.
	// use both a VALID and INVALID edge here, for testing purposes. 
	// don't have to worry about [indexing/cardinality] 

	// #TODO :: get some debug help from Dr. Vouga here. Will need it later.
	// #TODO :: ASSERT that no face has an index of (-1).
	// #TODO :: run ASSERTION that all faces are unique, and not (-1) BTW. Is this expected?
	// Eigen::Matrix3i <vFIndices>, of form (v_i, v_j, f_i), will serve as a useful index map.
	Eigen::MatrixXi vfIndices;	
	vfIndices.resize(remeshBoundaryOne.size(),3);	
	for ( int i = 0; i < remeshBoundaryOne.size(); ++i)
	{
		int v_i = remeshBoundaryOne[mod(i,bvOne_num)];
		int v_i_plus_1 = remeshBoundaryOne[mod((i+1), bvOne_num)];
		Eigen::Vector2i e_i = Eigen::Vector2i(v_i,v_i_plus_1);
		// let us detect corresponding face here
		int f_i = detectFace(e_i, remeshed.F);
		if(f_i == -1)
			cout << "ERROR :: INVALID FACE INDEX\n";
		Eigen::Vector3i myData = Eigen::Vector3i(v_i,v_i_plus_1,f_i);
		vfIndices.row(i) = myData;
	}
	//cout << vfIndices << endl;

	// [3] Construct normalized edge vectors and face normals.
	// TAKE their cross product :: create new matrix of f_i_ext. 
	// #TODO :: be careful of ordering of vectors in X-product. Depends on way boundry is iterated over too.
	Eigen::MatrixXd f_i_ext_Mat;
	f_i_ext_Mat.resize(vfIndices.rows(), 3);
	Eigen::MatrixXd V_norm;
	V_norm.resize(vfIndices.rows(), 3);
	for(int i = 0; i < vfIndices.rows(); ++i)
	{
		Eigen::Vector3i myData = vfIndices.row(i);
		int v_i = myData(0);
		int v_i_plus_1 = myData(1);
		int f_i = myData(2);
		// #TODO :: fix this later
		if(f_i == -1) 
			continue;
		Eigen::Vector3d f_norm = N_face.row(f_i);	
		Eigen::Vector3d e_i = (remeshed.V.row(v_i_plus_1) - remeshed.V.row(v_i));
		double e_i_len = e_i.norm();
		e_i.normalize();
		Eigen::Vector3d f_i_ext = (f_norm.cross(e_i)).normalized() * e_i_len;
		f_i_ext_Mat.row(i) = f_i_ext;

		// extraneuous
		Eigen::Vector3d c = 0.5 * (remeshed.V.row(v_i) + remeshed.V.row(v_i_plus_1));
		V_norm.row(i) = c;
	}

	// you know what ... just refactor and include this code back in your old alignment code. That's easier! 
	// I don't want to rewrite this logic again!
}
