// ERR - in some cases, these bool opers seg fault? Perhaps invalid shapes being used under some of these settings? I need to check this!
// #TODO :: you're going to have to extend this codebase to namespace set up, becuase you'll be using KINNECT data too!
// #TODO :: REMEMBER that your inputs are 3D trimeshs's with CLEARLY DEFINED BOUNDARIES! you're not passing those 2d-mirrors here! That was an err in logic.

// YES - BOOL OPERS DO NOT WORK AS EXPECTED. you can't take the intersection of 2 different objs, even if they share the same plane, to get a boundary. In addition, all bool operations will result in water tight meshes, and the output of vertices and faces is 
// [1] not guaranteed to be in order as your inputs
// [2] will most likely possess floating point exception [ FPE ] errors

// #TODO :: we'd like to use a plane, but we lack representation of infinities in finite, discrete mesh structures. So we'll just "pretend" to have one by using a VERY huge object that possesses.
// [1] huge z-coords scaling ( try this first )
// [2] huge y-coords scaling 
// #TODO :: another next step, is to be able to generate some random perturbations of them. Add functionality for that  [ use <random_dir.h> to help with this. Translation + rotational components. ]. You're basically just projecting both boundaries somewhere into space.

// USER_DEFINED LIBRARIES
#include "bool_opers.h" // primary include
#include "glob_defs.h"
#include "interpSurface.h"
#include "remesh.h"
#include "vectormath.h" // NOTE :: needed for projecting a mesh randomly into space
#include "tutorial_shared_path.h"
#include "helpers.h"

// LIBIGL INCLUDES
// specific to CSG/bool opers + unit-sphere scaling
#include <igl/readOFF.h>
#include <igl/writeOFF.h>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/viewer/Viewer.h>

// includes for CSG/bool opers AND unit-sphere scaling
#include <igl/copyleft/cgal/mesh_boolean.h>
#include <igl/centroid.h> 
#include <igl/random_dir.h>
// EIGEN INCLUDES
#include <Eigen/Core>

// C++ STATIC LIB INCLUDES
#include <iostream>
#include <iostream>
#include <fstream>

// C++ NAMESPACES 
using namespace Eigen;  
using namespace std;
using namespace igl;

namespace BOOL_OPERS 
{
	void testMeshThrow()
	{
	/*
		igl::readOBJ(GLOBAL::boolOpersMeshOne,obj1_preSc.V,obj1_preSc.F); // <xxx>
		igl::readOBJ(GLOBAL::boolOpersMeshTwo,obj2_preSc.V,obj2_preSc.F); // cube
		applyUnitSphereRescaling(obj1_preSc.V, obj1_preSc.F, obj1.V,obj1.F);
		obj2.V = obj2_preSc.V;
		obj2.F = obj2_preSc.F;

		// throw both meshes ( obj1, obj2)
		throwMeshInSpace(obj1.V,obj1.F, frag1.V, frag1.F);
		throwMeshInSpace(obj2.V,obj2.F, frag2.V, frag2.F);

		igl::writeOFF(GLOBAL::boolOpersFragOne,frag1.V,frag1.F);
		igl::writeOFF(GLOBAL::boolOpersFragTwo,frag2.V,frag2.F);
	*/
	}

	void throwMeshInSpace(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::MatrixXd& V_p, Eigen::MatrixXi& F_p)
	{
		// apply rotation
		Eigen::Vector3d axisAngle = igl::random_dir();
		Eigen::Matrix3d rotationMat = VectorMath::rotationMatrix(axisAngle);
		Eigen::MatrixXd perturbedByR = V * rotationMat;

		// apply translation
		double scale = 10.0; // #NOTE :: can set to desired value.
		Eigen::Vector3d randomDir = igl::random_dir();
		Eigen::Vector3d translationVec = scale * randomDir;
		Eigen::MatrixXd perturbedByT = perturbedByR.rowwise() + translationVec.transpose(); 

		// output results
		V_p = perturbedByT;
		F_p = F;
	}

	/*
	* To note about boolean operations 
	* YES, meshes must be water tight for bool operations to be applicable. 
	*    they must also be orientable.
	* Hence why they fail in cases such as <circle.off> but not <camelHead.off>
	* NOTE :: V1 = left, V2 = right #TODO update naming conventions.
	* We'll always use the unit-cube ( it's a basic cutter for us ). We'll vary the other mesh though.
	* #TODO :: randomized these offset vectors, but always to be unit vectors ( max dist = 1 ). 
	*/
	void generateBoolOperMeshes(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::MatrixXd& V1, Eigen::MatrixXi& F1, Eigen::MatrixXd& V2, Eigen::MatrixXi& F2)
	{
		// APPLY offsets to V, both [-|+] x-axis. Generate two additional meshes for cutting
		Eigen::Vector3d offset_left = Eigen::Vector3d(-0.5,0,0);
		Eigen::Vector3d offset_right = Eigen::Vector3d(0.5,0,0);

		V1 = V.rowwise() + offset_left.transpose();
		V2 = V.rowwise() + offset_right.transpose();
		F1 = F2 = F;
	}


	// #TODO:: include as a test method, run in a testing system ( 
	// Method responsibility - asserts preservation of order of vertex indexing between the interpolating and remeshed surfaces. 
	// equivalently, assserts that the n vertices (interp) = first n ( remesh )
	// ... ( e.g. akin to bb test, from Amazon ) 
	void runIndexTestInterpRemesh()
	{
	/*
		if(!readOFF(GLOBAL::pipelineScan1File,scan1.V,scan1.F)) {
			std::cout << "Failed to load partial scan one." << std::endl;
		} 

		if(!readOFF(GLOBAL::pipelineScan2File,scan2.V,scan2.F)) {
			std::cout << "Failed to load partial scan two." << std::endl;
		}

		std::cout << "Executing pipeline for the following meshes" << std::endl;
		std::cout << "Mesh one [" << GLOBAL::pipelineScan1File << "]" << std::endl;
		std::cout << "Mesh two [" << GLOBAL::pipelineScan2File << "]" << std::endl;

		INTERP_SURF::generateInterpolatingSurface(scan1.V,scan1.F,scan2.V,scan2.F, interp.V,interp.F);
		double rEL = REMESH::avgEdgeLenInputMeshes(scan1.V,scan1.F,scan2.V,scan2.F);
		rEL = rEL / 5.0; 
		bool remSucc = REMESH::remeshSurface(interp.V,interp.F,remeshed.V,remeshed.F, rEL);
		if(!remSucc)
		{
			std::cout << "BAD REMESHING!" << std::endl;
			exit(0);
		}

		igl::writeOFF("interpolating_surface.off", interp.V,interp.F);
		igl::writeOFF("remeshed_surface.off", remeshed.V,remeshed.F);
	*/
	}


	/* 
	* Idea taken from calculating smallest bounding sphere for a mesh
	* and scaling according to its radius
	* should be noted that for meshes contained in the unit sphere, it'll still work
	* just based on a fractional radius
	* Note that operations are limited to watertight meshes
	*/
	void applyUnitSphereRescaling(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::MatrixXd& V_scaled, Eigen::MatrixXi& F_scaled)
	{
		/*
	 	 * [1] calcualte COM ( = centroid ) 
		 * We use the centedoid here since we can assume our meshes are of constant density. Note 
		 *    that centroid calculatison work only for CLOSED meshes.
		 */
		Eigen::MatrixXd V_shifted;
		Eigen::Vector3d cen = Eigen::Vector3d(0,0,0);
		double vol;
		igl::centroid(V, F, cen,vol); 

		// [2] apply COM to origin translation
		V_shifted = V.rowwise() - cen.transpose();

		// [3] calculate max distance of all points in V_shifted to origin (0,0,0)
		Eigen::Vector3d origin = Eigen::Vector3d(0,0,0); 
		double max_dist = std::numeric_limits<double>::min();
		for(int i = 0; i < V_shifted.rows(); ++i)
		{
			Eigen::Vector3d cur_point = V.row(i);
			double dist = (cur_point - origin).norm();
			max_dist = std::max(max_dist,dist);
		}

		// [4] scale every point by max_dist 
		V_scaled = V_shifted / max_dist;
		F_scaled = F;
	}

	void generateMeshFragment(Eigen::MatrixXd& VA, Eigen::MatrixXi& FA, 
							Eigen::MatrixXd& VB, Eigen::MatrixXi& FB,
							Eigen::MatrixXd& VC, Eigen::MatrixXi& FC)
	{
		// Bool opers = { "Union", "Intersect", "Minus", "XOR", "Resolve" }
		Eigen::VectorXi J;
		igl::MeshBooleanType boolean_type(igl::MESH_BOOLEAN_TYPE_MINUS); 

		Eigen::MatrixXd tempV;
		Eigen::MatrixXi tempF;
		igl::copyleft::cgal::mesh_boolean(VA,FA,VB,FB,boolean_type,tempV,tempF,J);

		Eigen::MatrixXd v_boundary;
		Eigen::MatrixXi f_boundary;
		extractBoundary(VA, FA, J, VC, FC);
		return;

	}

	/*
	* NOTE :: J is of form [FA:FA.rows() + FB] = birth facet order
	* NOTE :: J is an index into a face, but does NOT actually correspond to a face of interest. 
	*/
	void extractBoundary(Eigen::MatrixXd& VA, Eigen::MatrixXi& FA, Eigen::VectorXi& J, 
						Eigen::MatrixXd& bV, Eigen::MatrixXi& bF)
	{
		// loop over and get total # of faces w/corresponding vertices
		int faceCnt = 0;
		for(int i = 0; i < J.size(); ++i)
		{	
			if(J(i) <= FA.rows()) // yah, since vertex is 1-based indexing
				faceCnt++;
		}	

		// Keep a portion of the mesh in A, and chuck out portion in B
		bF = Eigen::MatrixXi(faceCnt,3);
		int idx = 0;
		for(int i = 0; i < J.size(); ++i)
		{	
			int faceIdx = J(i);
			if(faceIdx <= FA.rows())
			{
				Eigen::Vector3i myFace = FA.row(faceIdx);
				bF.row(idx) = myFace;
				idx++;
			}
		}

		// #TODO :: handle dangling vertices. 
		// Assert that they DONT crash existing alignment algortihm.
		// IF so, then you need to fix this.
		bV = VA;
	}
}

/* 
* 
* EXTRANEUOUS INFORMATION NOTED 
* #include <igl/circumradius.h> - not desired. This is a per triangle measure.
* #include <igl/inradius.h> - not desired. Also a per triangle measure.
* #include <igl/all_pairs_distances.h> // could use this to cheat ( take max of vector? )
* #include <igl/min.h> // #TODO :: use this in place of for-loop dist calculations?
* #include <igl/max.h>
* #include <igl/boundary_facets.h> - this determines boundary faces, but for tetrahedra; boundary edges, for triangles. We could use this, but it's too much work.
* #include <igl/boundary_edges.h> - Akin to <igl/boundary_facets.h>. Also too much work to use.
*/


/*
* OLD testing code for <unit_sphere_scaling>. Place elsewhere.
* void unitSphereTests()
* {
std::string rescaled_mesh_a_file_name = "rescaled_A.obj";
std::string rescaled_mesh_b_right_file_name = "rescaled_right_B.obj";
std::string rescaled_mesh_b_left_file_name = "rescaled_left_B.obj";
igl::writeOBJ(rescaled_mesh_a_file_name, VA_scaled,FA);
igl::writeOBJ(rescaled_mesh_b_left_file_name, VB_left_scaled,FB);
igl::writeOBJ(rescaled_mesh_b_right_file_name, VB_right_scaled,FB);
}
*/

// #TODO :: NOTE THIS CODE TOO!
// note that here, <C> is actually for coloring purposes. 
// # NOTE :: this is the technique Vouga was talking about for the remeshed vs interp thing. ... oh, coloring would be really helpful, now that I think about it.
/*
Eigen::MatrixXd C(FC.rows(),3);
for(size_t f = 0;f<C.rows();f++)
{
if(J(f)<FA.rows())
{
	C.row(f) = Eigen::RowVector3d(1,0,0);
}
else
{
	C.row(f) = Eigen::RowVector3d(0,1,0);
}
}
*/

/*
 * NEW CODING PATTERNS :: 
 * 	[1] dir for shared inputs
 * 	[2] dir for your specific outputs
 * 	[3] One method handles all IO operations. Do not leave IO strewn across everywhere.
 */

/* A WAY TO EXTRACT PLANAR BOUNDARIES! WIL NO LONGER BE USED
std::set<int> onPlane;   
	for(int i = 0; i < VC.rows(); ++i)
	{
		Vector3d myRow = VC.row(i);
		if(myRow(0) == 0)
			onPlane.insert(i);
	}

	for(int i = 0; i < FC.rows(); ++i)
	{
		Vector3i myFace = FC.row(i);
		int v1 = myFace(0);
		int v2 = myFace(1);
		int v3 = myFace(2);
		const bool is_in_v1 = onPlane.find(v1) != onPlane.end();
		const bool is_in_v2 = onPlane.find(v2) != onPlane.end();
		const bool is_in_v3 = onPlane.find(v3) != onPlane.end();

		if(is_in_v1 && is_in_v2 && is_in_v3) {
			faceCnt++;
		}
	}


	// so let's just generate the mesh from this instead
	// need to set up our mapping here too!
	// or the cheap solution - hey, let's just create isolated points :-P
	Eigen::MatrixXi bF = Eigen::MatrixXi(faceCnt,3);
	int idx = 0;
	for(int i = 0; i < FC.rows(); ++i)
	{
		Vector3i myFace = FC.row(i);
		int v1 = myFace(0);
		int v2 = myFace(1);
		int v3 = myFace(2);
		const bool is_in_v1 = onPlane.find(v1) != onPlane.end();
		const bool is_in_v2 = onPlane.find(v2) != onPlane.end();
		const bool is_in_v3 = onPlane.find(v3) != onPlane.end();

		if(is_in_v1 && is_in_v2 && is_in_v3) {
			bF.row(idx) = myFace;
			idx++;
		}
	}
*/

/*
 * OLD CODE - RECTANGUALR PRISM SCALING
	//obj2.V = obj2_preSc.V;
	//obj2.F = obj2_preSc.F;

	// We're basically making a really huge block; emulating an \inf plane, in a sense
	for(int i = 0; i < obj2.V.rows(); ++i)
	{
		Vector3d temp = obj2.V.row(i);
		temp(2) = temp(2) * 100;
		temp(1) = temp(1) * 100;
		obj2.V.row(i) = temp;
 *	}
 */
