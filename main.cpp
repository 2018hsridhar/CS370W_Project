// NEW pattern
	// [1] dir for shared inputs
	// [1] dir for your specific outputs
// a good pattern to estabslih
// ... have one method handle all IO operations
// ... it can grow confusing to have too many methods handle IO operations

// #TODO :: so with CSG operations, I can generate seperate meshes. BUT ... I need to get their boundary. And I'm not sure how to do this. I.e. I'm really just interested in a specific plane.
	// ... wait a sec, isn't there some sort of simple vertex test here, or what not? I feel like there is. 
	// ... well, at least I know what to look for next now.
	// ... or one could apply a projection ( onto a 2D plane ), right?
	// ... wait, isn't there a physical sim project where we did this a while ago?

// #TODO :: we'll apply code logic for 2 or 3 splits. we need 5 generic vars for 2 splits, 7 generic vars for 3 splits [ 1 = original mesh, 1-1 = cutting meshes/resultant cut mesh ].

// USER_DEFINED LIBRARIES
#include "bool_opers.h" // primary include
#include "glob_defs.h"
#include "interpSurface.h"
#include "remesh.h"
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

struct Mesh
{
	Eigen::MatrixXd V; 
	Eigen::MatrixXi F;
	Eigen::MatrixXd bV;  // boundary vertices
	Eigen::VectorXi bI;  // boundary indices - for specific faces [ NOT a {0,1} boolean thing ] 
	Eigen::MatrixXd N;   // Normals ( for Faces? or Vertices? )
	int numBV; 		 	 // number of boundary vertices
} meshPreScaled, mesh, meshL, meshR, frag1, frag2, scan1,scan2, interp, remeshed;

const char * MESH_BOOLEAN_TYPE_NAMES[] =
{
	"Union",
	"Intersect",
	"Minus",
	"XOR",
	"Resolve",
};

int main(int argc, char *argv[])
{
    igl::readOBJ(GLOBAL::boolOpersMesh,meshPreScaled.V,meshPreScaled.F);
    applyUnitSphereRescaling(meshPreScaled.V, meshPreScaled.F, mesh.V);
	mesh.F = meshPreScaled.F;
	generateBoolOperMeshes(mesh.V, mesh.F, meshL.V, meshL.F, meshR.V, meshR.F);
	generateMeshFragment(mesh.V,mesh.F,meshL.V,meshL.F, frag1.V, frag1.F);
	generateMeshFragment(mesh.V,mesh.F,meshR.V,meshR.F, frag2.V, frag2.F);
    igl::writeOBJ(GLOBAL::boolOpersFragOne,frag1.V,frag1.F);
    igl::writeOBJ(GLOBAL::boolOpersFragTwo,frag2.V,frag2.F);
}

/*
 * To note about boolean operations 
 * YES, meshes must be water tight for bool operations to be applicable. 
 *    they must also be orientable.
 * Hence why they fail in cases such as <circle.off> but not <camelHead.off>
 * NOTE :: V1 = left, V2 = right #TODO update naming conventions.
 */
void generateBoolOperMeshes(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::MatrixXd& V1, Eigen::MatrixXi& F1, Eigen::MatrixXd& V2, Eigen::MatrixXi& F2)
{
    // Rescale both objects to unit sphere
    Eigen::MatrixXd VA_scaled;
    Eigen::MatrixXd VB_pre_offset_scaled;
    applyUnitSphereRescaling(V, F, VB_pre_offset_scaled);

	// APPLY offsets to VB, both [-|+] x-axis. Generate two additional meshes for cutting
	Eigen::Vector3d offset_left = Eigen::Vector3d(-0.5,0,0);
	Eigen::Vector3d offset_right = Eigen::Vector3d(0.5,0,0);

	V1 = VB_pre_offset_scaled.rowwise() + offset_left.transpose();
	V2 = VB_pre_offset_scaled.rowwise() + offset_right.transpose();
	F1 = F2 = F;
}


// #TODO:: include as a test method, run in a testing system ( 
// Method responsibility - asserts preservation of order of vertex indexing between the interpolating and remeshed surfaces. 
// equivalently, assserts that the n vertices (interp) = first n ( remesh )
// ... ( e.g. akin to bb test, from Amazon ) 
void runIndexTestInterpRemesh()
{
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
}


/* Idea taken from calculating smallest bounding sphere for a mesh
 *     and scaling according to its radius
 *     should be noted that for meshes contained in the unit sphere, it'll still work
 *     just based on a fractional radius
 *     Note that operations are limited to watertight meshes
 */

void applyUnitSphereRescaling(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::MatrixXd& V_scaled)
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
}

// NOTE :: this method really can be placed back elsewhere, but for now, kep it in case we have to extend it later on!
void generateMeshFragment(Eigen::MatrixXd& VA, Eigen::MatrixXi& FA, 
							Eigen::MatrixXd& VB, Eigen::MatrixXi& FB,
							Eigen::MatrixXd& VC, Eigen::MatrixXi& FC)
{
	Eigen::VectorXi J;
	igl::MeshBooleanType boolean_type(igl::MESH_BOOLEAN_TYPE_INTERSECT);
	// #TODO :: why is seg fault happening here?
	igl::copyleft::cgal::mesh_boolean(VA,FA,VB,FB,boolean_type,VC,FC,J);
	
	// note that here, <C> is actually for coloring purposes. We can eliminate this later, if need be!
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
}

/*
bool key_down(igl::viewer::Viewer &viewer, unsigned char key, int mods)
{
  switch(key)
  {
    default:
      return false;
    case '.':
      boolean_type =
        static_cast<igl::MeshBooleanType>(
          (boolean_type+1)% igl::NUM_MESH_BOOLEAN_TYPES);
      break;
    case ',':
      boolean_type =
        static_cast<igl::MeshBooleanType>(
          (boolean_type+igl::NUM_MESH_BOOLEAN_TYPES-1)%
          igl::NUM_MESH_BOOLEAN_TYPES);
      break;
    case '[':
      viewer.core.camera_dnear -= 0.1;
      return true;
    case ']':
      viewer.core.camera_dnear += 0.1;
      return true;
  }
  applyBooleanOperations()
  return true;
}
*/


/* 
 * 
 * EXTRANEUOUS INFORMATION NOTED 
 * #include <igl/circumradius.h> - not desired. This is a per triangle measure.
 * #include <igl/inradius.h> - not desired. Also a per triangle measure.
 * #include <igl/all_pairs_distances.h> // could use this to cheat ( take max of vector? )
 * #include <igl/min.h> // #TODO :: use this in place of for-loop dist calculations?
 * #include <igl/max.h>
 *
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


