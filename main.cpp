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
} scan1,scan2,scan1_rest,scan2_rest,interp,remeshed, flowed, debug_result, result, normals, normals_info;

Eigen::MatrixXd VA,VB,VC;
Eigen::VectorXi J,I;
Eigen::MatrixXi FA,FB,FC;
igl::MeshBooleanType boolean_type(igl::MESH_BOOLEAN_TYPE_INTERSECT);

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
	// note :: your input meshes MUST be orientable. Else err
	// so I can get the intersection ... but how to get two pieces not part of intersection?A
	// wiith the plane ... intersection and unions do NOT work as expected. Kinda weird!

	/*
     ***************** To note about boolen operations *******************
	 * YES, meshes must be water tight for bool operations to be applicable. 
	 * Hence why they fail in cases such as "circle.off" but not "camelHead.off"
     */

/*
    igl::readOBJ(TUTORIAL_SHARED_PATH "/cube.obj",VA,FA);
    igl::readOBJ(TUTORIAL_SHARED_PATH "/cube.obj",VB,FB);

    // use code to rescale both objects to unit sphere
    Eigen::MatrixXd VA_scaled;
    Eigen::MatrixXd VB_left_scaled;
    Eigen::MatrixXd VB_right_scaled;
    Eigen::MatrixXd VB_pre_offset_scaled;
    applyUnitSphereRescaling(VA, FA, VA_scaled);
    applyUnitSphereRescaling(VB, FB, VB_pre_offset_scaled);

	// APPLY offset to VB
	Eigen::Vector3d offset_right = Eigen::Vector3d(0.5,0,0);
	Eigen::Vector3d offset_left = Eigen::Vector3d(-0.5,0,0);
	VB_right_scaled = VB_pre_offset_scaled.rowwise() + offset_right.transpose();
	VB_left_scaled = VB_pre_offset_scaled.rowwise() + offset_left.transpose();

    // write off values
	// #TODO :: why is this rescaling_code failing? that's weird? Why is this null?
	std::string rescaled_mesh_a_file_name = "rescaled_A.obj";
	std::string rescaled_mesh_b_right_file_name = "rescaled_right_B.obj";
	std::string rescaled_mesh_b_left_file_name = "rescaled_left_B.obj";
    igl::writeOBJ(rescaled_mesh_a_file_name, VA_scaled,FA);
    igl::writeOBJ(rescaled_mesh_b_left_file_name, VB_left_scaled,FB);
    igl::writeOBJ(rescaled_mesh_b_right_file_name, VB_right_scaled,FB);
	// #TODO :: this code failed to scale the sphere. I'm not exactly sure why. MUST FIX THIS!!

  // Initialize
  //applyBooleanOperations(viewer);
  //applyBooleanOperations();
  //igl::writeOFF(TUTORIAL_SHARED_PATH "/boolean.off", VC, FC);
*/
/*
  viewer.core.show_lines = true;
  viewer.callback_key_down = &key_down;
  viewer.core.camera_dnear = 3.9;
  cout<<
    "Press '.' to switch to next boolean operation type."<<endl<<
    "Press ',' to switch to previous boolean operation type."<<endl<<
    "Press ']' to push near cutting plane away from camera."<<endl<<
    "Press '[' to pull near cutting plane closer to camera."<<endl<<
    "Hint: investigate _inside_ the model to see orientation changes."<<endl;
  viewer.launch();
*/
	return 0;
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


// #TODO :: remove viewer logic later. Convert to stand-alone method
// Handling <Input mesh is not orientable!> errors
//    note that operations can be applied only on watertight meshes.
//void applyBooleanOperations(igl::viewer::Viewer &viewer)
void applyBooleanOperations()
{
	igl::copyleft::cgal::mesh_boolean(VA,FA,VB,FB,boolean_type,VC,FC,J);
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
	//viewer.data.clear();
	//viewer.data.set_mesh(VC,FC);
	//viewer.data.set_colors(C);
	std::cout<<"A "<<MESH_BOOLEAN_TYPE_NAMES[boolean_type]<<" B."<<std::endl;
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
 * EXTRANEUOUS INFORMATION NOTED 
 * #include <igl/circumradius.h> - not desired. This is a per triangle measure.
 * #include <igl/inradius.h> - not desired. Also a per triangle measure.
 * #include <igl/all_pairs_distances.h> // could use this to cheat ( take max of vector? )
 * #include <igl/min.h> // #TODO :: use this in place of for-loop dist calculations?
 * #include <igl/max.h>
 */


