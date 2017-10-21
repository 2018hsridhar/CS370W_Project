// #TODO :: so with CSG operations, I can generate seperate meshes. BUT ... I need to get their boundary. And I'm not sure how to do this. I.e. I'm really just interested in a specific plane.
	// ... wait a sec, isn't there some sort of simple vertex test here, or what not? I feel like there is. 
	// ... well, at least I know what to look for next now.
	// ... or one could apply a projection ( onto a 2D plane ), right?
	// ... wait, isn't there a physical sim project where we did this a while ago?

// #TODO :: we'll apply code logic for 2 or 3 splits. we need 5 generic vars for 2 splits, 7 generic vars for 3 splits [ 1 = original mesh, 1-1 = cutting meshes/resultant cut mesh ].
// #TODO :: refactor. This code is a mess of tests.
// ... okay, do review those parts later. BUT realize that at the moment, getting these CSGs generated is priority at the moment!

#include <igl/readOFF.h>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/viewer/Viewer.h>

#include <Eigen/Core>
#include <iostream>

#include "tutorial_shared_path.h"

// includes for CSG/bool opers AND unit-sphere scaling
// #TODO ... wait don't I have a useful max distance library? ... might be useful for later!
#include <igl/copyleft/cgal/mesh_boolean.h>
#include <igl/doublearea.h>
 #include <igl/centroid.h> // actually, closed != water tight. This should work
// #include <igl/circumradius.h> - not desired. This is a per triangle measure.
// #include <igl/inradius.h> - not desired. Also a per triangle measure.
#include <igl/all_pairs_distances.h> // could use this to cheat ( take max of vector? )
#include <igl/min.h> // #TODO :: use this to easily cheat too?
#include <igl/max.h>

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


/*
 ********************************************************************************
 ********************************************************************************
 ********************************************************************************
 ********************************************************************************
 ********************************************************************************
 ********************************************************************************
 */

using namespace Eigen;  
using namespace std;
using namespace igl;

struct Mesh
{
	Eigen::MatrixXd V; 
	Eigen::MatrixXd bV;  // boundary vertices
	Eigen::MatrixXi F;
	Eigen::VectorXi bI; // boundary indices; just indxs to specific faces ( not a 1,0 boolean thing ) 
	Eigen::MatrixXd N;
	int numBV; 			// number of boundary vertices
} scan1,scan2,scan1_rest,scan2_rest,interp,remeshed, flowed, debug_result, result, normals, normals_info;

// #TODO :: look up <circumradius.cpp> library? 

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

void applyUnitSphereRescaling(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::MatrixXd& V_scaled);
void applyBooleanOperations();

// #TODO :: take out viewer logic later. For now it's useful since we want to see these cuts getting generated.
//void applyBooleanOperations()
void applyBooleanOperations(igl::viewer::Viewer &viewer)
{
  //std::cout << "before bool operations app" << std::endl;
	// getting thrown a <Input mesh is not orientable!> error here! 
    // error is being thrown under /cgal/propagate_winding_numbers.cpp. Winding number is off for the planes case. Need to try cubes instead? 
  igl::copyleft::cgal::mesh_boolean(VA,FA,VB,FB,boolean_type,VC,FC,J);
  //std::cout << "after bool operations app" << std::endl;
  Eigen::MatrixXd C(FC.rows(),3);
  for(size_t f = 0;f<C.rows();f++)
  {
    if(J(f)<FA.rows())
    {
      C.row(f) = Eigen::RowVector3d(1,0,0);
    }else
    {
      C.row(f) = Eigen::RowVector3d(0,1,0);
    }
  }
  /*
  viewer.data.clear();
  viewer.data.set_mesh(VC,FC);
  viewer.data.set_colors(C);
  */ 
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
  applyBooleanOperations(viewer);
  return true;
}
*/


void runIndexTestInterpRemesh()
{
// check indexing structure of meshes ( simple test ) 
	if(!readOFF(GLOBAL::pipelineScan1File,scan1.V,scan1.F)) {
		std::cout << "Failed to load partial scan one." << std::endl;
	} 

	if(!readOFF(GLOBAL::pipelineScan2File,scan2.V,scan2.F)) {
		std::cout << "Failed to load partial scan two." << std::endl;
	}

	std::cout << "Executing pipeline for the following meshes" << std::endl;
	std::cout << "Mesh one [" << GLOBAL::pipelineScan1File << "]" << std::endl;
	std::cout << "Mesh two [" << GLOBAL::pipelineScan2File << "]" << std::endl;

// check if first n ( interp ) = first n ( remesh )  

	INTERP_SURF::generateInterpolatingSurface(scan1.V,scan1.F,scan2.V,scan2.F, interp.V,interp.F);
	double rEL = REMESH::avgEdgeLenInputMeshes(scan1.V,scan1.F,scan2.V,scan2.F);
	rEL = rEL / 5.0; // #NOTE :: can actually visualize pinching occuring here
	bool remSucc = REMESH::remeshSurface(interp.V,interp.F,remeshed.V,remeshed.F, rEL);
	if(!remSucc)
	{
		std::cout << "BAD REMESHING!" << std::endl;
		exit(0);
	}

	igl::writeOFF("interpolating_surface.off", interp.V,interp.F);
	igl::writeOFF("remeshed_surface.off", remeshed.V,remeshed.F);
}


// so the cool thing, is that with unit sphere scaling
// and knowing that I can apply row-wise translations
// I should be able to easily generate some nice mesh data to use here :-)
// ... so let's see if I can get both the cube cuts and sphere cuts tonight :-). That'll be a good start and motivator! 

// need to rename some of these variables -> see support structures above!

int main(int argc, char *argv[])
{
    using namespace Eigen;
    using namespace std;
	// note :: your input meshes MUST be orientable. Else err
	// so I can get the intersection ... but how to get two pieces not part of intersection?A
	// wiith the plane ... intersection and unions do NOT work as expected. Kinda weird!

	/*
     ***************** To note about boolen operations *******************
	 * YES, meshes must be water tight for bool operations to be applicable. 
	 * Hence why they fail in cases such as "circle.off" but not "camelHead.off"
     */

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
}

// first check if libigl already provides you [unit-sphere] scaling funcs. Those might prove useful.
// ... unfortunately, this is NOT provided at all! need to rethink approach here!
// getting the mesh's COM is decently easy
// can I get a <bounding_sphere>. Need to find out more!
// NOTE :: we can keep faces the same. We care only for vertices

// #TODO :: make this a boolean to indicate <success/failure>? 
// NOTE :: thsi method will also account if our object is contained WITHIN the unit sphere. Division by a fractional amount will be performed here.
// NOTE ::these operations are being applied on WATERTIGHT meshes. so we expect to be able to calculate a centroid here.

void applyUnitSphereRescaling(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::MatrixXd& V_scaled)
{
	// [1] calcualte COM ( = centroid ) 
	// ... wait a sec, can't we just use the centroid since this is a const density obj?
	// ... not :: centroid calculations work for CLOSED MESHES ONLY! They won't work for your two inputs, unfortunately, since you deal with boundaries!
	Eigen::MatrixXd V_shifted;
	Eigen::Vector3d cen(0,0,0); // centroid = (3x1) dimensional vector  ... might've needed this
	double vol;
	igl::centroid(V, F, cen,vol); // why is this centroid full of Nans? 

	// [2] apply transformation
	V_shifted = V.rowwise() - cen.transpose();
   
	// [3] calculate max distance of all points in V_shifted to origin (0,0,0)
	// honestly, it's easier for me to just code up a for loop over the entirety of V_shifted.
	Eigen::Vector3d origin = Eigen::Vector3d(0,0,0); 
	double max_dist = std::numeric_limits<double>::min();
	for(int i = 0; i < V_shifted.rows(); ++i)
	{
		Eigen::Vector3d cur_point = V.row(i);
		double dist = (cur_point - origin).norm();
		max_dist = std::max(max_dist,dist);
	}

	// [4] scale every point by max_dist #TODO :: a functional way to write this up?
	// #TODO :: ensure correctness of unit sphere rescaling operations here BTW!
	// 		--- would be a good idea to just straight up visualize this!
	// #TODO :: is this the current point of failure? possibly? need to find out more though!
	V_scaled = V_shifted / max_dist;
}
