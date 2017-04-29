#ifndef J_ALIGN_H
#define J_ALIGN_H

#include "tutorial_shared_path.h"
#include <set>
#include <time.h> 

// LIBIGL includes


// possibly needed ... will have to check on this later
#include <igl/random_dir.h> 
#include <igl/axis_angle_to_quat.h>
#include <igl/quat_to_mat.h>
#include <igl/doublearea.h>

#include <igl/per_vertex_normals.h>
// what the heck is a weighting scheme, on the vertex normals, anyways??

// J_ALIGN = Impulse based alignment of rigid bodies
using namespace Eigen; 

/* @ each timestep ... solve for energy, then apply rigid motion
 * calculations ( for entities, such as normals ) must be performed @ each timestep
 * ... should the method be broken up for each body ... or fufill for both?
 * ... idea ... each body ( feed in the rigidBodyInst and rigidBodyTemplates )  
 * ... @ every step, I expect to recieve an outputted "transformation" matrix too!
 * ... do I need to store the initial config ( to get resulting transformation matrix too? ) ... probably!
 * ... any useful debug methods? or debug statements to put, if Vouga needs to look at this? ... possibly!
 * ... templating and instance  code must be updated ... my files are not tetrahedrons!
 * I should consider my entire set of faces, as boundary faces ... for my two inpuits, is 3D data, "rolled out" as a 2D plane
 * I wonder how useful ~<...>() destructors are, in the process of testing or rigid body manipulation? Perhaps you can see things that shouldn't be getting deleted, get deleted ( or overwritten possibly, based on implemented safety measures)?
 * Do I still need to sort those faces, for the inertia tensor constructor? unsure atm
 * It should also be noted, that my "energy" value, is only for one trans mat @ each timestep!
 * there is no need to descend on the energy values here! ... but you do WANT to assert that the energy is decreasing over time ( if not ... that is a bug! ) 
* --- cutoff @ a low enough energy ( for now, do < 0.1)
* --- the set up for this algorithm REMAINS HIGHLY CRITICAL!!! BE CAREFUL HERE!!! ... need to update thy configuration vector, over time!
 */

namespace J_ALIGN 
{
	void calculateExternalForces();
	void calculateConfigForce();
	void applyRigidMotion();	
	void convertConfigToTransformationMatrix();
	void updateConfiguration();
	void assembleConfiguration();
	void dissassembleConfiguration();	
}

#endif
