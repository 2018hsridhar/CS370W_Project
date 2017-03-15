// Q1 :: What is the order of convergence of this method? not sure atm
// purpose of this code :: given a (4x4) transformation matrix/alignemnt, iteratively improve the alignemnt, via stoichastic gradient descent !  try 100 of these things ... take min. it will most likely converge ! 

#include "tutorial_shared_path.h"
#include <set>
#include <time.h> // #TODO :: do I even need this  ??

// LIBIGL includes
#include <igl/random_dir.h> 
#include <igl/axis_angle_to_quat.h>
#include <igl/quat_to_mat.h>
#include <igl/doublearea.h>

// MY LIBS
#include "helpers.h"
#include "glob_defs.h"


// SGD = Stoichastic Gradient Descent
namespace SGD
{
	double calculateSurfaceEnergy(Eigen::VectorXd& V, Eigen::VectorXi& F);
	void generateTransMats(Eigen::Matrix4d& input, std::vector<Eigen::Matrix4d>& transMats);
	void findOptimalTransMat(std::vector<Eigen::Matrix4d>& transMats, std::vector<double>& energies, Eigen::Matrix4d& opt);
	// assert that the size of the two vectors are the same
	// assert that energy is never not positive -- it is a func of S.A. atm, so > 0, should hold true ( or is it >= 0 ?? )
}
