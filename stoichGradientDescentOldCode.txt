// Q1 :: What is the order of convergence of this method? not sure atm
// purpose of this code :: given a (4x4) transformation matrix/alignemnt, iteratively improve the alignemnt, via stoichastic gradient descent !  try 100 of these things ... take min. it will most likely converge ! 

#include <igl/readSTL.h>
#include <igl/readOFF.h>
#include <igl/readOBJ.h>

#include <igl/writeSTL.h>
#include <igl/viewer/Viewer.h>
#include "tutorial_shared_path.h"
#include <igl/rotation_matrix_from_directions.h>
#include <igl/writeOFF.h>
#include <igl/per_vertex_normals.h>
#include <igl/point_mesh_squared_distance.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h> 
#include <igl/doublearea.h>


// other include files, at the moment
#include <igl/random_dir.h> 
#include <igl/axis_angle_to_quat.h>
#include <igl/quat_to_mat.h>
#include <set>
#include <time.h> // \#TODO :: do I even need this  ??
//#include <random.h> // c++11 rand val generator  ( note :: do not actually use rand() ! ) 

using namespace Eigen; 
using namespace std;

int main(int argc, char *argv[])
{
	// Develop an initial transformation matrix, filled with random data 
	Eigen::Matrix4d T = Eigen::Matrix4d::Identity();

	// REFERNCE THIS FOR RANDOMIZER :: http://stackoverflow.com/questions/288739/generate-random-numbers-uniformly-over-an-entire-range/20136256#20136256
	// generate 100 random (3x1) vector of small {e_x,e_y,e_z} values

	const double range_from  = -0.001;
	const double range_to    = 0.001;
	std::random_device                  rand_dev;
	std::mt19937                        generator(rand_dev());
	std::uniform_real_distribution<double>  distr(range_from, range_to);

	int i;
	//for(i = 0; i < 100; i++)
	for(i = 0; i < 1; i++)
	{

		// generator values for an epsilon vector 
		double e_x = distr(generator);
		double e_y = distr(generator);
		double e_z = distr(generator);

		Eigen::Vector3d epsilon_vector;
		epsilon_vector[0] = e_x;
		epsilon_vector[1] = e_y;
		epsilon_vector[2] = e_z;

		double epsilon_angle = distr(generator); // is there adifferent wahy to generator the angle? 
		double epsilon_axis[] = { e_x,e_y,e_z};
        
        // note :: i might want to work with TrackBall.cpp's approach ( see /include/igl )
		// generate (4x4) transformation matrix, via ( axis-angle) => quat => matrix pipeline
		// or should I actually try this out? 
		Eigen::Vector3d epsilon_axis2 = igl::random_dir(); // this might be a better idea ! and then convert to an arr?? 
		double qrot[4];
		igl::axis_angle_to_quat(epsilon_axis,epsilon_angle,qrot);
		double mat[16];
		igl::quat_to_mat(qrot, mat);

        Eigen::Matrix4d rotatePerturbation = Eigen::Matrix4d::Identity();  // note :: translate vector is zeroed out here! 
        for ( unsigned i = 0; i < 4; ++i )
            for ( unsigned j = 0; j < 4; ++j )
               rotatePerturbation(i,j) = mat[i+4*j]; 

        Eigen::Matrix4d translatePerturbation = Eigen::Matrix4d::Zero();
        // I should be able to rewrite this! 
        translatePerturbation(0,3) = e_x;
        translatePerturbation(1,3) = e_y;
        translatePerturbation(2,3) = e_z;
    
		// generate perturbed ( Rotation, Translation, RotateTrans ) matrices. 
        Eigen::Matrix4d petrubedRotation =  T * rotatePerturbation;
        Eigen::Matrix4d perturbedTranslation = T + translatePerturbation;
        Eigen::Matrix4d perturbedRotAndTrans =  (T * rotatePerturbation) + translatePerturbation;

       // 

	}
}

/* initialize random seed: 
randVals = rand(3,1)*7 - 10;
randVals = 10.^randVals;
 * srand (time(NULL));
 * generate secret number between -7 and 2: 
 * int randVal = rand() % 10 - 8;
 * double epsilonVal = 10 ^ randVal;
 */












