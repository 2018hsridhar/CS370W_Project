/*
#include <igl/readOFF.h>
#include <igl/viewer/Viewer.h>
#include "tutorial_shared_path.h"
#include <igl/rotation_matrix_from_directions.h>
#include <igl/writeOFF.h>
#include <igl/per_vertex_normals.h>
using namespace Eigen; // I can only wonder why this managed to fix the issue ??? weird ??  
using namespace std;

// I guess this is how you pass matrices as arguements AROUND in EIGEN !!!

// note :: this is ONLY FOR POINTS ... if they were vectors ... then I need a different method !
Eigen::MatrixXf convertToHomogenousForm(const Ref<const MatrixXf>& mat )
{
  Eigen::MatrixXf homogenoizedMatrix;
  homogenoizedMatrix = mat.transpose().colwise().homogeneous().transpose(); 
  return homogenoizedMatrix;
}

Eigen::MatrixXf normalizeHomogenousMatrix (const Ref<const MatrixXf>& mat )
{
  Eigen::MatrixXf normalizedMatrix;
  normalizedMatrix = mat.transpose().colwise().hnormalized().transpose(); 
  return normalizedMatrix;
}

int main(int argc, char *argv[])
{

  // note :: do store a couple of test matrices ... they seem to be useful for l8r purposes :-) 
  Eigen::Matrix3f m = Matrix3f::Random();
  std::cout << m << std::endl;

  Eigen::MatrixXf homogenous = convertToHomogenousForm(m);
  std::cout << "\n\n\n" << std::endl;
  std::cout << homogenous << std::endl;

  Eigen::MatrixXf deHomogenous = normalizeHomogenousMatrix(homogenous);
  std::cout << "\n\n\n" << std::endl;
  std::cout << deHomogenous << std::endl;
}
*/




