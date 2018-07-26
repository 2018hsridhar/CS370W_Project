#include <igl/readOFF.h>
#include <igl/viewer/Viewer.h>
#include "path.h"
#include "glob_defs.h"
#include <igl/rotation_matrix_from_directions.h>
#include <igl/writeOFF.h>

// initial [ the 1st ] scan
Eigen::MatrixXd V;
Eigen::MatrixXi F;

// reflecting [ the 2nd ] scan
Eigen::MatrixXd V_rotated;
Eigen::MatrixXi F_rotated;

// centering the camera. Said translation will not be visble in the viewer.
bool key_down_distinct( igl::viewer::Viewer& viewer, unsigned char key, int modifier)
{
  std::cout << "Key : " << key << (unsigned int) key << std::endl;
  if ( key == '1' )
  {
    viewer.data.clear();
    viewer.data.set_mesh(V,F);
    viewer.core.align_camera_center(V,F); 
  }
  else if ( key == '2' ) 
  {
    viewer.data.clear();
    viewer.data.set_mesh(V_rotated,F_rotated);
    viewer.core.align_camera_center(V_rotated,F_rotated); 
  } 
  return false;
}

//int main(int argc, char *argv[])
void generateTestData()
{
  // Load the first mesh, in <OFF> format
  igl::readOFF(GLOBAL::pipelineScan1File, V, F);
  std::cout << R"(
1 switch to identity view
2 Switch to rotated view 1
    )";

  // generate rotation along x-axis ( multiply -z * y ) 
  Eigen::Vector3d e1 ( 0, 0, -1 ); 
  Eigen::Vector3d e2 ( 0, 1, 0 ); 

  // one particular issue atm :: this rotation always works by 90-degree angles !!
  Eigen::MatrixXd rot_matrix = igl::rotation_matrix_from_directions (e1, e2 );
  Eigen::MatrixXd transformMat = (rot_matrix * rot_matrix);
  Eigen::VectorXd translate_yaxis = Eigen::Vector3d(0,-1,0);  

  // V_rotated = V * transformMat;
  V_rotated = (V*transformMat).rowwise() + translate_yaxis.transpose(); //
  F_rotated = F; 
  std::string secondScan = "camelhead2.off";
  igl::writeOFF(secondScan, V_rotated,F_rotated);

  // Plot the ROTATED mesh 
  igl::viewer::Viewer viewer;
  viewer.callback_key_down = &key_down_distinct;
  viewer.data.set_mesh(V, F);
  viewer.launch();

}


