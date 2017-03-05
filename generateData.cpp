#include <igl/readOFF.h> 
#include <igl/viewer/Viewer.h>
#include "tutorial_shared_path.h"
#include <igl/rotation_matrix_from_directions.h>
#include <igl/writeOFF.h>
#include "glob_defs.h"

// #TODO :: should I generate a couple of nice VIEWERS, to ook at all the different stuf? not sure? 

Eigen::MatrixXd V;
Eigen::MatrixXi F;

Eigen::MatrixXd Vt;
Eigen::MatrixXi Ft;

/*
bool key_down( igl::viewer::Viewer& viewer, unsigned char key, int modifier)
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
    viewer.data.set_mesh(Vt,Ft);
    viewer.core.align_camera_center(Vt,Ft); 
  } 
  return false;
}
*/

 /// tkae note that the viewer.core.align_camera_center, will make it look like both plane scans are exactly the same ( u got horribly scorched on this a while ago! TAKE NOTE ) 

/*
int main(int argc, char *argv[])
{
  // Load a mesh in OFF format
  igl::readOFF(GLOBAL::interpSurfGenScan1File,V,F);
  std::cout << R"(
1 switch to identity view
2 Switch to rotated view 1
    )";

	Eigen::Vector3d e1 (1,0,0);
	Eigen::Vector3d e2 (0,0,1);
	Eigen::Matrix3d rot_matrix_yAxis = igl::rotation_matrix_from_directions(e1,e2);
    Eigen::MatrixXd Vt_temp;
  Vt_temp = (V * rot_matrix_yAxis) * rot_matrix_yAxis;

	// apply translation vector operation
	// gonna take a stab here, and guess that the mirror plane, has to be moved in the positive z-axis direction ( say, by 0.5 units )
	// otherwise, it'll just straight up be in the negative dir 
  
	// I was hoping that this would work 
	//double tz = 0.5;
  //Vt.col(3) += tz;  

	Eigen::Vector3d translateZ (0,0,-0.5);
	Vt = Vt_temp.rowwise() + translateZ.transpose(); 

  Ft = F; 

  igl::writeOFF(GLOBAL::interpSurfGenScan2File,Vt,Ft);

  // Plot the ROTATED mesh 
  igl::viewer::Viewer viewer;
  viewer.callback_key_down = &key_down;
  viewer.data.set_mesh(V, F);
  viewer.launch();
}
	*/




