#include <igl/readOFF.h>
#include "tutorial_shared_path.h"
#include <igl/viewer/Viewer.h>

// MY LIBS 
#include "meanCurvatureFlow.h"

using namespace Eigen; 
using namespace std;

Eigen::MatrixXd V;
Eigen::MatrixXi F;

int main(int argc, char *argv[])
{

    // LOAD mesh data ( OFF format )
    //igl::readOFF(TUTORIAL_SHARED_PATH "/camelhead.off",V,F);
    igl::readOFF(TUTORIAL_SHARED_PATH "/cow.off",V,F);
    bool hasBndry = MCF::meshHasBoundary(V,F);
    if(hasBndry)
    {
        std::cout << "INPUT Mesh does have a boundary.\n";
    }
    // PLOT initial mesh 
    igl::viewer::Viewer viewer;

	// ... how to update viewer, over a timestep ... that would be nice !
	// ... see Clement's approach to MCF code ... sure, that works, I guess ... 

    //viewer.callback_key_down = &key_down;
    //viewer.data.set_mesh(V, F);
	Eigen::MatrixXd Vc;
	MCF::computeMeanCurvatureFlow(V,F,0.001, Vc);
	viewer.data.set_mesh(Vc,F);
    viewer.launch();

    return 0;
}
