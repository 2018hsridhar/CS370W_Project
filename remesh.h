#ifndef REMESH_H
#define REMESH_h

#include "tutorial_shared_path.h"

// LibIGL includes
#include <igl/readOFF.h>
#include <igl/writeOFF.h>
#include <igl/viewer/Viewer.h>
#include <igl/avg_edge_length.h>
#include <igl/cat.h>

// CGAL includes
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <boost/function_output_iterator.hpp>

// C++ includes
#include <fstream>
#include <vector>

// using namespace ( should not exist! ) 
using namespace Eigen; 

namespace REMESH 
{
	bool remeshSurface(Eigen::MatrixXd& V, Eigen::MatrixXi& F, 
						Eigen::MatrixXd& Vr, Eigen::MatrixXi& Fr, 
						double target_edge_length);
	double avgEdgeLenInputMeshes(Eigen::MatrixXd& V1, Eigen::MatrixXi& F1, Eigen::MatrixXd& V2, Eigen::MatrixXi& F2);

}

#endif
