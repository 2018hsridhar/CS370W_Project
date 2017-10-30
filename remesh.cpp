#include "remesh.h"
#include "glob_defs.h"
using namespace Eigen;

// CGAL typedefs
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;
typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<Mesh>::edge_descriptor     edge_descriptor;
namespace PMP = CGAL::Polygon_mesh_processing;

//const char* filename =  TUTORIAL_SHARED_PATH "/pig.off";
//const char* outfilename = TUTORIAL_SHARED_PATH "/pigOutput.off";

// HalfEdge struct
struct halfedge2edge
{
  halfedge2edge(const Mesh& m, std::vector<edge_descriptor>& edges)
    : m_mesh(m), m_edges(edges)
  {}
  void operator()(const halfedge_descriptor& h) const
  {
    m_edges.push_back(edge(h, m_mesh));
  }
  const Mesh& m_mesh;
  std::vector<edge_descriptor>& m_edges;
}; 
// Preallocation of V,F
/*
Eigen::MatrixXd V;
Eigen::MatrixXi F;

Eigen::MatrixXd Vt;
Eigen::MatrixXi Ft;

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

// I will write a function, assuming that {V,F} are given, and return a {Vt,Ft}
// crud ... is there like, a temporary file ( for storage ), or a way in CGAl for me to output Eigen matrices?? that would prove convenient!
/*
int main(int argc, char* argv[])
{
// the CGAL portion!
// Load a mesh in OFF format
  igl::readOFF(GLOBAL::remeshInputFile, V, F); 
  std::cout << R"(
1 Switch to original view.
2 Switch to triagnulated poly mesh view.
    )";
	REMESH::remeshSurface(V,F,Vt,Ft);

	// Plot the original mesh 
	igl::viewer::Viewer viewer;
	viewer.callback_key_down = &key_down;
	viewer.data.set_mesh(V, F);
	viewer.launch();

  return 0;
}
*/

namespace REMESH
{ 
	bool remeshSurface(Eigen::MatrixXd& V, Eigen::MatrixXi& F, 
							Eigen::MatrixXd& Vr, Eigen::MatrixXi& Fr,
							double target_edge_length)
	{
		bool badInterp = true;
		igl::writeOFF(GLOBAL::remeshInputFile,V,F);
		std::ifstream input(GLOBAL::remeshInputFile);
	//	std::ifstream input(TUTORIAL_SHARED_PATH "/badOffset.off");
		Mesh mesh;
		if (!input || !(input >> mesh)) {
			std::cerr << "Not a valid off file." << std::endl;
			return false;
		}
		unsigned int nb_iter = 3;
		//std::cout << "Split border...";
		std::vector<edge_descriptor> border;
		PMP::border_halfedges(faces(mesh),
			mesh,
			boost::make_function_output_iterator(halfedge2edge(mesh, border)));
		PMP::split_long_edges(border, target_edge_length, mesh);
		//std::cout << "done." << std::endl;
		//std::cout << "Start remeshing of " << GLOBAL::remeshInputFile 
		//	<< " (" << num_faces(mesh) << " faces)..." << std::endl;
		if(num_faces(mesh) == 0)
		{
			std::cout << "REMESHING WILL FAIL HERE" << std::endl;
			badInterp = false;
		}
		// something in isotropic_remeshing seems to fail ... occasional floating_point_exception?
		// and why is the number of faces so low ( 116 vs 232 ) in the camel case? That too, is weird!!
		PMP::isotropic_remeshing(
			faces(mesh),
			target_edge_length,
			mesh,
			PMP::parameters::number_of_iterations(nb_iter)
			.protect_constraints(true)//i.e. protect border, here
			);
	//	std::cout << "Remeshing done." << std::endl;
		std::ofstream cube_off(GLOBAL::remeshOutputFile);
		cube_off << mesh;
		cube_off.close();
		igl::readOFF(GLOBAL::remeshOutputFile, Vr,Fr);
		return badInterp;
	}

	double avgEdgeLenInputMeshes(Eigen::MatrixXd& V1, Eigen::MatrixXi& F1,
							 	Eigen::MatrixXd& V2, Eigen::MatrixXi& F2)
	{
		Eigen::MatrixXd V_cat;
		Eigen::MatrixXi F_cat;
		igl::cat(1,V1,V2,V_cat);
		igl::cat(1,F1,MatrixXi(F2.array() + V1.rows()), F_cat);	
		return igl::avg_edge_length(V_cat,F_cat);
	}

	/*
	* Gets boundary vertices for both parts of the remeshed interpolated surface
	* Based on solving for <boundary_loop> of this remeshed interpolated surface 
	*/
	//void getRemeshedBoundaryVerts(const int numBV_one, const int numBV_two, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::MatrixXd& remeshBoundaryOne, Eigen::MatrixXd& remeshBoundaryTwo)
	void getRemeshedBoundaryVerts(const int numBV_one, const int numBV_two, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, std::vector<int>& remeshBoundaryOne, std::vector<int>& remeshBoundaryTwo)
	{
		// ITERATE over <resmeshed.V>
		// KNOW that the first [boundaryOne,boundaryTwo] vertices are already set.
		// CONSTRUCT boundary loops for the remeshsed surface 
		Eigen::MatrixXd remeshBndry;
		std::vector<std::vector<int>> bndryLoop; // #NOTE :: thsi is the appropriate call
		// Refer to <hessian_energy.cpp> for your example.

		// GET boundary_loops for the remeshed vertex
		// The indexing scheme of both files makes our lives super easy! 
		// Alongside the fact that this useful method exists too!
		igl::boundary_loop(F, bndryLoop);

		// GET both boundary loop vectors
		std::vector<int> firstBL = bndryLoop[0];
		std::vector<int> secondBL = bndryLoop[1];

		// ASSERT that all verts in boundaryOne are in the firstBL 
		// [scan1 = boundary1], and [scan2 = boundary2] here
		int b1_s_ptr = 0;
		int b2_s_ptr = numBV_one;
		while(b1_s_ptr < numBV_one)
		{
			// assert vector contains said index
			if(std::find(firstBL.begin(), firstBL.end(), b1_s_ptr) == firstBL.end()) 
			{
				std::cout << "Fail to contain a vertex in boundary 1" << std::endl;
				break;
			}
			b1_s_ptr++;
		}

		// ASSERT that all verts in boundaryTwo are in the secondBL 
		while(b2_s_ptr < numBV_two)
		{
			// assert vector contains said index
			if(std::find(secondBL.begin(), secondBL.end(), b2_s_ptr) == secondBL.end()) 
			{
				std::cout << "Fail to contain a vertex in boundary 2" << std::endl;
				break;
			}
			b2_s_ptr++;
		}

		remeshBoundaryOne = firstBL;
		remeshBoundaryTwo = secondBL;

		// CONVERT type :: firstBL -> Eigen::VectorXi ( it is a vector, after all )
/*
		Eigen::VectorXi bl_one;
		bl_one.resize(firstBL.size());
		for(int i = 0; i < firstBL.size(); ++i)
			bl_one(i) = firstBL[i];

		Eigen::VectorXi bl_two;
		bl_two.resize(secondBL.size());
		for(int i = 0; i < secondBL.size(); ++i)
			bl_two(i) = secondBL[i];

		// NOW get those two boundaries!
		igl::slice(V,bl_one,1,remeshBoundaryOne);
		igl::slice(V,bl_two,1,remeshBoundaryTwo);
*/		
	}
}
