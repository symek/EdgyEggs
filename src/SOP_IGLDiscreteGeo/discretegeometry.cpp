#include <igl/cotmatrix.h>
#include <igl/invert_diag.h>
#include <igl/massmatrix.h>
#include <igl/principal_curvature.h>
#include <igl/grad.h>
#include <igl/parula.h>
#include <igl/eigs.h>
#include <igl/exact_geodesic.h>

#include <GU/GU_Detail.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "discretegeometry.hpp"

namespace SOP_IGL {

void compute_curvature(GU_Detail *gdp, Eigen::MatrixXd &V, Eigen::MatrixXi &F, const int false_colors)
{
    // Alternative discrete mean curvature
    Eigen::MatrixXd HN;
    Eigen::SparseMatrix<double> L,M,Minv;
    igl::cotmatrix(V,F,L);
    igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_VORONOI,M);
    igl::invert_diag(M,Minv);
    // Laplace-Beltrami of position
    HN = -Minv*(L*V);
    // Extract magnitude as mean curvature
    Eigen::VectorXd H = HN.rowwise().norm();

    // Compute curvature directions via quadric fitting
    Eigen::MatrixXd PD1,PD2;
    Eigen::VectorXd PV1,PV2;
    igl::principal_curvature(V,F,PD1,PD2,PV1,PV2);
    // mean curvature
    H = 0.5*(PV1+PV2);
   

    GA_RWHandleF      curvature_h(gdp->addFloatTuple(GA_ATTRIB_POINT, "curvature", 1));
    GA_RWHandleV3     tangentu_h(gdp->addFloatTuple(GA_ATTRIB_POINT, "tangentu", 3));
    GA_RWHandleV3     tangentv_h(gdp->addFloatTuple(GA_ATTRIB_POINT, "tangentv", 3));

    GA_Offset ptoff;
    if (curvature_h.isValid() && tangentu_h.isValid() && tangentv_h.isValid()) {
        GA_FOR_ALL_PTOFF(gdp, ptoff) {
            const GA_Index ptidx = gdp->pointIndex(ptoff);
            UT_ASSERT((uint)ptidx < H.rows());
            UT_ASSERT((uint)ptidx < PD1.rows());
            UT_ASSERT((uint)ptidx < PD2.rows());
            const float curv = H((uint)ptidx, 0);
            const UT_Vector3 tnu(PD1((uint)ptidx, 0), PD1((uint)ptidx, 1), PD1((uint)ptidx, 2));
            const UT_Vector3 tnv(PD2((uint)ptidx, 0), PD2((uint)ptidx, 1), PD2((uint)ptidx, 2));
            curvature_h.set(ptoff, curv);
            tangentu_h.set(ptoff, tnu);
            tangentv_h.set(ptoff, tnv);
        }
    }

    if (false_colors) {
        // Pseudo Colors:
        Eigen::MatrixXd C;
        igl::parula(H,true,C);
        GA_RWHandleV3   color_h(gdp->addFloatTuple(GA_ATTRIB_POINT, "Cd", 3));
        GA_FOR_ALL_PTOFF(gdp, ptoff) { 
            const GA_Index ptidx = gdp->pointIndex(ptoff);  
            UT_ASSERT((uint)ptidx < C.rows());
            const UT_Vector3 cd(C((uint)ptidx, 0), C((uint)ptidx, 1), C((uint)ptidx, 2));
            color_h.set(ptoff, cd);
        }
    }
}

void compute_gradient(GU_Detail *gdp, const GA_ROHandleF &sourceAttrib, Eigen::MatrixXd &V, Eigen::MatrixXi &F)
{
    const uint numPoints = gdp->getNumPoints();
    Eigen::VectorXd U(numPoints);
    GA_Offset ptoff;
    GA_FOR_ALL_PTOFF(gdp, ptoff) {
        const float val      = sourceAttrib.get(ptoff);
        const GA_Index ptidx = gdp->pointIndex(ptoff);
        U((uint)ptidx) = val;
    }

    // Compute gradient operator: #F*3 by #V
    Eigen::SparseMatrix<double> G;
    igl::grad(V,F,G);

    // Compute gradient of U
    Eigen::MatrixXd GU = Eigen::Map<const Eigen::MatrixXd>((G*U).eval().data(),F.rows(),3);
    // Compute gradient magnitude
    const Eigen::VectorXd GU_mag = GU.rowwise().norm();


    // Copy attributes to Houdini
    { 
        GA_RWHandleV3  gradAttrib_h(gdp->addFloatTuple(GA_ATTRIB_POINT, "gradientAttrib", 3));
        GA_Offset ptoff;
        GA_FOR_ALL_PTOFF(gdp, ptoff) { 
            const GA_Index ptidx = gdp->pointIndex(ptoff);  
            UT_ASSERT((uint)ptidx < GU_mag.rows());
            UT_ASSERT((uint)ptidx < GU.rows());
            UT_Vector3 grad(GU((uint)ptidx, 0),  GU((uint)ptidx, 1), GU((uint)ptidx, 2));
            const float gmag = GU_mag((uint)ptidx, 0);
            grad *= gmag;
            grad.normalize();//?
            gradAttrib_h.set(ptoff, grad);
        }
    }
}

int compute_laplacian(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, 
                      const Eigen::SparseMatrix<double> &L, Eigen::MatrixXd &U)
{
    Eigen::SparseMatrix<double> M;
    igl::massmatrix(U,F,igl::MASSMATRIX_TYPE_BARYCENTRIC, M);

    // Solve (M-delta*L) U = M*U
    const auto & S = (M - 0.001*L);
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double > > solver(S);

    if(solver.info() != Eigen::Success)
      return solver.info();
    U = solver.solve(M*U).eval();
    return Eigen::Success;
}


bool compute_exact_geodesic_distance(
	const Eigen::MatrixXd & points, 
	const Eigen::VectorXi & faces,
	const int anchor,
	      Eigen::VectorXd & distances)
{
	Eigen::VectorXi VS,FS,VT,FT;
    VS.resize(1);
    VS << anchor;
    // All vertices are the targets
    VT.setLinSpaced(points.rows(), 0, points.rows()-1);
    igl::exact_geodesic(points, faces, VS, FS, VT, FT, distances);
    return true;
}


} // end of SOP_IGL namespace