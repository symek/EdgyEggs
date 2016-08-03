#pragma once

#include <GA/GA_SplittableRange.h>
#include <GA/GA_Range.h>
#include <GA/GA_PageIterator.h>
#include <GA/GA_PageHandle.h>
#include <SOP/SOP_Node.h>

#include <time.h>

#include "converters.hpp"

namespace SOP_IGL {

class Timer 
{
    private:
        double begTime;
    public:
        void start()
        {
            begTime = clock();
        }

        double current() 
        {
        //int threads = UT_Thread::getNumProcessors();
            return (difftime(clock(), begTime) / CLOCKS_PER_SEC);// / (double) threads;
        }

        bool isTimeout(double seconds)
        {
            return seconds >= current();
        }
};


void compute_curvature(GU_Detail *gdp, Eigen::MatrixXd &V, Eigen::MatrixXi &F, const int false_colors=0)
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
        GA_RWHandleV3  gradAttrib_h(gdp->addFloatTuple(GA_ATTRIB_POINT, "gradient", 3));
        GA_Offset ptoff;
        GA_FOR_ALL_PTOFF(gdp, ptoff) { 
            const GA_Index ptidx = gdp->pointIndex(ptoff);  
            UT_ASSERT((uint)ptidx < GU_mag.rows());
            UT_ASSERT((uint)ptidx < GU.rows());
            UT_Vector3 grad(GU((uint)ptidx, 0),  GU((uint)ptidx, 1), GU((uint)ptidx, 2));
            const float gmag = GU_mag((uint)ptidx, 0);
            grad *= gmag;
            gradAttrib_h.set(ptoff, grad);
        }
    }
}

class SOP_IGLDiscreteGeometry : public SOP_Node
{
public:
         SOP_IGLDiscreteGeometry(OP_Network *net, const char *name, OP_Operator *op);
    virtual ~SOP_IGLDiscreteGeometry();

    static PRM_Template      myTemplateList[];
    static OP_Node      *myConstructor(OP_Network*, const char *,
                                OP_Operator *);

    /// This method is created so that it can be called by handles.  It only
    /// cooks the input group of this SOP.  The geometry in this group is
    /// the only geometry manipulated by this SOP.
    virtual OP_ERROR         cookInputGroups(OP_Context &context, 
                        int alone = 0);

protected:
    /// Method to cook geometry for the SOP
    virtual OP_ERROR         cookMySop(OP_Context &context);

private:
    void    getGroups(UT_String &str)        { evalString(str, "group", 0, 0); }
    int     CURVATURE(fpreal t)              { return evalInt("curvature", 0, t); }
    int     FALSE_CURVE_COLORS(fpreal t)     { return evalInt("false_curve_colors", 0, t); }
    int     GRAD_ATTRIB(fpreal t)            { return evalInt("grad_attrib", 0, t); }
    void    GRAD_ATTRIB_NAME(UT_String &str) { evalString(str,"grad_attrib_name", 0, 0); }
    fpreal  LAPLACIAN(fpreal t)              { return evalFloat("laplacian", 0, t); }
    // fpreal  QCOEF(fpreal t)          { return evalFloat("qcoef", 0, t); }
    // fpreal  ZCOEF(fpreal t)           { return evalFloat("zcoef", 0, t); }
    // fpreal  RADIUS(fpreal t)    { return evalFloat("radius", 0, t); }
    // int     LAYERS(fpreal t)    { return evalInt("layers", 0, t); }
    // fpreal  LAMBDA(fpreal t)    { return evalFloat("lambda", 0, t); }
    // int     TANGENT(fpreal t)   { return evalInt("tangent", 0, t); }

    /// This is the group of geometry to be manipulated by this SOP and cooked
    /// by the method "cookInputGroups".
    const GA_PointGroup *myGroup;
};


} // End SOP_IGL namespace


// class op_IGLDiscreteGeometry {
// public:
//     op_IGLDiscreteGeometry(const std::string &str_model, const bool tnSpace, GU_Detail *gdp)
//         : mystr_model(str_model),  myGdp(gdp), myTnSpace(tnSpace) {};
//             // Take a SplittableRange (not a GA_Range)
//     void    operator()(const GA_SplittableRange &r) const
//             {
//                 alglib::rbfmodel model;
//                 GA_RWPageHandleV3 handle_P(myGdp->getP());
//                 GA_ROPageHandleV3 handle_U(myGdp, GA_ATTRIB_POINT, "tangentu");
//                 GA_ROPageHandleV3 handle_V(myGdp, GA_ATTRIB_POINT, "tangentv");
//                 GA_ROPageHandleV3 handle_N(myGdp, GA_ATTRIB_POINT, "N");

//                 // Execute storage:
//                 alglib::real_1d_array coord("[0,0,0]");
//                 alglib::real_1d_array result("[0,0,0]");
//                 std::string copy_serialized(mystr_model);
//                 alglib::rbfunserialize(copy_serialized, model);
                
//                 // Iterate over pages in the range
//                 for (GA_PageIterator pit = r.beginPages(); !pit.atEnd(); ++pit)
//                 {
//                     GA_Offset start, end;
//                     // iterate over the elements in the page.
//                     for (GA_Iterator it(pit.begin()); it.blockAdvance(start, end); )
//                     {
//                         // Perform any per-page setup required, then
//                         handle_P.setPage(start); handle_U.setPage(start);
//                         handle_V.setPage(start); handle_N.setPage(start);
//                         for (GA_Offset i = start; i < end; ++i)
//                         {
//                             const UT_Vector3 pos = handle_P.get(i);
//                             const double dp[] = {pos.x(), pos.y(), pos.z()};
//                             coord.setcontent(3, dp);
//                             alglib::rbfcalc(model, coord, result);
//                             UT_Vector3 displace = UT_Vector3(result[0], result[1], result[2]);
//                             if (myTnSpace)
//                             {
//                                 UT_Vector3 u = handle_U.get(i); 
//                                 UT_Vector3 v = handle_V.get(i);
//                                 UT_Vector3 n = handle_N.get(i);
//                                 u.normalize(); v.normalize(); n.normalize();
//                                 UT_Matrix3  b(u.x(), u.y(), u.z(),
//                                               v.x(), v.y(), v.z(),
//                                               n.x(), n.y(), n.z());

//                                 b = b.transposedCopy() * b;
//                                 UT_Vector3 a1(u * b); a1.normalize();
//                                 UT_Vector3 a2(v * b); a2.normalize();
//                                 const float da1 = displace.dot(a1);
//                                 const float da2 = displace.dot(a2);
//                                 displace        = UT_Vector3(a1 * da1 + a2 * da2);
        
//                             }

//                             handle_P.set(i, pos+displace);
//                         }
//                     }
//                 }
//             }
//     private:
//             const std::string &mystr_model;
//             GU_Detail         *myGdp;
//             const bool        myTnSpace; 
// };
// void
// IGLDiscreteGeometryThreaded(const GA_Range &range, const std::string &str_model, \
//     const bool tnSpace, GU_Detail *gdp)
// {
//     // Create a GA_SplittableRange from the original range
//     GA_SplittableRange split_range = GA_SplittableRange(range);
//     UTparallelFor(split_range, op_IGLDiscreteGeometry(str_model, tnSpace, gdp));
// }
