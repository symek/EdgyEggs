#include <igl/cotmatrix.h>
#include <igl/invert_diag.h>
#include <igl/massmatrix.h>
#include <igl/principal_curvature.h>
#include <igl/grad.h>
#include <igl/parula.h>
#include <igl/eigs.h>
#include <igl/exact_geodesic.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <UT/UT_DSOVersion.h>
#include <GU/GU_Detail.h>
#include <GU/GU_PrimPoly.h>
#include <OP/OP_Operator.h>
#include <OP/OP_AutoLockInputs.h>
#include <OP/OP_OperatorTable.h>
#include <PRM/PRM_Include.h>
#include <UT/UT_Matrix3.h>
#include <UT/UT_Matrix4.h>
#include <SYS/SYS_Math.h>

#include "../converters.hpp"
#include "SOP_IGLDiscreteGeo.hpp"
#include "discretegeometry.hpp"


using namespace SOP_IGL;

void
newSopOperator(OP_OperatorTable *table)
{
    table->addOperator(new OP_Operator(
        "igldiscretegeo",
        "EE IGL Discrete Geometry",
        SOP_IGLDiscreteGeometry::myConstructor,
        SOP_IGLDiscreteGeometry::myTemplateList,
        1,
        2,
        0));
}


static PRM_Name names[] = {
    PRM_Name("curvature",          "Principal Curvature"),
    PRM_Name("false_curve_colors", "False Curve Colors"),
    PRM_Name("grad_attrib",        "Gradient of Attribute (scalar)"),    
    PRM_Name("grad_attrib_name",   "Scalar Attribute Name"),
    PRM_Name("laplacian",          "Laplacian (smoothing/sharpening)"),
    PRM_Name("eigenvectors",       "Eigen Decomposition "),
    PRM_Name("geodesicdistance",   "Geodesic distance"),
};

static PRM_Range  laplaceRange(PRM_RANGE_PRM, 0, PRM_RANGE_PRM, 10);

PRM_Template
SOP_IGLDiscreteGeometry::myTemplateList[] = {
    PRM_Template(PRM_TOGGLE, 1, &names[0], PRMzeroDefaults),
    PRM_Template(PRM_TOGGLE, 1, &names[1], PRMzeroDefaults),
    PRM_Template(PRM_TOGGLE, 1, &names[2], PRMzeroDefaults),
    PRM_Template(PRM_STRING, 1, &names[3], 0),
    PRM_Template(PRM_INT_J,  1, &names[4], PRMzeroDefaults, 0, &laplaceRange),
    PRM_Template(PRM_INT_J , 1, &names[5], PRMzeroDefaults, 0, &laplaceRange),
    PRM_Template(PRM_TOGGLE, 1, &names[6], PRMzeroDefaults),
    PRM_Template(),
};


OP_Node *
SOP_IGLDiscreteGeometry::myConstructor(OP_Network *net, const char *name, OP_Operator *op)
{
    return new SOP_IGLDiscreteGeometry(net, name, op);
}

SOP_IGLDiscreteGeometry::SOP_IGLDiscreteGeometry(OP_Network *net, const char *name, OP_Operator *op)
    : SOP_Node(net, name, op), myGroup(NULL)
{
   
    mySopFlags.setManagesDataIDs(true);
}

SOP_IGLDiscreteGeometry::~SOP_IGLDiscreteGeometry() {}

OP_ERROR
SOP_IGLDiscreteGeometry::cookInputGroups(OP_Context &context, int alone)
{
    
    return cookInputPointGroups(
        context, // This is needed for cooking the group parameter, and cooking the input if alone.
        myGroup, // The group (or NULL) is written to myGroup if not alone.
        alone,   // This is true iff called outside of cookMySop to update handles.
                 // true means the group will be for the input geometry.
                 // false means the group will be for gdp (the working/output geometry).
        true,    // (default) true means to set the selection to the group if not alone and the highlight flag is on.
        0,       // (default) Parameter index of the group field
        -1,      // (default) Parameter index of the group type field (-1 since there isn't one)
        true,    // (default) true means that a pointer to an existing group is okay; false means group is always new.
        false,   // (default) false means new groups should be unordered; true means new groups should be ordered.
        true,    // (default) true means that all new groups should be detached, so not owned by the detail;
                 //           false means that new point and primitive groups on gdp will be owned by gdp.
        0        // (default) Index of the input whose geometry the group will be made for if alone.
    );
}

OP_ERROR
SOP_IGLDiscreteGeometry::cookMySop(OP_Context &context)
{
    
    OP_AutoLockInputs inputs(this);
    if (inputs.lock(context) >= UT_ERROR_ABORT)
        return error();

    fpreal t = context.getTime();
    duplicateSource(0, context);

    // Copy to eigen.
    gdp->convex(); // only triangles for now, but point count will match
    uint numPoints = gdp->getNumPoints();
    uint numPrims  = gdp->getNumPrimitives();
    Eigen::MatrixXd V(numPoints, 3); // points
    Eigen::MatrixXi F(numPrims, 3); // faces

    SOP_IGL::detail_to_eigen(*gdp, V, F);

    /* Curvature */
    if(CURVATURE(t)) 
    {
        SOP_IGL::compute_curvature(gdp, V, F, FALSE_CURVE_COLORS(t));
    }

    /* Gradient of an attribute */
    if (GRAD_ATTRIB(t)) 
    {

        // Fetch our attribute names
        UT_String sourceGradAttribName;
        GRAD_ATTRIB_NAME(sourceGradAttribName);

        if (!sourceGradAttribName.isstring())
            return error();
      
        // Make sure a name is valid.
        sourceGradAttribName.forceValidVariableName();
        GA_ROHandleF  sourceGradAttrib_h(gdp->findAttribute(GA_ATTRIB_POINT, sourceGradAttribName));
        if (sourceGradAttrib_h.isValid()) {
            SOP_IGL::compute_gradient(gdp, sourceGradAttrib_h, V, F);
        } else {
            addWarning(SOP_MESSAGE, "Can't compute a gradient from a given attribute.");
        }
    }


    /*    Laplacian smoothing   */
    const float laplacian = LAPLACIAN(t);

    if (laplacian != 0)
    {
        int laplacian_iterations = (int)ceil(laplacian);
        float laplacian_ratio    = laplacian - floorf(laplacian);
        laplacian_ratio = laplacian_ratio != 0.f ? laplacian_ratio : 1.f;

        // Start the interrupt server
        UT_AutoInterrupt boss("Laplacian smoothing...");
        Eigen::SparseMatrix<double> L;
        // Compute Laplace-Beltrami operator: #V by #V
        igl::cotmatrix(V,F,L);
        // Smoothing:
        Eigen::MatrixXd U; U = V;
        Eigen::MatrixXd T;
        T = Eigen::MatrixXd::Zero(V.rows(), V.cols());

        while(laplacian_iterations) 
        {
            // User interaption/
            if (boss.wasInterrupted())
                return error();

            if (SOP_IGL::compute_laplacian(V, F, L, U) != Eigen::Success) {
                addWarning(SOP_MESSAGE, "Can't compute laplacian with current geometry.");
                return error();
            }
            laplacian_iterations--;

            // if (laplacian_iterations > 0)
            //     T += (U - T);
            // else
            //     T += (U - T) * laplacian_ratio;
            T = U;
        }

        // Copy back to Houdini:
        GA_Offset ptoff;
        GA_FOR_ALL_PTOFF(gdp, ptoff) 
        {
            const GA_Index ptidx = gdp->pointIndex(ptoff);
            if ((uint)ptidx < T.rows()) 
            {
                const UT_Vector3 pos(T((uint)ptidx, 0),
                                     T((uint)ptidx, 1),
                                     T((uint)ptidx, 2));
                gdp->setPos3(ptoff, pos);
            }
        }
    }
    
    // Thisize_ts won't compile with Eigen > 3.2.8
    #if 1
    /*  Eigen decompositon*/
    const size_t eigenvectors = EIGENVECTORS(t);
    if (eigenvectors != 0) 
    {
        Eigen::SparseMatrix<double> L, M;
        igl::cotmatrix(V,F,L);
        L = (-L).eval();

        igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);
        int c = 0;
        bool  twod = V.col(2).minCoeff()==V.col(2).maxCoeff();
        double bbd = (V.colwise().maxCoeff()-V.colwise().minCoeff()).norm();

        const size_t k = eigenvectors;
        Eigen::VectorXd D;
        Eigen::MatrixXd U;
        if(!igl::eigs(L, M, k+1, igl::EIGS_TYPE_SM, U, D)) {
            addWarning(SOP_MESSAGE, "Can't compute eigen decomposition.");
            return error();
        }

        U = ((U.array()-U.minCoeff())/(U.maxCoeff()-U.minCoeff())).eval();
        U = U.rightCols(k).eval();
        std::cout << U << '\n';

        // Rescale eigen vectors for visualization
        Eigen::VectorXd Z = bbd*0.5*U.col(c);
        Eigen::MatrixXd C;
        igl::parula(U.col(c).eval(),false,C);
        GA_Attribute * cd_a = gdp->addFloatTuple(GA_ATTRIB_POINT, "Cd", 3);
        GA_RWHandleV3 cd_h(cd_a);
        GA_Offset ptoff;
        GA_FOR_ALL_PTOFF(gdp, ptoff) {
            const GA_Index ptidx = gdp->pointIndex(ptoff);
            if ((uint)ptidx < C.rows()) 
            {
                const UT_Vector3 clr(C((uint)ptidx, 0),
                                     C((uint)ptidx, 1),
                                     C((uint)ptidx, 2));
                cd_h.set(ptoff, clr);
            }
        }
    }

    #endif

    if (evalInt("geodesicdistance", 0, t))  {
        const GA_Attribute * anchors_attr = gdp->findFloatTuple(GA_ATTRIB_POINT, "anchors", 1);
        if (!anchors_attr) {
            addWarning(SOP_MESSAGE, "Create float attribute 'anchors' \
                and set 1 for points in question.");
            return error();
        }

        GA_ROHandleF  anchors_h(anchors_attr);

        int anchor_index = -1;
        if (anchors_h.isValid()) {
            GA_Offset ptoff;
            GA_FOR_ALL_PTOFF(gdp, ptoff) {
                const float is_anchor = anchors_h.get(ptoff);
                if (is_anchor > 0.0) {
                    const GA_Index i = gdp->pointIndex(ptoff);
                    anchor_index = static_cast<int>(i);
                }
            }
        }

        if (anchor_index >= 0) {
            Eigen::VectorXd distances;
            if (!SOP_IGL::compute_exact_geodesic_distance(V, F, anchor_index, distances)) {
                addWarning(SOP_MESSAGE, "Can't compute geodesic distance.");
                return error();
            }

            GA_Attribute * distance_attr = gdp->addFloatTuple(GA_ATTRIB_POINT, "distance", 1);
            GA_RWHandleF distance_h(distance_attr);
            if (distance_h.isInvalid()) {
                addWarning(SOP_MESSAGE, "Can't add distance attribute.");
                return error();
            }

            GA_Offset ptoff;
            GA_FOR_ALL_PTOFF(gdp, ptoff) {
                const GA_Index i = gdp->pointIndex(ptoff);
                const double distance = distances(i);
                distance_h.set(ptoff, (float)distance);
            }
        }
    }


    gdp->getP()->bumpDataId();
    return error();
}
