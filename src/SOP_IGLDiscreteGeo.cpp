#include <igl/avg_edge_length.h>
#include <igl/cotmatrix.h>
#include <igl/invert_diag.h>
#include <igl/massmatrix.h>
#include <igl/parula.h>

#include <igl/barycenter.h>
#include <igl/grad.h>
#include <igl/jet.h>

#include <igl/per_corner_normals.h>
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/principal_curvature.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <GU/GU_Detail.h>
#include <GU/GU_PrimPoly.h>

#include "SOP_IGLDiscreteGeo.hpp"

#include <UT/UT_DSOVersion.h>
#include <GU/GU_Detail.h>
#include <OP/OP_Operator.h>
#include <OP/OP_AutoLockInputs.h>
#include <OP/OP_OperatorTable.h>
#include <PRM/PRM_Include.h>
#include <UT/UT_Matrix3.h>
#include <UT/UT_Matrix4.h>
#include <SYS/SYS_Math.h>


// #include "stdafx.h"
#include <stddef.h>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


using namespace SOP_IGL;

void
newSopOperator(OP_OperatorTable *table)
{
    table->addOperator(new OP_Operator(
        "IGLDiscreteGeo",
        "IGL Discrete Geometry Operator",
        SOP_IGLDiscreteGeometry::myConstructor,
        SOP_IGLDiscreteGeometry::myTemplateList,
        1,
        1,
        0));
}


static PRM_Name names[] = {
    PRM_Name("curvature",          "Add Principal Curvature"),
    PRM_Name("false_curve_colors", "Add False Curve Colors"),
    PRM_Name("grad_attrib",        "Add Gradient of Attribute (scalar)"),    
    PRM_Name("grad_attrib_name",   "Scalar Attribute Name"),
    PRM_Name("laplacian",          "Laplacian (Smoothing)"),
    PRM_Name("eigenvectors",       "Eigen Decomposition"),
};



PRM_Template
SOP_IGLDiscreteGeometry::myTemplateList[] = {
    PRM_Template(PRM_TOGGLE, 1, &names[0], PRMzeroDefaults),
    PRM_Template(PRM_TOGGLE, 1, &names[1], PRMzeroDefaults),
    PRM_Template(PRM_TOGGLE, 1, &names[2], PRMzeroDefaults),
    PRM_Template(PRM_STRING, 1, &names[3], 0),
    PRM_Template(PRM_INT_J,  1, &names[4], PRMzeroDefaults),
    PRM_Template(PRM_TOGGLE, 1, &names[5], PRMzeroDefaults),
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
    gdp->convex(); // only triangles for now.
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
    int laplacian_iterations = LAPLACIAN(t);

    if (laplacian_iterations != 0)
    {
        // Start the interrupt server
        UT_AutoInterrupt boss("Laplacian smoothing...");
        Eigen::SparseMatrix<double> L;
        // Compute Laplace-Beltrami operator: #V by #V
        igl::cotmatrix(V,F,L);
        // Smoothing:
        Eigen::MatrixXd U; U = V;

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
        }

        // Copy back to Houdini:
        GA_Offset ptoff;
        GA_FOR_ALL_PTOFF(gdp, ptoff) 
        {
            const GA_Index ptidx = gdp->pointIndex(ptoff);
            if ((uint)ptidx < U.rows()) 
            {
                const UT_Vector3 pos(U((uint)ptidx, 0),
                                     U((uint)ptidx, 1),
                                     U((uint)ptidx, 2));
                gdp->setPos3(ptoff, pos);
            }
        }
    }
    
    // FIXME: igs reports error here: gs::eigs(L, M, k+1, igs::EIGS_TYPE_SM, U, D)
    #if 0 
     /*  Eigen decompositon*/
    if (EIGENVECTORS(t)) 
    {
        int c=0;
        double bbd = 1;
        bool twod = 0;
  
        twod = V.col(2).minCoeff() == V.col(2).maxCoeff();
        bbd = (V.colwise().maxCoeff() - V.colwise().minCoeff()).norm();

        Eigen::SparseMatrix<double> L, M;
        igs::cotmatrix(V,F,L);
        L = (-L).eval();

        igs::massmatrix(V, F, igs::MASSMATRIX_TYPE_DEFAULT, M);

        const size_t k = 5;
        Eigen::VectorXd D;
        if(!igs::eigs(L, M, k+1, igs::EIGS_TYPE_SM, U, D)) {
             addWarning(SOP_MESSAGE, "Can't compute eigen decomposition.");
            return error();
        }

        U = ((U.array()-U.minCoeff())/(U.maxCoeff()-U.minCoeff())).eval();

        //....
    }

    #endif



    

    gdp->getP()->bumpDataId();
    return error();
}
