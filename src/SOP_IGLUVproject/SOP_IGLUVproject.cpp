#include <igl/arap.h>
#include <igl/boundary_loop.h>
#include <igl/harmonic.h>
#include <igl/map_vertices_to_circle.h>

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
#include "SOP_IGLUVproject.hpp"

using namespace SOP_IGL;

void
newSopOperator(OP_OperatorTable *table)
{
    table->addOperator(new OP_Operator(
        "igluvproject",
        "EE IGLUVProject",
        SOP_IGLUVproject::myConstructor,
        SOP_IGLUVproject::myTemplateList,
        1,
        1,
        0));
}

static PRM_Name names[] = {
    PRM_Name("term",     "UV Projection"),
    PRM_Name("maxiter",  "Max interations (ARAP)"),
};

static PRM_Name  termChoices[] =
{
    PRM_Name("0", "ARAP"),
    PRM_Name("1", "LSCM"),
    PRM_Name(0)
};

static PRM_Default     iterDefault(5);
static PRM_ChoiceList  termMenu(PRM_CHOICELIST_SINGLE,  termChoices);
static PRM_Range       maxiterRange(PRM_RANGE_PRM, 0, PRM_RANGE_UI, 25);

PRM_Template
SOP_IGLUVproject::myTemplateList[] = {
    PRM_Template(PRM_ORD,   1, &names[0], 0, &termMenu, 0, 0),
    PRM_Template(PRM_INT_J, 1, &names[1], &iterDefault, 0, &maxiterRange),
    PRM_Template(),
};


OP_Node *
SOP_IGLUVproject::myConstructor(OP_Network *net, const char *name, OP_Operator *op)
{
    return new SOP_IGLUVproject(net, name, op);
}

SOP_IGLUVproject::SOP_IGLUVproject(OP_Network *net, const char *name, OP_Operator *op)
    : SOP_Node(net, name, op), myGroup(NULL)
{
   
    mySopFlags.setManagesDataIDs(true);
}

SOP_IGLUVproject::~SOP_IGLUVproject() {}

OP_ERROR
SOP_IGLUVproject::cookInputGroups(OP_Context &context, int alone)
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
SOP_IGLUVproject::cookMySop(OP_Context &context)
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


    const int uvterm  = TERM(t);
    const int maxiter = MAXITER(t);

    if (uvterm == 0)
    {
        Eigen::MatrixXd V_uv;
        Eigen::MatrixXd initial_guess;

        // Compute the initial solution for ARAP (harmonic parametrization)
        Eigen::VectorXi bnd;
        igl::boundary_loop(F, bnd);

        if (bnd.rows() == 0) {
            addWarning(SOP_MESSAGE, "No boundaries, can't proceed.");
            return error();
        }

        Eigen::MatrixXd bnd_uv;
        igl::map_vertices_to_circle(V, bnd, bnd_uv);

        // This won't compile with Eigen > 3.2.8
        #if 1

        if (!igl::harmonic(V,F,bnd,bnd_uv,1,initial_guess))
        {
            addWarning(SOP_MESSAGE, "Can't compute harmonics.");
            return error();
        }

        // Add dynamic regularization to avoid to specify boundary conditions
        igl::ARAPData arap_data;
        arap_data.with_dynamics = true;
        Eigen::VectorXi b  = Eigen::VectorXi::Zero(0);
        Eigen::MatrixXd bc = Eigen::MatrixXd::Zero(0,0);

        // Initialize ARAP
        arap_data.max_iter = maxiter;
        // 2 means that we're going to *solve* in 2d
        if (!arap_precomputation(V,F,2,b,arap_data))
        {
            addWarning(SOP_MESSAGE, "Can't precompute ARAP.");
            return error();
        }

        // Solve arap using the harmonic map as initial guess
        V_uv = initial_guess;

        if (!arap_solve(bc,arap_data,V_uv))
        { 
            addWarning(SOP_MESSAGE, "Can't solve ARAP.");
            return error();
        }

        GA_RWHandleV3 uv_h(gdp->addFloatTuple(GA_ATTRIB_POINT, "uv", 3));
        GA_Offset ptoff;
        GA_FOR_ALL_PTOFF(gdp, ptoff) 
        { 
            const GA_Index ptidx = gdp->pointIndex(ptoff);  
            UT_ASSERT((uint)ptidx < V_uv.rows());
            const UT_Vector3 uv(V_uv((uint)ptidx, 0), V_uv((uint)ptidx, 1), 0.f);
            uv_h.set(ptoff, uv);
        }
    }

    // etc
    #endif

   

    gdp->getP()->bumpDataId();
    return error();
}
