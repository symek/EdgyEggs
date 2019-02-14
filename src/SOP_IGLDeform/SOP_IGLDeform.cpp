#include <igl/arap.h>
#include <igl/biharmonic_coordinates.h>
#include <igl/cat.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/matrix_to_list.h>
#include <igl/parula.h>
#include <igl/point_mesh_squared_distance.h>
#include <igl/remove_unreferenced.h>
#include <igl/slice.h>


#include <Eigen/Geometry>
#include <Eigen/StdVector>


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

#include <vector>
#include <algorithm>
#include <iostream>


#include "../converters.hpp"
#include "SOP_IGLDeform.hpp"


using namespace SOP_IGL;

void
newSopOperator(OP_OperatorTable *table)
{
    table->addOperator(new OP_Operator(
        "igldeform",
        "EE IGLDeform",
        SOP_IGLDeform::myConstructor,
        SOP_IGLDeform::myTemplateList,
        2,
        3,
        0));
}

static PRM_Name names[] = {
    PRM_Name("method",   "Deformation Method"),
    PRM_Name("energy",   "ARAP Energy"),
    PRM_Name("maxiter",  "Max interations"),
};

static PRM_Name  methodChoices[] =
{
    PRM_Name("0", "Biharmonic"),
    PRM_Name("1", "ARAP"),
    PRM_Name(0)
};

static PRM_Name  energyChoices[] =
{
    PRM_Name("0", "Spokes"),
    PRM_Name("1", "Spokes and rims"),
    PRM_Name("2", "Elements"),
    PRM_Name("3", "Auto"),
    PRM_Name(0)
};

static PRM_Default     iterDefault(5);
static PRM_ChoiceList  methodMenu(PRM_CHOICELIST_SINGLE,  methodChoices);
static PRM_ChoiceList  energyMenu(PRM_CHOICELIST_SINGLE,  energyChoices);
static PRM_Range       maxiterRange(PRM_RANGE_UI, 0, PRM_RANGE_UI, 50);

PRM_Template
SOP_IGLDeform::myTemplateList[] = {
    PRM_Template(PRM_ORD,   1, &names[0], 0, &methodMenu, 0, 0),
    PRM_Template(PRM_ORD,   1, &names[1], 0, &energyMenu, 0, 0),
    PRM_Template(PRM_INT_J, 1, &names[2], &iterDefault, 0, &maxiterRange),
    PRM_Template(),
};

OP_Node *
SOP_IGLDeform::myConstructor(OP_Network *net, const char *name, OP_Operator *op)
{
    return new SOP_IGLDeform(net, name, op);
}

SOP_IGLDeform::SOP_IGLDeform(OP_Network *net, const char *name, OP_Operator *op)
    : SOP_Node(net, name, op), myGroup(NULL)
{
   
    mySopFlags.setManagesDataIDs(true);
}

SOP_IGLDeform::~SOP_IGLDeform() {}

OP_ERROR
SOP_IGLDeform::cookInputGroups(OP_Context &context, int alone)
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
SOP_IGLDeform::cookMySop(OP_Context &context)
{
    OP_AutoLockInputs inputs(this);
    if (inputs.lock(context) >= UT_ERROR_ABORT)
        return error();

    fpreal t = context.getTime();
    duplicateSource(0, context);

    // Copy to eigen.
    gdp->convex(); // only triangles for now.

    const int method = METHOD(t);
    if (method == AS_RIGID_AS_POSSIBLE) // I'll keep it simple for now.
    {
        const GU_Detail *deform_gdp = inputGeo(1);
        if (gdp->getNumPoints() != deform_gdp->getNumPoints())
        {
            addError(SOP_ERR_MISMATCH_POINT, "Rest and deform geometry should match.");
            return error();
        }
        
        uint numPoints = gdp->getNumPoints();
        uint numPrims  = gdp->getNumPrimitives();

        Eigen::MatrixXd V, V2;
        Eigen::MatrixXi F, F2; // they will be allocated in detail_to_eigen()

        SOP_IGL::detail_to_eigen(*gdp, V, F);
        SOP_IGL::detail_to_eigen(*deform_gdp, V2, F2);

        Eigen::MatrixXd U;
        Eigen::VectorXi S(numPoints);
        Eigen::VectorXi b;
        Eigen::RowVector3d mid;
        U = V;
        GA_ROHandleF  pintoanimation_h(gdp->findAttribute(GA_ATTRIB_POINT, "pintoanimation"));
        if (!pintoanimation_h.isValid()) {
            addWarning(SOP_MESSAGE, "Can do anything without pintoanimation attribute.");
            return error();
        }

        bool pinned_vertex_exists = false;
        bool free_vertex_exists = false;
        GA_Offset ptoff;
        GA_FOR_ALL_PTOFF(gdp, ptoff)
        {
            const int pin = pintoanimation_h.get(ptoff);

            if (pin == -1)
                free_vertex_exists = true;
            if (pin == 1)
                pinned_vertex_exists = true;

            S[static_cast<int>(ptoff)] = pin;
        }

        if (!pinned_vertex_exists || !free_vertex_exists) {
             addWarning(SOP_MESSAGE, "At least one vertex has to be pinned and some should be free.");
            return error();
        }
        // vertices in selection
        igl::colon<int>(0, V.rows()-1, b);
        b.conservativeResize(std::stable_partition(b.data(), b.data()+b.size(),\
            [&](int i)->bool{return S(i)>=0;})-b.data());
        // Centroid
        mid = 0.5*(V.colwise().maxCoeff() + V.colwise().minCoeff());
        // Precomputation
        igl::ARAPData arap_data;
        arap_data.energy = static_cast<igl::ARAPEnergyType>(ENERGY(t));
        arap_data.max_iter = MAXITER(t);
        if (!igl::arap_precomputation(V, F, V.cols(), b, arap_data)) {
            addWarning(SOP_MESSAGE, "Arap precomputation failed.");
            return error();   
        }
        // Copy deformed coordinates. 
        Eigen::MatrixXd bc(b.size(),V.cols());
        for(int i = 0;i<b.size();i++) {
          bc.row(i) = V2.row(b(i));
        }
        // Solve

        try {
            igl::arap_solve(bc, arap_data, U);
        }
        catch (...)
        {
            addWarning(SOP_MESSAGE, "Can't solve arap...");
            return error();
        }

        {
            GA_Offset ptoff;
            GA_FOR_ALL_PTOFF(gdp, ptoff) { 
                const int pin = pintoanimation_h.get(ptoff);
                // if (pin == 1)
                //     continue;
                const GA_Index ptidx = gdp->pointIndex(ptoff);  
                UT_ASSERT((uint)ptidx < U.rows());
                const UT_Vector3 pos(U((uint)ptidx, 0), U((uint)ptidx, 1), U((uint)ptidx, 2));
                gdp->setPos3(ptoff, pos);
            }
        }
    } 

    else if(method == BIHARMONIC_COORDINATES)
    {
        Eigen::MatrixXd W;
        igl::ARAPData arap_data;
        SOP_IGL::Mesh low, high, scene;

        const GU_Detail *rest_gdp   = inputGeo(1);
        const GU_Detail *deform_gdp = inputGeo(2);

        if (!rest_gdp || !deform_gdp) {
            addWarning(SOP_MESSAGE, "Needs second and thrid input.");
            return error();
        }

        if (rest_gdp->getNumPoints() != deform_gdp->getNumPoints())
        {
            addError(SOP_ERR_MISMATCH_POINT, "Rest and deform geometry should match.");
            return error();
        }

        std::cout << "before detail_to_eigen() \n";
        SOP_IGL::detail_to_eigen(*rest_gdp, low.V, low.F, low.T);
        SOP_IGL::detail_to_eigen(*gdp, high.V, high.F, high.T);
        std::cout << "after detail_to_eigen() \n";
        // std::cout << low.V.rows() << "x" << low.V.cols() << "\n";
        // std::cout << low.F.rows() << "x" << low.F.cols() << "\n";
        // std::cout << low.T.rows() << "x" << low.T.cols() << "\n";

        // std::cout << high.V.rows() << "x" << high.V.cols() << "\n";
        // std::cout << high.F.rows() << "x" << high.F.cols() << "\n";
        // std::cout << high.T.rows() << "x" << high.T.cols() << "\n";
       
        // Precomputation
        {
            std::cout << "beforeEigen::VectorXi::LinSpaced() \n";
            Eigen::VectorXi b;
            {
                Eigen::VectorXi J = Eigen::VectorXi::LinSpaced(high.V.rows(),0,high.V.rows()-1);
                Eigen::VectorXd sqrD;
                Eigen::MatrixXd _2;
                std::cout<<"Finding closest points...\n";
                igl::point_mesh_squared_distance(low.V, high.V, J, sqrD, b, _2);
                assert(sqrD.minCoeff() < 1e-7 && "low.V should exist in high.V");
            }
            // force perfect positioning, rather have popping in low-res than high-res.
            // The correct/elaborate thing to do is express original low.V in terms of
            // linear interpolation (or extrapolation) via elements in (high.V,high.F)
            igl::slice(high.V, b, 1, low.V);
            // list of points --> list of singleton lists
            std::vector<std::vector<int> > S;
            igl::matrix_to_list(b, S);
            std::cout << "Computing weights for " << b.size() <<
            " handles at " << high.V.rows() << " vertices...\n";
            // Technically k should equal 3 for smooth interpolation in 3d, but 2 is
            // faster and looks OK
            const int k = 2;
            std::cout << " Before biharmonic_coordinates \n";
            std::cout << low.V.rows() << "x" << low.V.cols() << "\n";
            std::cout << low.T.rows() << "x" << low.T.cols() << "\n";
            std::cout << high.V.rows() << "x" << high.V.cols() << "\n";
            std::cout << high.T.rows() << "x" << high.T.cols() << "\n";
            std::cout << W.rows() << "x" << W.cols() << "\n";
            std::cout << S.size() << "\n";
            igl::biharmonic_coordinates(high.V,high.T,S,k,W);
            std::cout << W.rows() << "x" << W.cols() << "\n";
            // igl::biharmonic_coordinates(high.V, high.T, S, k, W);
            std::cout<<"Reindexing...\n";
            // Throw away interior tet-vertices, keep weights and indices of boundary
            Eigen::VectorXi I, J;
            igl::remove_unreferenced(high.V.rows(), high.F, I, J);
            std::for_each(high.F.data(), high.F.data() + high.F.size(), [&I](int & a){a=I(a);});
            std::for_each(b.data(), b.data() + b.size(), [&I](int & a){a=I(a);});
            igl::slice(Eigen::MatrixXd(high.V), J, 1, high.V);
            igl::slice(Eigen::MatrixXd(W), J, 1, W);
        }

        // Resize low res (high res will also be resized by affine precision of W)
        low.V.rowwise() -= low.V.colwise().mean();
        low.V /= (low.V.maxCoeff()-low.V.minCoeff());
        low.V.rowwise() += Eigen::RowVector3d(0,1,0);
        low.U = low.V;
        high.U = high.V;

        arap_data.max_iter = MAXITER(t);
        arap_data.energy = static_cast<igl::ARAPEnergyType>(ENERGY(t));
        if(!igl::arap_precomputation(low.V, low.T, 3, Eigen::VectorXi(), arap_data))
        {
            addWarning(SOP_MESSAGE, "Arap precomputation failed.");
            return error();
        }

       { 
            GA_Offset ptoff;
           GA_FOR_ALL_PTOFF(deform_gdp, ptoff) {
               const UT_Vector3 pos = deform_gdp->getPos3(ptoff);
               low.U(static_cast<uint>(ptoff), 0) = pos.x();
               low.U(static_cast<uint>(ptoff), 1) = pos.y(); 
               low.U(static_cast<uint>(ptoff), 2) = pos.z();
           }
        }


        arap_solve(Eigen::MatrixXd(0,3), arap_data, low.U);


        {
            GA_Offset ptoff;
            GA_FOR_ALL_PTOFF(gdp, ptoff) { 
                const GA_Index ptidx = gdp->pointIndex(ptoff);  
                UT_ASSERT((uint)ptidx < high.U.rows());
                const UT_Vector3 pos(high.U((uint)ptidx, 0), high.U((uint)ptidx, 1), high.U((uint)ptidx, 2));
                gdp->setPos3(ptoff, pos);
            }
        }



    }


    // 
    gdp->getP()->bumpDataId();
    return error();
}
