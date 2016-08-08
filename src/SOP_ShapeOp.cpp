#define SHAPEOP_HEADER_ONLY
#include <libShapeOp/src/Solver.h>
#include <libShapeOp/src/Constraint.h>


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
#include <set>

#include "SOP_ShapeOp.hpp"

using namespace SOP_SHAPEOP;


namespace SOP_SHAPEOP {

void getPointNeighbours(const GU_Detail *gdp, const GA_Offset ptoff, std::set<GA_Offset> &pointList)
{
    GA_OffsetArray pointVertices;
    gdp->getVerticesReferencingPoint(pointVertices, ptoff);
    for(int i=0; i < pointVertices.size(); ++i) 
    {
        const GA_Offset primoff   = gdp->vertexPrimitive(pointVertices(i));
        const GA_Primitive *prim  = gdp->getPrimitive(primoff);
        const GA_Range vertices   = prim->getPointRange();
        GA_Range::const_iterator it;
        for (it=vertices.begin(); !it.atEnd(); ++it) 
        {
            GA_Offset v1, v2;
            if (prim->findEdgePoints(ptoff, *it, v1, v2))
                pointList.insert(*it);
        }
    }
}

} // end of SOP_SHAPEOP space


static PRM_Name names[] = {
    PRM_Name("maxiter",   "Max Solver interations"),
    PRM_Name("closeness", "Closeness"),
    PRM_Name("edgestrain","EdgeStrain"),
    PRM_Name("plane",     "Plane"),
    PRM_Name("similarity", "Similarity"), 
    PRM_Name("laplacian",  "Laplacian"),
};

static PRM_Name  termChoices[] =
{
    PRM_Name("0", "Edge Strain"),
    PRM_Name("1", "Triangle Strain"),
    PRM_Name(0)
};


static PRM_ChoiceList  termMenu(PRM_CHOICELIST_SINGLE,  termChoices);
static PRM_Range       maxiterRange(PRM_RANGE_UI, 0, PRM_RANGE_UI, 100);
static PRM_Range       laplaceRange(PRM_RANGE_UI, -1, PRM_RANGE_UI, 1);

PRM_Template
SOP_ShapeOp::myTemplateList[] = {
    PRM_Template(PRM_INT_J, 1, &names[0], PRMzeroDefaults, 0, &maxiterRange),
    PRM_Template(PRM_FLT_J, 1, &names[1], PRMzeroDefaults),
    PRM_Template(PRM_FLT_J, 1, &names[2], PRMzeroDefaults),
    PRM_Template(PRM_FLT_J, 1, &names[3], PRMzeroDefaults),
    PRM_Template(PRM_FLT_J, 1, &names[4], PRMzeroDefaults),
    PRM_Template(PRM_FLT_J, 1, &names[5], PRMzeroDefaults, 0, &laplaceRange),
    PRM_Template(),
};


OP_Node *
SOP_ShapeOp::myConstructor(OP_Network *net, const char *name, OP_Operator *op)
{
    return new SOP_ShapeOp(net, name, op);
}

SOP_ShapeOp::SOP_ShapeOp(OP_Network *net, const char *name, OP_Operator *op)
    : SOP_Node(net, name, op), myGroup(NULL), mySolver(NULL)
{
   
    mySopFlags.setManagesDataIDs(true);
    mySolver = new ShapeOp::Solver(); //shapeop_create();
}

SOP_ShapeOp::~SOP_ShapeOp() 
{
     // shapeop_delete(mySolver);
     if (mySolver)
        delete mySolver;
}

OP_ERROR
SOP_ShapeOp::cookInputGroups(OP_Context &context, int alone)
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
SOP_ShapeOp::cookMySop(OP_Context &context)
{
    
    OP_AutoLockInputs inputs(this);
    if (inputs.lock(context) >= UT_ERROR_ABORT)
        return error();

    fpreal t = context.getTime();
    duplicateSource(0, context);

     // Get rest and deform geometry:
    const GU_Detail *rest_gdp   = inputGeo(1);
    const GU_Detail *shape_gdp  = inputGeo(2);

    // Point count should match:
    if (rest_gdp->getNumPoints() != gdp->getNumPoints())
    {
        addError(SOP_ERR_MISMATCH_POINT, "Rest and deform geometry should match.");
        return error();
    }

    // Copy to eigen.
    //gdp->convex(); // only triangles for now.
    uint numPoints = gdp->getNumPoints();
    uint numPrims  = gdp->getNumPrimitives();
    uint numConstraints = 0; // we need to book keeping this as solver crashes with #constr. == 0
    std::vector<ShapeOp::Matrix3X> shapes;
    

    const int maxiterations    = MAXITER(t);

    // Eearly quite on 0 iterations
    if (!maxiterations) {
        addWarning(SOP_MESSAGE, "Increase iterations to run a solver.");
        return error();
    }

    GA_Offset ptoff;
    // TODO: place for using GA_AIFNumericArray?
    // std::vector<double> pos_vector;
    ShapeOp::Matrix3X positions(3, numPoints);
    std::vector<int> indices; 
    GA_FOR_ALL_PTOFF(gdp, ptoff) 
    {
        const UT_Vector3 pos = gdp->getPos3(ptoff);
        const int pidx       = gdp->pointIndex(ptoff);
        positions(0, pidx) = (double)pos.x();
        positions(1, pidx) = (double)pos.y();
        positions(2, pidx) = (double)pos.z();
        indices.push_back(pidx); // TODO: Do we need this?

    }

    mySolver->setPoints(positions);

    const double closeness = CLOSENESS(t);
    if(closeness)
    {
        GA_FOR_ALL_PTOFF(gdp, ptoff) 
        {
            const int pidx = gdp->pointIndex(ptoff);
            const UT_Vector3 pos = gdp->getPos3(ptoff);
            std::vector<int> indices(1, pidx);
            std::shared_ptr<ShapeOp::ClosenessConstraint> \
            constraint(new ShapeOp::ClosenessConstraint(indices, closeness, positions));
            if(constraint) {
                ShapeOp::Vector3 v(pos.x(), pos.y(), pos.z());
                constraint->setPosition(v);
                mySolver->addConstraint(constraint);
                numConstraints++;
            } else {
                addWarning(SOP_MESSAGE, "Can't setup some constraints...");
            }
        }
    }
        
 
    const double edgestrain = EDGESTRAIN(t);
    if (edgestrain) 
    {

        UniqueEdges edges; // unique edges.
        GA_Offset ptoff;
        GA_FOR_ALL_PTOFF(rest_gdp, ptoff)
        {
            const int pidx = (int)rest_gdp->pointIndex(ptoff);
            std::set<int> uniques;
            edges.insert(std::pair<int, std::set<int> >(pidx, uniques));
            std::set<GA_Offset> neighbours;
            getPointNeighbours(rest_gdp, ptoff, neighbours);
            std::set<GA_Offset>::const_iterator it;
            for (it=neighbours.begin(); it!=neighbours.end(); ++it)
            {
                const int nidx = (int)rest_gdp->pointIndex(*it);
                // Ommit a reapting edges.
                UniqueEdges::iterator jt = edges.find(nidx); 
                if(jt != edges.end()) 
                {
                    if (jt->second.find(pidx) != jt->second.end()) 
                        continue; //?
                } else {
                    edges.find(pidx)->second.insert(nidx);
                    std::vector<int> indices(2);
                    indices[0] = pidx; indices[1] = nidx;

                    const UT_Vector3 vp1 = rest_gdp->getPos3(ptoff);
                    const UT_Vector3 vp2 = rest_gdp->getPos3(*it);
                    const double edge_length = (double)UT_Vector3(vp2 - vp1).length();

                    std::shared_ptr<ShapeOp::EdgeStrainConstraint> \
                    constraint(new ShapeOp::EdgeStrainConstraint(indices, edgestrain, positions));
                    if(constraint) {
                        constraint->setEdgeLength(edge_length);
                        mySolver->addConstraint(constraint);
                        numConstraints++;
                    } else {
                        addWarning(SOP_MESSAGE, "Can't setup some constraints...");
                    }
                }  
            }

        }

    }

    const double plane = PLANE(t);
    if (plane) {
        std::shared_ptr<ShapeOp::PlaneConstraint> \
        constraint(new ShapeOp::PlaneConstraint(indices, plane, positions));
        if(constraint) {
            mySolver->addConstraint(constraint);
            numConstraints++;
        }
    }

    const float similarity = SIMILARITY(t);

    if (similarity && !shape_gdp) {
        addWarning(SOP_MESSAGE, "No shape geometry connected to third input.");
    }

    if (similarity && shape_gdp) {
        const int shapeNumPoints = shape_gdp->getNumPoints();
        ShapeOp::Matrix3X shape(3, shapeNumPoints);
        std::vector<int> indices;
        GA_FOR_ALL_PTOFF(shape_gdp, ptoff) {
            const UT_Vector3 pos = shape_gdp->getPos3(ptoff);
            const int pidx       = shape_gdp->pointIndex(ptoff);
            shape(0, pidx) = (double)pos.x();
            shape(1, pidx) = (double)pos.y();
            shape(2, pidx) = (double)pos.z();
            indices.push_back(pidx);
        }
        shapes.push_back(shape); // 'shapes' declared above in cookMySop() space.
        std::shared_ptr<ShapeOp::SimilarityConstraint> \
        constraint(new ShapeOp::SimilarityConstraint(indices, similarity, positions));
        if(constraint) {
            constraint->setShapes(shapes);
            mySolver->addConstraint(constraint);
            numConstraints++;
        } else {
            addWarning(SOP_MESSAGE, "Can't setup some constraints...");
        }

    }


    const float laplacian = LAPLACIAN(t);
    if (laplacian) {
        const float abs_laplacian = SYSabs(laplacian);
        bool disp_lap = false;
        if (laplacian < 0.f)
            disp_lap = true;

        GA_Offset ptoff;
        GA_FOR_ALL_PTOFF(gdp, ptoff) {
            std::vector<int> indices;
            const int ptidx = gdp->pointIndex(ptoff);
            indices.push_back(ptidx);
            std::set<GA_Offset> neighbours;
            getPointNeighbours(gdp, ptoff, neighbours);
            std::set<GA_Offset>::const_iterator it;
            for (it=neighbours.begin(); it!=neighbours.end(); ++it) {
                const int nidx = (int)gdp->pointIndex(*it);
                indices.push_back(nidx);
            }

            std::shared_ptr<ShapeOp::UniformLaplacianConstraint> \
            constraint(new ShapeOp::UniformLaplacianConstraint(indices, abs_laplacian, positions, disp_lap));
            if(constraint) {
                mySolver->addConstraint(constraint);
                numConstraints++;
            } else {
                addWarning(SOP_MESSAGE, "Can't setup some constraints...");
            }
        }
    }



    // Solver will crash with numConstraints == 0
    if (!numConstraints) {
        addWarning(SOP_MESSAGE, "No constrained specified.");
        return error();
    }

     if(!mySolver->initialize()) {
            addWarning(SOP_MESSAGE, "Can't initialize solver.");
            return error();
        }

    if(!mySolver->solve(maxiterations)){
        addWarning(SOP_MESSAGE, "Can't solve.");
        return error();
    }

    UT_ASSERT(positions.rows() == numPoints);
    positions = mySolver->getPoints();

    {
        GA_Offset ptoff;
        GA_FOR_ALL_PTOFF(gdp, ptoff) {
            const uint ptidx = gdp->pointIndex(ptoff);
            const UT_Vector3 pos(positions(0, ptidx), positions(1, ptidx), positions(2, ptidx));
            gdp->setPos3(ptoff, pos);
        }
    }
    
    gdp->getP()->bumpDataId();
    return error();
}
