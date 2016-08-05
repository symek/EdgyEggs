#define SHAPEOP_HEADER_ONLY
#include <libShapeOp/api/API.cpp>


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
};

static PRM_Name  termChoices[] =
{
    PRM_Name("0", "Edge Strain"),
    PRM_Name("1", "Triangle Strain"),
    PRM_Name(0)
};


static PRM_Default     pointNine(0.9);
static PRM_ChoiceList  termMenu(PRM_CHOICELIST_SINGLE,  termChoices);
static PRM_Range       maxiterRange(PRM_RANGE_UI, 0, PRM_RANGE_UI, 100);

PRM_Template
SOP_ShapeOp::myTemplateList[] = {
    PRM_Template(PRM_INT_J, 1, &names[0], PRMzeroDefaults, 0, &maxiterRange),
    PRM_Template(PRM_FLT_J, 1, &names[1], &pointNine),
    PRM_Template(PRM_FLT_J, 1, &names[2], PRMzeroDefaults),
    PRM_Template(PRM_FLT_J, 1, &names[3], PRMzeroDefaults),
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
    mySolver = shapeop_create();
}

SOP_ShapeOp::~SOP_ShapeOp() 
{
     shapeop_delete(mySolver);
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
    

    const int maxiterations    = MAXITER(t);

    // Eearly quite on 0 iterations
    if (!maxiterations) {
        addWarning(SOP_MESSAGE, "Increase iterations to run a solver.");
        return error();
    }

    GA_Offset ptoff;
    // TODO: place for using GA_AIFNumericArray?
    std::vector<double> pos_vector;
    GA_FOR_ALL_PTOFF(gdp, ptoff) 
    {
        const UT_Vector3 pos = gdp->getPos3(ptoff);
        pos_vector.push_back((double)pos.x());
        pos_vector.push_back((double)pos.y());
        pos_vector.push_back((double)pos.z());

    }

    shapeop_setPoints(mySolver, static_cast<ShapeOpScalar*>(&pos_vector[0]), numPoints*3);


    const double closeness = CLOSENESS(t);
    if(closeness)
    {
        // Obligatory? 
        GA_Offset ptoff;
        GA_FOR_ALL_PTOFF(gdp, ptoff)
        {
            const int ptidx = gdp->pointIndex(ptoff);
            const UT_Vector3 pos = gdp->getPos3(ptoff);
            const double p[3]    = {pos.x(), pos.y(), pos.z()};
            const int const_id   = shapeop_addConstraint(mySolver, "Closeness", (int*)&ptidx, 1, closeness);
            if (const_id == -1)
                continue;
            if (SO_SUCCESS != shapeop_editConstraint(mySolver, "Closeness", const_id, p, 3)) {
                addWarning(SOP_MESSAGE, "Can't setup some constraints..."); // TODO: how to handle errors.
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
                        continue;
                } else {
                    edges.find(pidx)->second.insert(nidx);
                    const int indices[2] = {pidx, nidx};
                    const int constraint_id = shapeop_addConstraint(mySolver, "EdgeStrain", (int*)indices, 2,  edgestrain);
                    if (constraint_id == -1) {
                        addWarning(SOP_MESSAGE, "Some errors in constraints occured.");
                        continue; // TODO?
                    }
                    const UT_Vector3 vp1 = rest_gdp->getPos3(ptoff);
                    const UT_Vector3 vp2 = rest_gdp->getPos3(*it);
                    const double edge_length = (double)UT_Vector3(vp2 - vp1).length();
                    const double fraction = edge_length / 2.f;
                    const double parms[3] = {edge_length, edge_length-fraction, edge_length+fraction};
                    if (SO_SUCCESS != shapeop_editConstraint(mySolver, "EdgeStrain", constraint_id, parms, 3)) 
                    {
                        addWarning(SOP_MESSAGE, "Can't setup some constraints..."); // TODO: how to handle errors.
                    }
                    #if 0
                    std::cout << "Edge: " << indices[0] << "--" << indices[1] << ". L/Min/Max" << parms[0];
                    std::cout << "/" << parms[1] << "/" << parms[2] << "\n";
                    #endif
                }
            }

        }

    }


    const double plane = PLANE(t);
    if (plane)
    {
        GA_Offset ptoff;
        std::vector<int> indices;
        GA_FOR_ALL_PTOFF(gdp, ptoff) {
            const int pidx = gdp->pointIndex(ptoff);
            indices.push_back(pidx);
        }

        shapeop_addConstraint(mySolver, "Plane", (int*)&indices[0], numPoints, plane);

    }


     if(shapeop_init(mySolver)) {
            addWarning(SOP_MESSAGE, "Can't initialize solver.");
            return error();
        }

    if(shapeop_solve(mySolver, maxiterations)){
        addWarning(SOP_MESSAGE, "Can't solve.");
        return error();
    }

    UT_ASSERT(pos_vector.size() == numPoints*3);
    shapeop_getPoints(mySolver, &pos_vector[0], numPoints*3);

    {
        GA_Offset ptoff;
        GA_FOR_ALL_PTOFF(gdp, ptoff) {
            const uint ptidx = gdp->pointIndex(ptoff);
            const UT_Vector3 pos(pos_vector[3*ptidx], pos_vector[1+3*ptidx], pos_vector[2+3*ptidx]);
            gdp->setPos3(ptoff, pos);
        }
    }
    
    gdp->getP()->bumpDataId();
    return error();
}
