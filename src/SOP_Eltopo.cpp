
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

#include "SOP_Eltopo.hpp"


using namespace SOP_ELTOPO;


namespace SOP_ELTOPO {

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

} // end of SOP_ELTOPO space


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
SOP_Eltopo::myTemplateList[] = {
    // PRM_Template(PRM_INT_J, 1, &names[0], PRMzeroDefaults, 0, &maxiterRange),
    // PRM_Template(PRM_FLT_J, 1, &names[1], PRMzeroDefaults),
    // PRM_Template(PRM_FLT_J, 1, &names[2], PRMzeroDefaults),
    // PRM_Template(PRM_FLT_J, 1, &names[3], PRMzeroDefaults),
    // PRM_Template(PRM_FLT_J, 1, &names[4], PRMzeroDefaults),
    // PRM_Template(PRM_FLT_J, 1, &names[5], PRMzeroDefaults, 0, &laplaceRange),
    PRM_Template(),
};


OP_Node *
SOP_Eltopo::myConstructor(OP_Network *net, const char *name, OP_Operator *op)
{
    return new SOP_Eltopo(net, name, op);
}

SOP_Eltopo::SOP_Eltopo(OP_Network *net, const char *name, OP_Operator *op)
    : SOP_Node(net, name, op), myGroup(NULL)
{
   
    mySopFlags.setManagesDataIDs(true);
}

SOP_Eltopo::~SOP_Eltopo() 
{
  
}



OP_ERROR
SOP_Eltopo::cookMySop(OP_Context &context)
{
    
    OP_AutoLockInputs inputs(this);
    if (inputs.lock(context) >= UT_ERROR_ABORT)
        return error();

    fpreal t = context.getTime();
    duplicateSource(0, context);

    // Only triangles
    gdp->convex();
    const int numPoints = gdp->getNumPoints();
    const int numPrims  = gdp->getNumPrimitives();

    std::vector<Vec3d> positions;
    std::vector<double> masses;
    GA_Offset ptoff;
    GA_FOR_ALL_PTOFF(gdp, ptoff) 
    {
        const UT_Vector3 pos = gdp->getPos3(ptoff);
        positions.push_back(Vec3d(pos.x(), pos.y(), pos.z()));
        masses.push_back(1.f);
    }

    std::vector<Vec3st> faces;
    GA_Iterator it(gdp->getPrimitiveRange());
    for (; !it.atEnd(); ++it)
    {
        const GEO_Primitive *prim = gdp->getGEOPrimitive(*it);
        GA_Primitive::const_iterator vt;
        Vec3st v;
        int vertex_index = 0;
        for (prim->beginVertex(vt); !vt.atEnd(); ++vt) {
            const GA_Offset voff = vt.getPointOffset();
            const int ptidx      = gdp->pointIndex(voff);
            v[SYSmin(vertex_index,2)] = static_cast<int>(ptidx);
            vertex_index++;
        }
        faces.push_back(v);
    }

    SurfTrackInitializationParameters p;
    p.m_max_volume_change = 1.f;   
    p.m_min_edge_length   = .01f;
    p.m_max_edge_length   = 1.f;
    p.m_collision_safety = false;
    SurfTrack surface_tracker(positions, faces, masses, p);
    // surface_tracker = std::make_shared<SurfTrack>(positions, faces, masses, p;
    surface_tracker.improve_mesh();
    
    // // do merging
    surface_tracker.topology_changes();
    
    surface_tracker.defrag_mesh();

    const int num_new_points = surface_tracker.get_num_vertices();
    const int num_new_prims  = surface_tracker.m_mesh.num_triangles();

    gdp->clearAndDestroy();
    gdp->appendPointBlock(num_new_points);

   { 
        GA_Offset ptoff;
        GA_FOR_ALL_PTOFF(gdp, ptoff) 
           {
              const Vec3d& v = surface_tracker.get_position((int)ptoff);
              const UT_Vector3 pos(v[0], v[1], v[2]);
              gdp->setPos3(ptoff, pos);
           }
   }

    for (int i=0; i< num_new_prims; ++i) 
    {
        const Vec3st& curr_tri = surface_tracker.m_mesh.get_triangle(i); 
        GU_PrimPoly *prim = GU_PrimPoly::build(gdp, 0, false, false);
        prim->appendVertex(curr_tri[0]);
        prim->appendVertex(curr_tri[1]);
        prim->appendVertex(curr_tri[2]); 
    }


    // gdp->getP()->bumpDataId();
    return error();
}
  