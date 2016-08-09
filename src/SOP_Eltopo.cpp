
#include <GU/GU_Detail.h>
#include <GU/GU_PrimPoly.h>
#include <OP/OP_Operator.h>
#include <OP/OP_Director.h>
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
    PRM_Name("collisionsafety",      "Assume non-intersecting input"),
    PRM_Name("allowtopologychange",  "Allow Topology Change"),
    PRM_Name("allownonnanifold",     "Allow Non Manifold"),
    PRM_Name("performimprovment",    "Perform Improvment"),
    PRM_Name("allowvertexmovement",  "Allow Vertex Movement"),
    PRM_Name("usefraction",          "Use Fraction"), 
    PRM_Name("usecurvaturewhenspliting",  "Use Curvature (spliting)"),
    PRM_Name("usecurvaturewhencollapsing","Use Curvature (collapsing)"),

    PRM_Name("minedgelength",        "Min Edge Length"),
    PRM_Name("maxedgelength",        "Max Edge Length"),
    PRM_Name("mintrianglearea",      "Min Triangle Area"),
    PRM_Name("maxvolumechange",      "Max Volume Change"),
    PRM_Name("edgeflipminlength",    "Edge Flip Min Length Change"),
    PRM_Name("mergeproximityeps",    "Merge Proximity Epsilon"),
    PRM_Name("subdivisionscheme",    "Subdivision Scheme"),
      
};



static PRM_Name  subdivideSchemeChoices[] =
{
    PRM_Name("0", "MidPoint"),
    PRM_Name("1", "Butterfly"),
    PRM_Name("2", "Modified Butterfly"),
    PRM_Name("3", "Quadratic Error Min"),
    PRM_Name(0)
};

static PRM_ChoiceList  subdivideSchemeMenu(PRM_CHOICELIST_SINGLE,  subdivideSchemeChoices);

PRM_Template
SOP_Eltopo::myTemplateList[] = {
    PRM_Template(PRM_TOGGLE, 1, &names[0], PRMzeroDefaults),
    PRM_Template(PRM_TOGGLE, 1, &names[1], PRMzeroDefaults),
    PRM_Template(PRM_TOGGLE, 1, &names[2], PRMzeroDefaults),
    PRM_Template(PRM_TOGGLE, 1, &names[3], PRMzeroDefaults),
    PRM_Template(PRM_TOGGLE, 1, &names[4], PRMzeroDefaults),
    PRM_Template(PRM_TOGGLE, 1, &names[5], PRMzeroDefaults),
    PRM_Template(PRM_TOGGLE, 1, &names[6], PRMzeroDefaults),
    PRM_Template(PRM_TOGGLE, 1, &names[7], PRMzeroDefaults),

    PRM_Template(PRM_FLT_J, 1, &names[8],  PRMpointOneDefaults),
    PRM_Template(PRM_FLT_J, 1, &names[9],  PRMoneDefaults),
    PRM_Template(PRM_FLT_J, 1, &names[10], PRMpointOneDefaults),
    PRM_Template(PRM_FLT_J, 1, &names[11], PRMoneDefaults),
    PRM_Template(PRM_FLT_J, 1, &names[12], PRMpointOneDefaults),
    PRM_Template(PRM_FLT_J, 1, &names[13], PRMpointOneDefaults),
    PRM_Template(PRM_ORD,   1, &names[14], 0, &subdivideSchemeMenu, 0, 0),

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
    const double delta = 1.f / OPgetDirector()->getChannelManager()->getSamplesPerSec();
    duplicateSource(0, context);
    const GU_Detail *advected_gdp = inputGeo(1);

    // Only triangles
    gdp->convex();
    const int numPoints = gdp->getNumPoints();
    const int numPrims  = gdp->getNumPrimitives();


    if(advected_gdp) {
        if (advected_gdp->getNumPoints() != numPoints)
        {
            addWarning(SOP_MESSAGE, "Advected geometry doesn't match.");
            advected_gdp = NULL;
        }
    }

    std::vector<Vec3d> positions;
    std::vector<double> masses;
    GA_Offset ptoff;
    GA_FOR_ALL_PTOFF(gdp, ptoff) 
    {
        const UT_Vector3 pos = gdp->getPos3(ptoff);
        positions.push_back(Vec3d(pos.x(), pos.y(), pos.z()));
        masses.push_back(0.1f);
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


    SurfTrackInitializationParameters parms;
    parms.m_collision_safety       = COLLISIONSAFETY(t);
    parms.m_allow_topology_changes = ALLOWTOPOLOGYCHANGE(t);
    parms.m_allow_non_manifold     = ALLOWNONNANIFOLD(t);
    parms.m_perform_improvement    = PERFORMIMPROVMENT(t);
    parms.m_allow_vertex_movement  = ALLOWVERTEXMOVEMENT(t);
    parms.m_use_fraction           = USEFRACTION(t);
    parms.m_use_curvature_when_splitting  = USECURVATUREWHENSPLITING(t);
    parms.m_use_curvature_when_collapsing = USECURVATUREWHENCOLLAPSING(t);
    parms.m_min_edge_length   = SYSmax(MINEDGELENGTH(t), 1e-7);
    parms.m_max_edge_length   = MAXEDGELENGTH(t);
    parms.m_min_triangle_area = SYSmax(MINTRIANGLEAREA(t),  1e-7);
    parms.m_max_volume_change = MAXVOLUMECHANGE(t);  
    parms.m_edge_flip_min_length_change = SYSmax(EDGEFLIPMINLENGTH(t), 1e-7);
    parms.m_merge_proximity_epsilon     = SYSmax(MERGEPROXIMITYEPS(t), 1e-7);

    const int subdivisionscheme = SUBDIVISIONSCHEME(t);
    std::shared_ptr<SubdivisionScheme> scheme;

    switch(subdivisionscheme)
    {
        case MIDPOINT:
            scheme = std::make_shared<MidpointScheme>(); break;
        case BUTTERFLY:
            scheme = std::make_shared<ButterflyScheme>(); break;
        case MODIFIED_BUTTERFLY:
            scheme = std::make_shared<ModifiedButterflyScheme>(); break;
        case QUADRATIC_ERROR_MIN:
            scheme = std::make_shared<QuadraticErrorMinScheme>(); break;
        default:
            scheme = std::make_shared<MidpointScheme>(); break;
    }


    parms.m_subdivision_scheme = scheme.get();
    SurfTrack surface_tracker(positions, faces, masses, parms);

    if (advected_gdp) {
        std::vector<Vec3d> new_positions;
        GA_Offset ptoff;
        GA_FOR_ALL_PTOFF(advected_gdp, ptoff) {
            const UT_Vector3 pos = advected_gdp->getPos3(ptoff);
            new_positions.push_back(Vec3d(pos.x(), pos.y(), pos.z()));
        }

        surface_tracker.set_all_newpositions(new_positions);
        double actual_delta = 0.f;
        surface_tracker.integrate(delta, actual_delta);
    }

   
    // Improve mesh.
    surface_tracker.improve_mesh();
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

    return error();
}
  