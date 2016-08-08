
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
    my_mesh         = std::make_shared<ElTopoMesh>();
    static_options  = std::make_shared<ElTopoStaticOperationsOptions>();
    general_options = std::make_shared<ElTopoGeneralOptions>();
}

SOP_Eltopo::~SOP_Eltopo() 
{
     // shapeop_delete(mySolver);
     // if (mySolver)
        // delete mySolver;
    // delete my_mesh;
    // delete static_options;
    // delete general_options;
}

OP_ERROR
SOP_Eltopo::cookInputGroups(OP_Context &context, int alone)
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

    std::vector<double> positions;
    std::vector<double> masses;
    GA_Offset ptoff;
    GA_FOR_ALL_PTOFF(gdp, ptoff) 
    {
        const UT_Vector3 pos = gdp->getPos3(ptoff);
        positions.push_back((double)pos.x());
        positions.push_back((double)pos.y());
        positions.push_back((double)pos.z());
        masses.push_back(1.f);
    }

    std::vector<int> faces;
    GA_Iterator it(gdp->getPrimitiveRange());
    for (; !it.atEnd(); ++it)
    {
        const GEO_Primitive *prim = gdp->getGEOPrimitive(*it);
        GA_Primitive::const_iterator vt;
        for (prim->beginVertex(vt); !vt.atEnd(); ++vt) {
            const GA_Offset voff = vt.getPointOffset();
            const int ptidx      = gdp->pointIndex(voff);
            faces.push_back(static_cast<int>(ptidx));
        }
    }

    my_mesh->num_vertices     = numPoints;
    my_mesh->vertex_locations = (double*)&positions[0];
    my_mesh->num_triangles    = numPrims;
    my_mesh->triangles        = (int*)&faces[0];
    my_mesh->vertex_masses    = (double*)&masses[0];


     // whether to perform mesh maintenance

    // m_proximity_epsilon( 1e-4 ),
    // m_friction_coefficient( 0.0 ),
    // m_min_triangle_area( 1e-7 ),
    // m_improve_collision_epsilon( 2e-6 ),
    // m_use_fraction( false ),
    // m_min_edge_length( UNINITIALIZED_DOUBLE ),     // <- Don't allow instantiation without setting these parameters
    // m_max_edge_length( UNINITIALIZED_DOUBLE ),     // <-
    // m_max_volume_change( UNINITIALIZED_DOUBLE ),   // <-
    // m_min_triangle_angle( 0.0 ),
    // m_max_triangle_angle( 180.0 ),
    // m_use_curvature_when_splitting( false ),
    // m_use_curvature_when_collapsing( false ),
    // m_min_curvature_multiplier( 1.0 ),
    // m_max_curvature_multiplier( 1.0 ),
    // m_allow_vertex_movement( true ),
    // m_edge_flip_min_length_change( 1e-8 ),
    // m_merge_proximity_epsilon( 1e-5 ),
    // m_subdivision_scheme(NULL),
    // m_collision_safety(true),
    // m_allow_topology_changes(true),
    // m_allow_non_manifold(true),
    // m_perform_improvement(true)     
    // // whether to allow merging and separation
    // static_options->m_perform_improvement = true;            
    // static_options->m_allow_topology_changes = true;
    // // maximum allowable change in volume when performing mesh maintenance
    // static_options->m_max_volume_change = 1.f;       
    // // edges shorter than this length will be collapsed
    // static_options->m_min_edge_length = 0.01f;              
    // // edges longer then this length will be subdivided
    // static_options->m_max_edge_length = 1.f;              
    // static_options->m_min_triangle_area = 0.001f;
    // static_options->m_min_triangle_angle = 0.01f;
    // static_options->m_max_triangle_angle = 120.f;   
    // static_options->m_use_curvature_when_splitting = false;
    // static_options->m_use_curvature_when_collapsing = false;

    // // Clamp curvature scaling to these values
    // static_options->m_min_curvature_multiplier = 1.f;
    // static_options->m_max_curvature_multiplier = 1.f;
    // static_options->m_allow_vertex_movement = true;
    // / Minimum edge length improvement in order to flip an edge
    // static_options->m_edge_flip_min_length_change = 0.001f;
    // /// Elements within this distance will trigger a merge attempt   
    // static_options->m_merge_proximity_epsilon = 0.01f;
    // /// Type of subdivision to use when collapsing or splitting (butterfly, quadric error minimization, etc.)
    // static_options->m_subdivision_scheme = NULL;
    // /// Whether to enforce collision-free surfaces (including during mesh maintenance operations)
    //  static_options->m_collision_safety;
    // /// Wether to allow non-manifold (edges incident on more than two triangles)
    static_options->m_allow_non_manifold  = false;
    static_options->m_allow_topology_changes = false;
    // static_options->m_collision_safety = true;

    general_options->m_verbose = 1;
    // general_options->m_proximity_epsilon =  1e-5 ;

    ElTopoDefragInformation defrag_info;
    ElTopoMesh outputs;


    el_topo_static_operations(my_mesh.get(), general_options.get(), \
        static_options.get(), &defrag_info, &outputs);

    std::cout << defrag_info.num_vertex_changes << ", " << defrag_info.num_triangle_changes << "\n";

    // el_topo_free_static_operations_results(&outputs, &defrag_info);

    gdp->getP()->bumpDataId();
    return error();
}
  