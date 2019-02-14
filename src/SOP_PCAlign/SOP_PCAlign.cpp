#include <UT/UT_DSOVersion.h>
#include <GU/GU_Detail.h>
#include <OP/OP_Operator.h>
#include <OP/OP_OperatorTable.h>
#include <PRM/PRM_Include.h>
#include <OP/OP_AutoLockInputs.h>
#include <SYS/SYS_Math.h>

#include <ICP.h>

#ifdef BUILD_WITH_INTEL_FGR
#include "pcl_fpfh.hpp"
#include "intel_fgr.hpp"
#endif

#include "SOP_PCAlign.hpp"


#include <unordered_map>
#include <memory>


using namespace pcalign;
typedef double Scalar;
typedef Eigen::Matrix<Scalar, 3, Eigen::Dynamic> Vertices;

void
newSopOperator(OP_OperatorTable *table)
{
    table->addOperator(new OP_Operator(
        "pcalign",
        "EE PointCloudAlign",
        SOP_PCAlign::myConstructor,
        SOP_PCAlign::myTemplateList,
        2,
        2,
        0));
}

static PRM_Name names[] = {
    PRM_Name("usepenalty",   "Use penalty"),
    PRM_Name("p",            "P norm"),
    PRM_Name("mu",           "Penelty weight"),
    PRM_Name("alpha",        "Penalty factor"),
    PRM_Name("maxmu",        "Max penalty"),
    PRM_Name("maxicp",       "Max ICP iteration"),
    PRM_Name("maxouter",     "Max outer iteration"),
    PRM_Name("maxinner",     "Max inner iteration"),
    PRM_Name("stop",         "Stopping criteria"),
    PRM_Name("method",       "Method"),
    PRM_Name("weightfunc",   "Reweight function"),
};

static PRM_Name  modelChoices[] =
{
    PRM_Name("0", "Rigid Motion"),
    PRM_Name("1", "Sparse ICP"),
    PRM_Name("2", "Reweighted ICP"),
    PRM_Name("3", "Intel FGR"),
    PRM_Name(0)
};

static PRM_Name weightFuncChoices[] = 
{
    PRM_Name("0", "PNorm"),
    PRM_Name("1", "Tukey"),
    PRM_Name("2", "Fair"),
    PRM_Name("3", "Logistic"),
    PRM_Name("4", "Trimmed"),
    PRM_Name("5", "None"),
    PRM_Name(0)
};

static PRM_ChoiceList  methodMenu(PRM_CHOICELIST_SINGLE, modelChoices);
static PRM_ChoiceList  weightMenu(PRM_CHOICELIST_SINGLE, weightFuncChoices);

static PRM_Default alphaDefault(1.2);
static PRM_Default maxmuDefault(1e5);
static PRM_Default maxicpDefault(10);
static PRM_Default stopDefault(1e-5);


PRM_Template
SOP_PCAlign::myTemplateList[] = {
    PRM_Template(PRM_ORD,   1, &names[9], 0, &methodMenu),
    PRM_Template(PRM_TOGGLE,1, &names[0], PRMzeroDefaults),
    PRM_Template(PRM_FLT_J, 1, &names[1], PRMoneDefaults),
    PRM_Template(PRM_FLT_J, 1, &names[2], PRMtenDefaults),
    PRM_Template(PRM_FLT_J, 1, &names[3], &alphaDefault),
    PRM_Template(PRM_FLT_J, 1, &names[4], &maxmuDefault),
    PRM_Template(PRM_INT_J, 1, &names[5], &maxicpDefault),
    PRM_Template(PRM_INT_J, 1, &names[6], &maxicpDefault),
    PRM_Template(PRM_INT_J, 1, &names[7], PRMoneDefaults),
    PRM_Template(PRM_FLT_LOG, 1, &names[8],     &stopDefault),
    PRM_Template(PRM_ORD,     1, &names[10], PRMfiveDefaults, &weightMenu),
    PRM_Template(),
};

bool
SOP_PCAlign::updateParmsFlags()
{
    bool    changed = SOP_Node::updateParmsFlags();
    bool    method   = evalInt("method", 0, 0);
    // if (method == ALIGN_METHOD::RIGID){
        // changed |= enableParm("maxicp", 0);
    // }

    // changed |= enableParm("copcolor", use_path);
    return changed;
}

OP_Node *
SOP_PCAlign::myConstructor(OP_Network *net, const char *name, OP_Operator *op)
{
    return new SOP_PCAlign(net, name, op);
}

SOP_PCAlign::SOP_PCAlign(OP_Network *net, const char *name, OP_Operator *op)
    : SOP_Node(net, name, op)
{
    mySopFlags.setManagesDataIDs(true);  
}

SOP_PCAlign::~SOP_PCAlign() {}

OP_ERROR
SOP_PCAlign::cookMySop(OP_Context &context)
{
    flags().timeDep = 1;
    OP_AutoLockInputs inputs(this);
    if (inputs.lock(context) >= UT_ERROR_ABORT)
        return error();

    fpreal t = context.getTime();
    duplicatePointSource(0, context);

    // Get second geometry:
    const GU_Detail * target_gdp = inputGeo(1);

    // any points?:
    if (target_gdp->getNumPoints() == 0 || gdp->getNumPoints() == 0) {
        addError(SOP_MESSAGE, "Needs two points clouds to align.");
        return error();
    }

    const int   method  = METHOD();
    const int   penalty = USE_PENALTY(t);
    const float p_norm  = P_NORM(t);   
    const float mu      = MU(t);         
    const float alpha   = ALPHA(t);         
    const float max_mu  = MAX_MU(t);         
    const int   max_icp = MAX_ICP(t);       
    const int   max_o   = MAX_OUTER(t);   
    const int   max_i   = MAX_INNER(t);   
    const float stop    = STOP(t);  
    const int   weightfunc = WEIGHTFUNC(t);    
   
    Vertices target, source;
    source.resize(Eigen::NoChange, gdp->getNumPoints());
    target.resize(Eigen::NoChange, target_gdp->getNumPoints());

    GA_Offset ptoff;
    GA_FOR_ALL_PTOFF(gdp, ptoff){
        const UT_Vector3 pos = gdp->getPos3(ptoff);
        const GA_Index   idx = gdp->pointIndex(ptoff);
        target(0, idx) = pos.x();
        target(1, idx) = pos.y();
        target(2, idx) = pos.z();
    }

    {
        GA_FOR_ALL_PTOFF(target_gdp, ptoff){
            const UT_Vector3 pos = target_gdp->getPos3(ptoff);
            const GA_Index   idx = target_gdp->pointIndex(ptoff);
            source(0, idx) = pos.x();
            source(1, idx) = pos.y();
            source(2, idx) = pos.z();
        }
    }

    if (method == ALIGN_METHOD::RIGID) {
        Eigen::Affine3d t = RigidMotionEstimator::point_to_point(source, target);
        const UT_Matrix4F m(t(0,0), t(0,1), t(0,2), t(0,3),
                            t(1,0), t(1,1), t(1,2), t(1,3),
                            t(2,0), t(2,1), t(2,2), t(2,3),
                            t(3,0), t(3,1), t(3,2), t(3,3));

        GA_RWHandleM4  xform_h(gdp->addFloatTuple(GA_ATTRIB_DETAIL, "transform", 16));
        xform_h.set(GA_Offset(0), m);

        GA_FOR_ALL_PTOFF(gdp, ptoff) {
            UT_Vector3 pos = gdp->getPos3(ptoff);
            pos *= m;
            gdp->setPos3(ptoff, pos);
        }

        gdp->getP()->bumpDataId();
        return error();

    } else if (method == ALIGN_METHOD::SPARSE_ICP) {
        SICP::Parameters parms;
        parms.use_penalty = static_cast<bool>(penalty);
        parms.p = p_norm;
        parms.mu = mu;
        parms.alpha = alpha;
        parms.max_mu = max_mu;
        parms.max_icp = max_icp;
        parms.max_outer = max_o;
        parms.max_inner = max_i;
        parms.stop = stop;
        SICP::point_to_point(source, target, parms);
        GA_FOR_ALL_PTOFF(gdp, ptoff){
            const GA_Index i = gdp->pointIndex(ptoff);
            const UT_Vector3 pos(source(0, i), source(1, i), source(2, i));
            gdp->setPos3(ptoff, pos);
        }

    } else if (method == ALIGN_METHOD::REWEIGHTED_ICP) {
        ICP::Parameters parms;
        // this crashes... atm
        // parms.f = static_cast<ICP::Function>(weightfunc);
        parms.p = p_norm;
        parms.max_icp = max_icp;
        parms.max_outer = max_o;
        parms.stop = stop;
        ICP::point_to_point(source, target, parms);
        GA_FOR_ALL_PTOFF(gdp, ptoff){
            const GA_Index i = gdp->pointIndex(ptoff);
            const UT_Vector3 pos(source(0, i), source(1, i), source(2, i));
            gdp->setPos3(ptoff, pos);
        }

    } 
    #ifdef BUILD_WITH_INTEL_FGR
    else if (method == ALIGN_METHOD::INTEL_FGR) {
        sop_pcl::Points::Ptr       points (new sop_pcl::Points);
        sop_pcl::Normals::Ptr      normals(new sop_pcl::Normals);
        sop_pcl::FeatureHists::Ptr fpfhs(new sop_pcl::FeatureHists);

        sop_pcl::gdp_to_pcl(gdp, points);
        sop_pcl::estimate_normals(points, normals);
        sop_pcl::compute_fpfh(points, normals, fpfhs);

        std::vector<Eigen::Vector3f> positions(points->size());
        std::vector<Eigen::VectorXf> features(fpfhs->size());
        // both are vectors of eigen Vectors, should be easy to
        // cheat and swap buffers...    
        for (int i=0; i<points->size(); ++i) {
            const pcl::PointXYZ & point = points->points[i];
            const pcl::FPFHSignature33 & feat  = fpfhs->points[i];
            const Eigen::Vector3f epoint(point.x, point.y, point.z);
            Eigen::VectorXf efeat; efeat.resize(33);
            for (int j=0; j<33; ++j) {
                efeat(j) = feat.histogram[j];
            }
            positions[i] = epoint;
            features[i]  = efeat;
        }

        FastGlobalRegistration ifgr; 
        ifgr.pointcloud_.push_back(positions);
        ifgr.features_.push_back(features);
        // ifgr.NormalizePoints();
        // ifgr.AdvancedMatching();
        // ifgr.OptimizePairwise(true, ITERATION_NUMBER); 
        // Eigen::Matrix4f t = ifgr.GetTrans();

        // const UT_Matrix4F m(t(0,0), t(0,1), t(0,2), t(0,3),
        //                     t(1,0), t(1,1), t(1,2), t(1,3),
        //                     t(2,0), t(2,1), t(2,2), t(2,3),
        //                     t(3,0), t(3,1), t(3,2), t(3,3));

        // GA_RWHandleM4  xform_h(gdp->addFloatTuple(GA_ATTRIB_DETAIL, "transform", 16));
        // xform_h.set(GA_Offset(0), m);

        // {
        //     GA_FOR_ALL_PTOFF(gdp, ptoff) {
        //             UT_Vector3 pos = gdp->getPos3(ptoff);
        //             pos *= m;
        //             gdp->setPos3(ptoff, pos);
        //         }
        // }

        // GA_RWHandleV3 norm_h(gdp->addFloatTuple(GA_ATTRIB_POINT, "N", 3));
        // GA_FOR_ALL_PTOFF(gdp, ptoff) {
        //     const GA_Index i = gdp->pointIndex(ptoff);
        //     const pcl::Normal & n = normals->at(i);
        //     const UT_Vector3 normal(n.normal_x, n.normal_y, n.normal_z);
        //     norm_h.set(ptoff, normal);
        // }
    

        // gdp->getP()->bumpDataId();
        // return error();
    }
    #endif


    gdp->getP()->bumpDataId();
    return error();
}
