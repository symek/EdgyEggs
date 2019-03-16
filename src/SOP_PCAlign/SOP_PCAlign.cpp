#include <UT/UT_DSOVersion.h>
#include <GU/GU_Detail.h>
#include <OP/OP_Operator.h>
#include <OP/OP_OperatorTable.h>
#include <PRM/PRM_Include.h>
#include <OP/OP_AutoLockInputs.h>
#include <SYS/SYS_Math.h>

#include <ICP.h>
// TODO: move it to Cmake
//#define BUILD_WITH_CPD 1
#ifdef BUILD_WITH_CPD
#include <cpd/gauss_transform_fgt.hpp>
#include <cpd/gauss_transform.hpp>
#include <cpd/rigid.hpp>
#include <cpd/nonrigid.hpp>
#endif


#ifdef BUILD_WITH_INTEL_FGR
#include "pcl_fpfh.hpp"
#include "intel_fgr.hpp"
#endif

#include "SOP_PCAlign.hpp"

#include <unordered_map>
#include <memory>


using namespace pcalign;


void
newSopOperator(OP_OperatorTable *table)
{
    table->addOperator(new OP_Operator(
        "pcalign",
        "EE Point Cloud Align",
        SOP_PCAlign::myConstructor,
        SOP_PCAlign::myTemplateList,
        2,
        2,
        0));
}

static PRM_Name names[] = {
    PRM_Name("alignmethod",   "Align Method"),
    PRM_Name("maxiterations", "Max iterations"),
    PRM_Name("maxouteriter",  "Max outer iteration"),
    PRM_Name("maxinneriter",  "Max inner iteration"),
    PRM_Name("stopcritera",   "Stopping criteria"),
    PRM_Name("usepenalty",    "Use penalty"),
    PRM_Name("pnorm",         "P norm"),
    PRM_Name("penaltyweight", "Penelty weight"),
    PRM_Name("penaltyfactor", "Penalty factor"),
    PRM_Name("maxpenalty",    "Max penalty"),
    PRM_Name("reweightfunc",  "Reweight function"),
    PRM_Name("doscaling",     "Allow Scaling"),
    PRM_Name("doreflections", "Allow reflections"),
    PRM_Name("docorrespondance", "Compute correspondence"),
    PRM_Name("outliers",         "Outliers factor"),
    PRM_Name("dotransform",       "Apply rigid transform"),
};

static PRM_Name  modelChoices[] =
{
    PRM_Name("0", "Rigid Motion"),
    PRM_Name("1", "Sparse ICP"),
    PRM_Name("2", "Reweighted ICP"),
    PRM_Name("3", "Intel FGR"),
    PRM_Name("4", "Coherent Point Drift Rigid"),
    PRM_Name("5", "Coherent Point Drift Non Rigid"),
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

static PRM_Default maxIterDefault(10);
static PRM_Default alphaDefault(1.2);
static PRM_Default maxmuDefault(1e5);
static PRM_Default stopDefault(1e-5);
static PRM_Default outliersDefault(0.1);


PRM_Template
SOP_PCAlign::myTemplateList[] = {
    PRM_Template(PRM_ORD,   1, &names[0], 0, &methodMenu), // align method
    PRM_Template(PRM_INT_J, 1, &names[1], &maxIterDefault), // max iterations (rigid, sparse, reweighted sparse)
    PRM_Template(PRM_INT_J, 1, &names[2], &maxIterDefault), // max outer iterations
    PRM_Template(PRM_INT_J, 1, &names[3], PRMoneDefaults),  // max inner iterations
    PRM_Template(PRM_FLT_LOG, 1, &names[4], &stopDefault),    // stop criteria
    PRM_Template(PRM_TOGGLE,1, &names[5], PRMzeroDefaults), // use penalty
    PRM_Template(PRM_FLT_J, 1, &names[6], PRMoneDefaults), // P norm 
    PRM_Template(PRM_FLT_J, 1, &names[7], PRMtenDefaults), // penalty weight (mu)
    PRM_Template(PRM_FLT_J, 1, &names[8], &alphaDefault), // penalty factor (alpha)
    PRM_Template(PRM_FLT_J, 1, &names[9], &maxmuDefault), // max penalty weight (max mu)
    PRM_Template(PRM_ORD,   1, &names[10], PRMfiveDefaults, &weightMenu), // weight function
    PRM_Template(PRM_TOGGLE,1, &names[11], PRMoneDefaults), // allow scaling rigid CPD
    PRM_Template(PRM_TOGGLE,1, &names[12], PRMzeroDefaults), // allow reflections rigid CPD
    PRM_Template(PRM_TOGGLE,1, &names[13], PRMzeroDefaults), // compute correspondance rigid CPD
    PRM_Template(PRM_FLT_LOG, 1,&names[14], &outliersDefault),    // outliers CPD
    PRM_Template(PRM_TOGGLE,1,  &names[15], PRMoneDefaults), //  Apply rigid transform
    PRM_Template(),
};


bool
SOP_PCAlign::updateParmsFlags()
{
    bool changed         = SOP_Node::updateParmsFlags();
    const bool method    = evalInt("alignmethod", 0, 0);
    const int usepenalty = evalInt("usepenalty", 0, 0);

    using ParmNameVector = std::vector<const std::string>;

    const ParmNameVector all_parms
    {
        "maxiterations","maxouteriter","maxinneriter", "stopcritera",  
        "usepenalty", "pnorm", "penaltyweight", "penaltyfactor",
        "maxpenalty", "reweightfunc", "doscaling", "doreflections", 
        "docorrespondance", "outliers", "dotransform"
    };

    const ParmNameVector sicp_parms{"maxiterations", "maxouteriter", 
        "maxinneriter", "stopcritera", "usepenalty"};

    const ParmNameVector penalty_parms{"pnorm", "penaltyweight", 
        "penaltyfactor", "maxpenalty"};

    const ParmNameVector reweighted_parms{"maxiterations", "maxouteriter", 
       "pnorm", "maxouteriter", "reweightfunc", "stopcritera"};

    if (method == ALIGN_METHOD::RIGID) {
        for(const auto & parm: all_parms) 
            changed |= enableParm(parm.c_str(), 0);
    } 
    else if (method == ALIGN_METHOD::SPARSE_ICP) 
    {
        for(const auto & parm: all_parms) 
            changed |= enableParm(parm.c_str(), 0);
        for(const auto & parm: sicp_parms) 
           changed |= enableParm(parm.c_str(), 1);

    } 
    // else if (method == ALIGN_METHOD::REWEIGHTED_ICP)

    //     if (usepenalty == 1) 
    //        for(const auto & parm: penalty_parms) 
    //             setVisibleState(parm.c_str(), true);
    //     else 
    //          for(const auto & parm: penalty_parms) 
    //             setVisibleState(parm.c_str(), false);
    // } else if (method == ALIGN_METHOD::REWEIGHTED_ICP)
    // {
    //     for(const auto & parm: reweighted_parms) 
    //             setVisibleState(parm.c_str(), true);
    // }

    #if 0
    #endif
    return changed;
}

void SOP_PCAlign::add_detail_array(const Eigen::MatrixXd & matrix, const char* attr_name)
{
    assert(matrix.rows()*matrix.cols() == 16);
    GA_Attribute * xform_detail_attrib = gdp->addFloatArray(GA_ATTRIB_DETAIL, attr_name , 1);
    const GA_AIFNumericArray * w_aif = xform_detail_attrib->getAIFNumericArray();
    UT_FprealArray array;
    array.setSize(matrix.rows()*matrix.cols());
    for(int i=0; i< matrix.rows(); ++i)
        for(int j=0; j<matrix.cols(); ++j)
            array(i*matrix.rows()+j) = matrix(i,j); 
    
    w_aif->set(xform_detail_attrib, 0, array);
    xform_detail_attrib->bumpDataId();
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

    const int   align_method    = ALIGNMETHOD();
    const int   use_penalty     = USEPENALTY(t);
    const float p_norm          = PNORM(t);   
    const float penalty_weight  = PENALTYWEIGHT(t);         
    const float penalty_factor  = PENALTYFACTOR(t);         
    const float max_penalty     = MAXPENALTY(t);         
    const int   max_iterations  = MAXITER(t);       
    const int   max_outer_iter  = MAXOUTERITER(t);   
    const int   max_inner_iter  = MAXINNERITER(t);   
    const float stop_critera    = STOPCRITERIA(t);  
    const int   weight_func     = WEIGHTFUNC(t);    
    const int   allow_scaling     = DOSCALING(t);    
    const int   allow_reflections = DOREFLECTIONS(t);    
    const int   allow_correspondance = evalInt("docorrespondance", 0, t);    
    const float outliers             = evalFloat("outliers", 0, t);
    const int   apply_transform      = evalInt("dotransform", 0, t);    
   
    Vertices source, target;
    copy_position_to_eigen(gdp, source);
    copy_position_to_eigen(target_gdp, target);

    UT_ASSERT(gdp->getNumPoints() == source.cols());
    UT_ASSERT(target_gdp->getNumPoints() == target.cols());
    UT_ASSERT(source.cols() == target.cols());

    if (align_method == ALIGN_METHOD::RIGID) 
    {
        Eigen::VectorXd weights;
        GA_Attribute * weights_attr = gdp->findFloatTuple(GA_ATTRIB_POINT, "confidence_weight", 1);
        if (!weights_attr) {
            weights = Eigen::VectorXd::Ones(gdp->getNumPoints());
        }
        else {
            if(!copy_float_to_eigen(weights_attr, weights)) {
                return error();
            }
        }

        UT_ASSERT(weights.rows() == gdp->getNumPoints());
        const Eigen::Affine3d xform = RigidMotionEstimator::point_to_point(source, target, weights);

        if(!xform.matrix().allFinite()) {
            addError(SOP_MESSAGE, "Can't compute rigid motion estimator.");
            return error();
        }

        GA_Attribute * xform_attr = gdp->addFloatTuple(GA_ATTRIB_DETAIL, "rigid_xform", 16);
        GA_RWHandleM4 xform_h(xform_attr);

        const UT_Matrix4F m(xform(0,0), xform(1,0), xform(2,0), xform(3,0),
                            xform(0,1), xform(1,1), xform(2,1), xform(3,1),
                            xform(0,2), xform(1,2), xform(2,2), xform(3,2),
                            xform(0,3), xform(1,3), xform(2,3), xform(3,3));
        
        if (xform_h.isValid()) {
            xform_h.set(0, m);
        }

        if (apply_transform) 
        {
            GA_Offset ptoff;
            GA_FOR_ALL_PTOFF(gdp, ptoff) 
            {
                UT_Vector3 pos = gdp->getPos3(ptoff);
                pos *= m;
                gdp->setPos3(ptoff, pos);
            }
        }

        gdp->getP()->bumpDataId();
        return error();

    }
    else if (align_method == ALIGN_METHOD::SPARSE_ICP) 
    {
        SICP::Parameters parms;
        parms.use_penalty = static_cast<bool>(use_penalty);
        parms.p = p_norm;
        parms.mu = penalty_weight;
        parms.alpha = penalty_factor;
        parms.max_mu = max_penalty;
        parms.max_icp = max_iterations;
        parms.max_outer = max_outer_iter;
        parms.max_inner = max_inner_iter;
        parms.stop = stop_critera;

        SICP::point_to_point(source, target, parms);

        GA_Offset ptoff;
        GA_FOR_ALL_PTOFF(gdp, ptoff)
        {   
            const GA_Index i = gdp->pointIndex(ptoff);
            const UT_Vector3 pos(source(0, i), source(1, i), source(2, i));
            gdp->setPos3(ptoff, pos);
        }
    } 
    else if (align_method == ALIGN_METHOD::REWEIGHTED_ICP) 
    {
        ICP::Parameters parms;
        // this crashes... atm
        parms.f = static_cast<ICP::Function>(weight_func);
        parms.p = p_norm;
        parms.max_icp = max_iterations;
        parms.max_outer = max_outer_iter;
        parms.stop = stop_critera;

        ICP::point_to_point(source, target, parms);

        GA_Offset ptoff;
        GA_FOR_ALL_PTOFF(gdp, ptoff)
        {
            const GA_Index i = gdp->pointIndex(ptoff);
            const UT_Vector3 pos(source(0, i), source(1, i), source(2, i));
            gdp->setPos3(ptoff, pos);
        }

    } 
    #if BUILD_WITH_CPD
    else if (align_method == ALIGN_METHOD::RIGID_CPD) 
    {
        Eigen::MatrixXd sourcet = source.transpose();
        Eigen::MatrixXd targett = target.transpose();

        cpd::Rigid rigid;

        rigid.scale(static_cast<bool>(allow_scaling));
        rigid.reflections(static_cast<bool>(allow_reflections));
        rigid.correspondence(static_cast<bool>(allow_correspondance));
        rigid.outliers(outliers);
        rigid.max_iterations(SYSmin(max_iterations,1));
        rigid.tolerance(stop_critera);

        // rigid.gauss_transform(std::move(
        // std::unique_ptr<cpd::GaussTransform>(new cpd::GaussTransformFgt())));

        cpd::RigidResult result = rigid.run(sourcet, targett);

        add_detail_array(result.matrix());

        if (apply_transform) 
        {   
            GA_Offset ptoff;
            GA_FOR_ALL_PTOFF(gdp, ptoff)
            {
                const GA_Index i = gdp->pointIndex(ptoff);
                const UT_Vector3 pos(result.points(i, 0), result.points(i, 1), result.points(i, 2));
                gdp->setPos3(ptoff, pos);
            }
        }
    }
    else if (align_method == ALIGN_METHOD::NONRIGID_CPD) 
    {

        Eigen::MatrixXd sourcet = source.transpose();
        Eigen::MatrixXd targett = target.transpose();

        cpd::Nonrigid nonrigid;

        nonrigid.correspondence(static_cast<bool>(allow_correspondance));
        nonrigid.outliers(outliers);
        nonrigid.max_iterations(SYSmin(max_iterations,1));
        nonrigid.tolerance(stop_critera);

        // nonrigid.gauss_transform(std::move(
        // std::unique_ptr<cpd::GaussTransform>(new cpd::GaussTransformFgt())));

        cpd::NonrigidResult result = nonrigid.run(sourcet, targett);

        GA_Offset ptoff;
        GA_FOR_ALL_PTOFF(gdp, ptoff)
        {
            const GA_Index i = gdp->pointIndex(ptoff);
            const UT_Vector3 pos(result.points(i, 0), result.points(i, 1), result.points(i, 2));
            gdp->setPos3(ptoff, pos);
        }
    }

    #endif

    // Intel FGR relies on PCL, which is to big to add as submodule
    // so we need this switch in case PCL was not found
    #ifdef BUILD_WITH_INTEL_FGR
    else if (align_method == ALIGN_METHOD::INTEL_FGR) 
    {
        FastGlobalRegistration ifgr; 
        std::vector<const GU_Detail*> fgr_sources{2};
        fgr_sources[0] = gdp; fgr_sources[1] = target_gdp;

        for (const auto geo: fgr_sources) 
        {
            // these are shared_ptr:
            sop_pcl::Points::Ptr       points (new sop_pcl::Points);
            sop_pcl::Normals::Ptr      normals(new sop_pcl::Normals);
            sop_pcl::FeatureHists::Ptr fpfhs(new sop_pcl::FeatureHists);

            sop_pcl::gdp_to_pcl(geo, points);
            sop_pcl::estimate_normals(points, normals);
            sop_pcl::compute_fpfh(points, normals, fpfhs);

            std::vector<Eigen::Vector3f> positions(points->size());
            std::vector<Eigen::VectorXf> features(fpfhs->size());

            // both are vectors of eigen Vectors, should be easy to
            // cheat and swap buffers...    
            for (int i=0; i<points->size(); ++i) 
            {
                const pcl::PointXYZ & point        = points->points[i];
                const pcl::FPFHSignature33 & feat  = fpfhs->points[i];
                const Eigen::Vector3f epoint(point.x, point.y, point.z);
                Eigen::VectorXf efeat; efeat.resize(33);
                for (int j=0; j<33; ++j) {
                    efeat(j) = feat.histogram[j];
                }
                positions[i] = epoint;
                features[i]  = efeat;
            }

            ifgr.pointcloud_.emplace_back(positions);
            ifgr.features_.emplace_back(features);
        }
        ifgr.NormalizePoints();
        ifgr.AdvancedMatching();
        ifgr.OptimizePairwise(true, ITERATION_NUMBER); 

        const Eigen::Matrix4f t = ifgr.GetTrans();
        const Eigen::MatrixXd xform = t.cast<double>();
        add_detail_array(xform);

        if (apply_transform)
        {
            const UT_Matrix4F m(t(0,0), t(0,1), t(0,2), t(0,3),
                                t(1,0), t(1,1), t(1,2), t(1,3),
                                t(2,0), t(2,1), t(2,2), t(2,3),
                                t(3,0), t(3,1), t(3,2), t(3,3));
            GA_Offset ptoff;
            GA_FOR_ALL_PTOFF(gdp, ptoff) {
                UT_Vector3 pos = gdp->getPos3(ptoff);
                pos *= m;
                gdp->setPos3(ptoff, pos);
            }
        
        }

        // GA_RWHandleV3 norm_h(gdp->addFloatTuple(GA_ATTRIB_POINT, "N", 3));
        // GA_FOR_ALL_PTOFF(gdp, ptoff) {
        //     const GA_Index i = gdp->pointIndex(ptoff);
        //     const pcl::Normal & n = normals->at(i);
        //     const UT_Vector3 normal(n.normal_x, n.normal_y, n.normal_z);
        //     norm_h.set(ptoff, normal);
        // }
    }
    #endif


    gdp->getP()->bumpDataId();
    return error();
}
