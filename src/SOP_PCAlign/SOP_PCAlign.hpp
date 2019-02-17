#pragma once 
#include <SOP/SOP_Node.h>

namespace pcalign {

typedef double Scalar;
typedef Eigen::Matrix<Scalar, 3, Eigen::Dynamic> Vertices;

enum ALIGN_METHOD {
    RIGID,
    SPARSE_ICP,
    REWEIGHTED_ICP,
    INTEL_FGR,
    RIGID_CPD,
    NONRIGID_CPD,
};

 enum WEIGHT_FUNC {
        PNORM,
        TUKEY,
        FAIR,
        LOGISTIC,
        TRIMMED,
        NONE,
    };

void copy_float_to_eigen(const GA_Attribute * attr, Eigen::VectorXd & weightsV) 
{
    const GA_Detail & gdp = attr->getDetail();
    GA_ROHandleF  weights_h(attr);

    if (weights_h.isValid())
    {   
        weightsV.resize(gdp.getNumPoints());
        GA_Offset ptoff;
        GA_FOR_ALL_PTOFF(&gdp, ptoff)
        {
            const float w = weights_h.get(ptoff);
            const GA_Index idx = gdp.pointIndex(ptoff);
            weightsV(idx) = static_cast<double>(w);
        }
    }
}

void copy_position_to_eigen(const GU_Detail * gdp, Vertices & matrix) 
{
    // Vertices matrix;
    matrix.resize(Eigen::NoChange, gdp->getNumPoints());

    GA_Offset ptoff;
    GA_FOR_ALL_PTOFF(gdp, ptoff)
    {
        const UT_Vector3 pos = gdp->getPos3(ptoff);
        const GA_Index   idx = gdp->pointIndex(ptoff);
        matrix(0, idx) = static_cast<double>(pos.x());
        matrix(1, idx) = static_cast<double>(pos.y());
        matrix(2, idx) = static_cast<double>(pos.z());
    }
    // return matrix;
}

Eigen::MatrixXd copy_position_to_eigen_rows(const GU_Detail * gdp) 
{
    Eigen::MatrixXd matrix;
    matrix.resize(gdp->getNumPoints(), 3);

    GA_Offset ptoff;
    GA_FOR_ALL_PTOFF(gdp, ptoff)
    {
        const UT_Vector3 pos = gdp->getPos3(ptoff);
        const GA_Index   idx = gdp->pointIndex(ptoff);
        matrix(idx, 0) = static_cast<double>(pos.x());
        matrix(idx, 1) = static_cast<double>(pos.y());
        matrix(idx, 2) = static_cast<double>(pos.z());
    }
    return matrix;
}

class SOP_PCAlign : public SOP_Node
{
public:
    SOP_PCAlign(OP_Network *net, const char *name, OP_Operator *op);
    virtual ~SOP_PCAlign();
    virtual bool updateParmsFlags() override;
    static PRM_Template      myTemplateList[];
    static OP_Node      *myConstructor(OP_Network*, const char *,
                                OP_Operator *);
protected:
    /// Method to cook geometry for the SOP
    virtual OP_ERROR         cookMySop(OP_Context &context) override;
private:

    int     ALIGNMETHOD()         { return evalInt("alignmethod", 0, 0); }
    int     MAXITER(fpreal t)     { return evalInt("maxiterations", 0, t); }
    int     MAXOUTERITER(fpreal t)   { return evalInt("maxouteriter", 0, t); }
    int     MAXINNERITER(fpreal t)   { return evalInt("maxinneriter", 0, t); }
    fpreal  STOPCRITERIA(fpreal t)   { return evalFloat("stopcritera", 0, t); }
    int     USEPENALTY(fpreal t)     { return evalInt("usepenalty", 0, t); }
    fpreal  PNORM(fpreal t)         { return evalFloat("pnorm", 0, t); }
    fpreal  PENALTYWEIGHT(fpreal t) { return evalFloat("penaltyweight", 0, t); }
    fpreal  PENALTYFACTOR(fpreal t) { return evalFloat("penaltyfactor", 0, t); }
    fpreal  MAXPENALTY(fpreal t)    { return evalFloat("maxpenalty", 0, t); }
    int     WEIGHTFUNC(fpreal t)    { return evalInt("reweightfunc", 0, t); }
    int     DOSCALING(fpreal t)     { return evalInt("doscaling", 0, t); }
    int     DOREFLECTIONS(fpreal t) { return evalInt("doreflections", 0, t); }

    void    add_detail_array(const Eigen::MatrixXd & , const char* attr_name="rigid_xform");

};

} // End pcallign namespace

