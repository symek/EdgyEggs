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
    CPD,
};

 enum WEIGHT_FUNC {
        PNORM,
        TUKEY,
        FAIR,
        LOGISTIC,
        TRIMMED,
        NONE,
    };


bool copy_position_to_eigen(const GU_Detail * gdp, Vertices & matrix ) 
{
    matrix.resize(Eigen::NoChange, gdp->getNumPoints());

    GA_Offset ptoff;
    GA_FOR_ALL_PTOFF(gdp, ptoff)
    {
        const UT_Vector3 pos = gdp->getPos3(ptoff);
        const GA_Index   idx = gdp->pointIndex(ptoff);
        matrix(0, idx) = pos.x();
        matrix(1, idx) = pos.y();
        matrix(2, idx) = pos.z();
    }
    return true;
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

};

} // End pcallign namespace

