#pragma once 
#include <SOP/SOP_Node.h>

namespace pcalign {


enum ALIGN_METHOD {
    RIGID,
    SPARSE_ICP,
    REWEIGHTED_ICP,
    INTEL_FGR,
};

 enum WEIGHT_FUNC {
        PNORM,
        TUKEY,
        FAIR,
        LOGISTIC,
        TRIMMED,
        NONE,
    };

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

    int     METHOD()              { return evalInt("method", 0, 0); }
    int     USE_PENALTY(fpreal t) { return evalInt("usepenalty", 0, t); }
    fpreal  P_NORM(fpreal t)      { return evalFloat("pnorm", 0, t); }
    fpreal  MU(fpreal t)          { return evalFloat("mu", 0, t); }
    fpreal  ALPHA(fpreal t)       { return evalFloat("alpha", 0, t); }
    fpreal  MAX_MU(fpreal t)      { return evalFloat("maxmu", 0, t); }
    int     MAX_ICP(fpreal t)     { return evalInt("maxicp", 0, t); }
    int     MAX_OUTER(fpreal t)   { return evalInt("maxouter", 0, t); }
    int     MAX_INNER(fpreal t)   { return evalInt("maxinner", 0, t); }
    fpreal  STOP(fpreal t)       { return evalFloat("stop", 0, t); }
    int     WEIGHTFUNC(fpreal t) { return evalInt("weightfunc", 0, t); }

};

} // End pcallign namespace

