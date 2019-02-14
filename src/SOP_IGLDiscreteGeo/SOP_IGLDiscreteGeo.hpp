#pragma once

#include <GA/GA_SplittableRange.h>
#include <GA/GA_Range.h>
#include <GA/GA_PageIterator.h>
#include <GA/GA_PageHandle.h>
#include <SOP/SOP_Node.h>
#include <UT/UT_Interrupt.h>

#include <time.h>

namespace SOP_IGL {


class SOP_IGLDiscreteGeometry : public SOP_Node
{
public:
         SOP_IGLDiscreteGeometry(OP_Network *net, const char *name, OP_Operator *op);
    virtual ~SOP_IGLDiscreteGeometry();

    static PRM_Template      myTemplateList[];
    static OP_Node      *myConstructor(OP_Network*, const char *,
                                OP_Operator *);

    /// This method is created so that it can be called by handles.  It only
    /// cooks the input group of this SOP.  The geometry in this group is
    /// the only geometry manipulated by this SOP.
    virtual OP_ERROR         cookInputGroups(OP_Context &context, 
                        int alone = 0);

protected:
    /// Method to cook geometry for the SOP
    virtual OP_ERROR         cookMySop(OP_Context &context);

private:
    void    getGroups(UT_String &str)        { evalString(str, "group", 0, 0); }
    int     CURVATURE(fpreal t)              { return evalInt("curvature", 0, t); }
    int     FALSE_CURVE_COLORS(fpreal t)     { return evalInt("false_curve_colors", 0, t); }
    int     GRAD_ATTRIB(fpreal t)            { return evalInt("grad_attrib", 0, t); }
    void    GRAD_ATTRIB_NAME(UT_String &str) { evalString(str,"grad_attrib_name", 0, 0); }
    fpreal  LAPLACIAN(fpreal t)              { return evalFloat("laplacian", 0, t);  }
    int     EIGENVECTORS(fpreal t)           { return evalInt("eigenvectors", 0, t); }

    /// This is the group of geometry to be manipulated by this SOP and cooked
    /// by the method "cookInputGroups".
    const GA_PointGroup *myGroup;
};

} // End SOP_IGL namespace

