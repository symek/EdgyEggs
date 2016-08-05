#pragma once

#include <GA/GA_SplittableRange.h>
#include <GA/GA_Range.h>
#include <GA/GA_PageIterator.h>
#include <GA/GA_PageHandle.h>
#include <SOP/SOP_Node.h>
#include <UT/UT_Interrupt.h>

#include <time.h>

class ShapeOpSolver;

namespace SOP_SHAPEOP {

typedef std::map<int, std::set<int> > UniqueEdges;

enum ShapeOp_ConstraintType 
{
    EDGE_STRAIN,
    TRIANGLE_STRAIN,
    LINE,
};

class SOP_ShapeOp : public SOP_Node
{
public:
         SOP_ShapeOp(OP_Network *net, const char *name, OP_Operator *op);
    virtual ~SOP_ShapeOp();

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
    int     MAXITER(fpreal t)                { return evalInt("maxiter", 0, t); }
    fpreal  CLOSENESS(fpreal t)              { return evalFloat("closeness", 0, t); }
    fpreal  EDGESTRAIN(fpreal t)             { return evalFloat("edgestrain", 0, t); }


    /// This is the group of geometry to be manipulated by this SOP and cooked
    /// by the method "cookInputGroups".
    const GA_PointGroup *myGroup;
    ShapeOpSolver       *mySolver;
     
};


} // End SOP_SHAPEOP namespace

