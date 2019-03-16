#pragma once

#include <GA/GA_SplittableRange.h>
#include <GA/GA_Range.h>
#include <GA/GA_PageIterator.h>
#include <GA/GA_PageHandle.h>
#include <SOP/SOP_Node.h>
#include <UT/UT_Interrupt.h>

// FIXME: remove from header...
#include <eltopo.h>
#include <surftrack.h>
#include <dynamicsurface.h>
#include "vec.h"
#include <subdivisionscheme.h>

#include <time.h>

namespace SOP_ELTOPO {

enum eltopo_subdiv_schemes {
    MIDPOINT = 0,
    BUTTERFLY,
    MODIFIED_BUTTERFLY,
    QUADRATIC_ERROR_MIN,
};

typedef std::map<int, std::set<int> > UniqueEdges;


class SOP_Eltopo : public SOP_Node
{
public:
         SOP_Eltopo(OP_Network *net, const char *name, OP_Operator *op);
    virtual ~SOP_Eltopo();

    static PRM_Template      myTemplateList[];
    static OP_Node      *myConstructor(OP_Network*, const char *,
                                OP_Operator *);

    /// This method is created so that it can be called by handles.  It only
    /// cooks the input group of this SOP.  The geometry in this group is
    /// the only geometry manipulated by this SOP.
    // virtual OP_ERROR         cookInputGroups(OP_Context &context, 
                        // int alone = 0);

protected:
    /// Method to cook geometry for the SOP
    virtual OP_ERROR         cookMySop(OP_Context &context);

private:
    // void    getGroups(UT_String &str)         { evalString(str, "group", 0, 0); }
    int     COLLISIONSAFETY(fpreal t)             { return evalInt("collisionsafety", 0, t); }
    int     ALLOWTOPOLOGYCHANGE(fpreal t)         { return evalInt("allowtopologychange", 0, t); }
    int     PERFORMIMPROVMENT(fpreal t)           { return evalInt("performimprovment", 0, t); }
    int     ALLOWVERTEXMOVEMENT(fpreal t)         { return evalInt("allowvertexmovement", 0, t); }
    int     USEFRACTION(fpreal t)                 { return evalInt("usefraction", 0, t); }
    int     USECURVATUREWHENSPLITING(fpreal t)    { return evalInt("usecurvaturewhenspliting", 0, t); }
    int     USECURVATUREWHENCOLLAPSING(fpreal t)  { return evalInt("usecurvaturewhencollapsing", 0, t); }
    int     ALLOWNONNANIFOLD(fpreal t)            { return evalInt("allownonnanifold", 0, t); }

    fpreal  MINEDGELENGTH(fpreal t)              { return evalFloat("minedgelength", 0, t); }
    fpreal  MAXEDGELENGTH(fpreal t)              { return evalFloat("maxedgelength", 0, t); }
    fpreal  MINTRIANGLEAREA(fpreal t)            { return evalFloat("mintrianglearea", 0, t); }
    fpreal  MAXVOLUMECHANGE(fpreal t)            { return evalFloat("maxvolumechange", 0, t); }
    fpreal  EDGEFLIPMINLENGTH(fpreal t)          { return evalFloat("edgeflipminlength", 0, t); }
    fpreal  MERGEPROXIMITYEPS(fpreal t)          { return evalFloat("mergeproximityeps", 0, t); }
    int     SUBDIVISIONSCHEME(fpreal t)          { return evalInt("subdivisionscheme", 0, t); }
    // fpreal  LAPLACIAN(fpreal t)              { return evalFloat("laplacian", 0, t); }


    /// This is the group of geometry to be manipulated by this SOP and cooked
    /// by the method "cookInputGroups".
    const GA_PointGroup *myGroup;
    // std::shared_ptr<SurfTrackInitializationParameters> parms;
    // std::shared_ptr<SurfTrack>                         surface_tracker;

    // std::shared_ptr<ElTopoGeneralOptions>          general_options;
    // std::shared_ptr<ElTopoStaticOperationsOptions> static_options;
    // std::shared_ptr<ElTopoIntegrationOptions>      integration_options;
    // std::shared_ptr<ElTopoDefragInformation>       defrag_options;
    // std::shared_ptr<ElTopoMesh>                    my_mesh;
};


} // End SOP_Eltopo namespace
