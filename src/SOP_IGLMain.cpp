#include <UT/UT_DSOVersion.h>
#include <OP/OP_OperatorTable.h>
#include "SOP_IGLUVproject.hpp"
#include "SOP_IGLDiscreteGeo.hpp"

using namespace SOP_IGL;

void
newSopOperator(OP_OperatorTable *table)
{
    table->addOperator(new OP_Operator(
        "IGLDiscreteGeo",
        "IGL Discrete Geometry",
        SOP_IGLDiscreteGeometry::myConstructor,
        SOP_IGLDiscreteGeometry::myTemplateList,
        1,
        1,
        0));


    table->addOperator(new OP_Operator(
        "IGLUVProject",
        "IGL UV Project",
        SOP_IGLUVproject::myConstructor,
        SOP_IGLUVproject::myTemplateList,
        1,
        1,
        0));
}