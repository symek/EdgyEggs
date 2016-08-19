#include <GU/GU_Detail.h>
#include <GU/GU_PrimPoly.h>
#include <GU/GU_PrimTetrahedron.h>
#include <GEO/GEO_PolyCounts.h>
#include <igl/readMESH.h>

#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <Eigen/Sparse>

// #define SHAPEOP_HEADER_ONLY
// #include <libShapeOp/api/API.cpp>

typedef std::map<int,std::set<int> > UniqueEdges;

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


struct Mesh
{
  Eigen::MatrixXd V,U;
  Eigen::MatrixXi T,F;
} low;

int main(int argc, const char* argv[])
{

    if (argc < 3)
        return 1;
   
    GU_Detail gdp;
    const char *mesh     = argv[1];
    const char *meshbgeo = argv[2];
  


  if(!igl::readMESH(mesh, low.V, low.T, low.F))
  {
    std::cout<<"failed to load mesh\n";
    return 1;
}
 

  std::cout << "Vertices:  " << low.V.rows() << "x" << low.V.cols() << "\n";
  std::cout << "Tetrahed:  " << low.T.rows() << "x" << low.T.cols() << "\n";
  std::cout << "Triangels: " << low.F.rows() << "x" << low.F.cols() << "\n";


    const uint npoints = low.V.rows();
    const uint ntets   = low.T.rows();
    const uint ntris   = low.F.rows();

    std::vector<int> tetpointnumbers;
    std::vector<int> tripointnumbers;


   
    for(int i=0; i< ntets; ++i) {    
        tetpointnumbers.push_back(low.T(i, 0));
        tetpointnumbers.push_back(low.T(i, 1));
        tetpointnumbers.push_back(low.T(i, 2));
        tetpointnumbers.push_back(low.T(i, 3));
    }

    for(int i=0; i< ntris; ++i) {    
            tripointnumbers.push_back(low.F(i, 0));
            tripointnumbers.push_back(low.F(i, 1));
            tripointnumbers.push_back(low.F(i, 2));
    }




    const GA_Offset startpt  = gdp.appendPointBlock((GA_Size)npoints);

    {
        GA_Offset ptoff;
        GA_FOR_ALL_PTOFF(&gdp, ptoff) 
        {
            const int ptidx = gdp.pointIndex(ptoff);
            UT_ASSERT(ptidx < npoints);
            const UT_Vector3 pos(low.V(ptidx, 0), low.V(ptidx, 2), low.V(ptidx, 1)); // note: switch axes.
            gdp.setPos3(ptoff, pos);
        }
    }

    GEO_PolyCounts polygonsizes; 
    polygonsizes.append(3, ntris);
    GU_PrimTetrahedron::buildBlock(&gdp, startpt, (GA_Size)npoints, ntets, &tetpointnumbers[0]);
    GU_PrimPoly::buildBlock(&gdp, startpt, (GA_Size)npoints, polygonsizes, &tripointnumbers[0]);

    gdp.save(meshbgeo, 0);

    return 0;

    }

