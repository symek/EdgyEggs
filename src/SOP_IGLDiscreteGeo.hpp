#pragma once

#include <GA/GA_SplittableRange.h>
#include <GA/GA_Range.h>
#include <GA/GA_PageIterator.h>
#include <GA/GA_PageHandle.h>

#include <SOP/SOP_Node.h>


#include <time.h>


namespace SOP_IGL {

class Timer 
{
    private:
        double begTime;
    public:
        void start()
        {
            begTime = clock();
        }

        double current() 
        {
        //int threads = UT_Thread::getNumProcessors();
            return (difftime(clock(), begTime) / CLOCKS_PER_SEC);// / (double) threads;
        }

        bool isTimeout(double seconds)
        {
            return seconds >= current();
        }
};


void detail_to_eigen(const GU_Detail &gdp, Eigen::MatrixXd &points, Eigen::MatrixXi &faces)
{
  GA_Offset ptoff;
  GA_FOR_ALL_PTOFF(&gdp, ptoff)
  {
      const UT_Vector3 pos = gdp.getPos3(ptoff);
      points(static_cast<uint>(ptoff), 0) = pos.x();
      points(static_cast<uint>(ptoff), 1) = pos.y(); 
      points(static_cast<uint>(ptoff), 2) = pos.z();
  }

  GA_Iterator it(gdp.getPrimitiveRange());
  for (; !it.atEnd(); ++it) {
        const GEO_Primitive *prim = gdp.getGEOPrimitive(*it);
        GA_Primitive::const_iterator vt;
        int vertex_index = 0;
        for (prim->beginVertex(vt); !vt.atEnd(); ++vt) {
            const GA_Offset voff = vt.getPointOffset();
            faces(prim->getMapIndex(), SYSmin(vertex_index, 2)) = static_cast<int>(voff);
            vertex_index++;
        }
    }

}

void eigen_to_detail(const Eigen::MatrixXd &points, const Eigen::MatrixXi &faces, GU_Detail &gdp)
{

  GA_Offset ptoff;
  UT_Vector3 pos; 
  gdp.appendPointBlock((GA_Size)points.rows());
  GA_FOR_ALL_PTOFF(&gdp, ptoff)
  {
      pos.x() = points(static_cast<uint>(ptoff), 0);
      pos.y() = points(static_cast<uint>(ptoff), 1); 
      pos.z() = points(static_cast<uint>(ptoff), 2);
      gdp.setPos3(ptoff, pos);
  }

  // TODO: optimize!
  for (uint i = 0; i < faces.rows(); ++i) {
    GU_PrimPoly *prim = GU_PrimPoly::build(&gdp, 0, false, false);
    prim->appendVertex((GA_Index)faces(i, 0));
    prim->appendVertex((GA_Index)faces(i, 1));
    prim->appendVertex((GA_Index)faces(i, 2)); 
  }

}


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
    void    getGroups(UT_String &str){ evalString(str, "group", 0, 0); }
    void    MODEL(UT_String &str)    { evalString(str, "model", 0, 0); }
    void    TERM(UT_String &str)     { evalString(str, "term", 0, 0); }
    fpreal  QCOEF(fpreal t)          { return evalFloat("qcoef", 0, t); }
    fpreal  ZCOEF(fpreal t)           { return evalFloat("zcoef", 0, t); }
    fpreal  RADIUS(fpreal t)    { return evalFloat("radius", 0, t); }
    int     LAYERS(fpreal t)    { return evalInt("layers", 0, t); }
    fpreal  LAMBDA(fpreal t)    { return evalFloat("lambda", 0, t); }
    int     TANGENT(fpreal t)   { return evalInt("tangent", 0, t); }

    /// This is the group of geometry to be manipulated by this SOP and cooked
    /// by the method "cookInputGroups".
    const GA_PointGroup *myGroup;
};


} // End SOP_IGL namespace


// class op_IGLDiscreteGeometry {
// public:
//     op_IGLDiscreteGeometry(const std::string &str_model, const bool tnSpace, GU_Detail *gdp)
//         : mystr_model(str_model),  myGdp(gdp), myTnSpace(tnSpace) {};
//             // Take a SplittableRange (not a GA_Range)
//     void    operator()(const GA_SplittableRange &r) const
//             {
//                 alglib::rbfmodel model;
//                 GA_RWPageHandleV3 handle_P(myGdp->getP());
//                 GA_ROPageHandleV3 handle_U(myGdp, GA_ATTRIB_POINT, "tangentu");
//                 GA_ROPageHandleV3 handle_V(myGdp, GA_ATTRIB_POINT, "tangentv");
//                 GA_ROPageHandleV3 handle_N(myGdp, GA_ATTRIB_POINT, "N");

//                 // Execute storage:
//                 alglib::real_1d_array coord("[0,0,0]");
//                 alglib::real_1d_array result("[0,0,0]");
//                 std::string copy_serialized(mystr_model);
//                 alglib::rbfunserialize(copy_serialized, model);
                
//                 // Iterate over pages in the range
//                 for (GA_PageIterator pit = r.beginPages(); !pit.atEnd(); ++pit)
//                 {
//                     GA_Offset start, end;
//                     // iterate over the elements in the page.
//                     for (GA_Iterator it(pit.begin()); it.blockAdvance(start, end); )
//                     {
//                         // Perform any per-page setup required, then
//                         handle_P.setPage(start); handle_U.setPage(start);
//                         handle_V.setPage(start); handle_N.setPage(start);
//                         for (GA_Offset i = start; i < end; ++i)
//                         {
//                             const UT_Vector3 pos = handle_P.get(i);
//                             const double dp[] = {pos.x(), pos.y(), pos.z()};
//                             coord.setcontent(3, dp);
//                             alglib::rbfcalc(model, coord, result);
//                             UT_Vector3 displace = UT_Vector3(result[0], result[1], result[2]);
//                             if (myTnSpace)
//                             {
//                                 UT_Vector3 u = handle_U.get(i); 
//                                 UT_Vector3 v = handle_V.get(i);
//                                 UT_Vector3 n = handle_N.get(i);
//                                 u.normalize(); v.normalize(); n.normalize();
//                                 UT_Matrix3  b(u.x(), u.y(), u.z(),
//                                               v.x(), v.y(), v.z(),
//                                               n.x(), n.y(), n.z());

//                                 b = b.transposedCopy() * b;
//                                 UT_Vector3 a1(u * b); a1.normalize();
//                                 UT_Vector3 a2(v * b); a2.normalize();
//                                 const float da1 = displace.dot(a1);
//                                 const float da2 = displace.dot(a2);
//                                 displace        = UT_Vector3(a1 * da1 + a2 * da2);
        
//                             }

//                             handle_P.set(i, pos+displace);
//                         }
//                     }
//                 }
//             }
//     private:
//             const std::string &mystr_model;
//             GU_Detail         *myGdp;
//             const bool        myTnSpace; 
// };
// void
// IGLDiscreteGeometryThreaded(const GA_Range &range, const std::string &str_model, \
//     const bool tnSpace, GU_Detail *gdp)
// {
//     // Create a GA_SplittableRange from the original range
//     GA_SplittableRange split_range = GA_SplittableRange(range);
//     UTparallelFor(split_range, op_IGLDiscreteGeometry(str_model, tnSpace, gdp));
// }
