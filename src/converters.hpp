class GU_Detail;
class GU_PrimPoly;

namespace SOP_IGL { 

void detail_to_eigen(const GU_Detail &gdp, Eigen::MatrixXd &points, Eigen::MatrixXi &faces);
void eigen_to_detail(const Eigen::MatrixXd &points, const Eigen::MatrixXi &faces, GU_Detail &gdp);
void eigen_to_detail_points(const Eigen::MatrixXd &points, GU_Detail &gdp);
void eigen_to_point_attribF(const Eigen::MatrixXd &M, GU_Detail *gdp, const char* attribName);
void eigen_to_point_attribV(const Eigen::MatrixXd &M, GU_Detail *gdp, const char* attribName);

} // end of namespace SOP_IGL