#pragma once

namespace SOP_IGL  {

void compute_curvature(GU_Detail *gdp, Eigen::MatrixXd &V, Eigen::MatrixXi &F, const int false_colors=0);

void compute_gradient(GU_Detail *gdp, const GA_ROHandleF &sourceAttrib, Eigen::MatrixXd &V, Eigen::MatrixXi &F);

int compute_laplacian(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, 
                      const Eigen::SparseMatrix<double> &L, Eigen::MatrixXd &U);

bool compute_exact_geodesic_distance(const Eigen::MatrixXd & points, const Eigen::VectorXi & faces,
									 const int anchor, Eigen::VectorXd & distances);


} // end of SOP_IGL namespace