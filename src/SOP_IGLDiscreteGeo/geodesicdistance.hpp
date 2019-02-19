namespace SOP_IGL {

bool compute_exact_geodesic_distance(
	const Eigen::MatrixXd & points, 
	const Eigen::VectorXi & faces,
	const int anchor,
	      Eigen::VectorXd & distances)
{
	Eigen::VectorXi VS,FS,VT,FT;
    VS.resize(1);
    VS << anchor;
    // All vertices are the targets
    VT.setLinSpaced(points.rows(), 0, points.rows()-1);
    igl::exact_geodesic(points, faces, VS, FS, VT, FT, distances);
    return true;
}


} // end of SOP_IGL namespace