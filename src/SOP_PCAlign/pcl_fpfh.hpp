#pragma once
#include <pcl/point_types.h>
#include <pcl/features/normal_3d.h>
#include <pcl/features/fpfh.h>

namespace pcalign 
{
namespace sop_pcl {
    
typedef pcl::PointCloud<pcl::Normal> Normals;
typedef pcl::PointCloud<pcl::PointXYZ> Points;
typedef pcl::PointCloud<pcl::FPFHSignature33> FeatureHists;

inline void gdp_to_pcl(const GU_Detail * gdp, Points::Ptr points) {

    GA_Offset ptoff;
    GA_FOR_ALL_PTOFF(gdp, ptoff) {
        const UT_Vector3 pos = gdp->getPos3(ptoff);
        pcl::PointXYZ point(pos.x(), pos.y(), pos.z());
        points->push_back(point);
    }
}

void estimate_normals(const Points::Ptr points, Normals::Ptr normals, 
                      const float radius = 0.03) {

    // Create the normal estimation class, and pass the input dataset to it
    pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> normalEstimator;
    normalEstimator.setInputCloud(points);
    // Create an empty kdtree representation, and pass it to the normal estimation object.
    // Its content will be filled inside the object, based on the given input dataset (as no other search surface is given).
    pcl::search::KdTree<pcl::PointXYZ>::Ptr tree(new pcl::search::KdTree<pcl::PointXYZ> ());
    normalEstimator.setSearchMethod(tree);
    // Use all neighbors in a sphere of radius 3cm
    normalEstimator.setRadiusSearch(radius);
    // Compute the features
    normalEstimator.compute (*normals);
}

void compute_fpfh(const Points::Ptr points, const Normals::Ptr normals, 
                  FeatureHists::Ptr fpfhs,  const float radius = 0.05) {

    // Create the FPFH estimation class, and pass the input dataset+normals to it
    pcl::FPFHEstimation<pcl::PointXYZ, pcl::Normal, pcl::FPFHSignature33> fpfh;
    fpfh.setInputCloud(points);
    fpfh.setInputNormals(normals);
    // Create an empty kdtree representation, and pass it to the FPFH estimation object.
    // Its content will be filled inside the object, based on the given input dataset (as no other search surface is given).
    pcl::search::KdTree<pcl::PointXYZ>::Ptr tree (new pcl::search::KdTree<pcl::PointXYZ>);
    fpfh.setSearchMethod(tree);
    // Use all neighbors in a sphere of radius 5cm
    // IMPORTANT: the radius used here has to be larger than the radius used to 
    // estimate the surface normals!!!
    fpfh.setRadiusSearch (radius);
    // Compute the features
    fpfh.compute (*fpfhs);
}
}// end of pcl space
}// end of pcalign space