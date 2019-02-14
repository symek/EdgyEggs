// ----------------------------------------------------------------------------
// -                       Fast Global Registration                           -
// ----------------------------------------------------------------------------
// The MIT License (MIT)
//
// Copyright (c) Intel Corporation 2016
// Qianyi Zhou <Qianyi.Zhou@gmail.com>
// Jaesik Park <syncle@gmail.com>
// Vladlen Koltun <vkoltun@gmail.com>
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
// IN THE SOFTWARE.
// ----------------------------------------------------------------------------
#pragma once
#include <vector>
#include <flann/flann.hpp>

#include <Eigen/Core>
#include <Eigen/Geometry>

#define DIV_FACTOR          1.4     // Division factor used for graduated non-convexity
#define USE_ABSOLUTE_SCALE  0       // Measure distance in absolute scale (1) or in scale relative to the diameter of the model (0)
#define MAX_CORR_DIST       0.025   // Maximum correspondence distance (also see comment of USE_ABSOLUTE_SCALE)
#define ITERATION_NUMBER    64      // Maximum number of iteration
#define TUPLE_SCALE         0.95    // Similarity measure used for tuples of feature points.
#define TUPLE_MAX_CNT       1000    // Maximum tuple numbers.


using namespace Eigen;
using namespace std;

namespace pcalign {

typedef vector<Vector3f> Points;
typedef vector<VectorXf> Feature;

class FastGlobalRegistration{
public:
    void LoadFeature(const Points& pts, const Feature& feat);
    void ReadFeature(const char* filepath);
    void NormalizePoints();
    void AdvancedMatching();
    void WriteTrans(const char* filepath);
    Matrix4f GetTrans();
    double OptimizePairwise(bool decrease_mu_, int numIter_);

    // containers
    vector<Points> pointcloud_;
    vector<Feature> features_;
private:
    Matrix4f TransOutput_;
    vector<pair<int, int>> corres_;

    // for normalization
    Points Means;
    float GlobalScale;
    float StartScale;

    // some internal functions
    void ReadFeature(const char* filepath, Points& pts, Feature& feat);

    void SearchFLANNTree(flann::Index<flann::L2<float>>* index,
        VectorXf& input,
        std::vector<int>& indices,
        std::vector<float>& dists,
        int nn);
};

} // end of pcalign space
