#ifndef _PHOTOMETRIC_ERROR_FUNCTIONS_H
#define _PHOTOMETRIC_ERROR_FUNCTIONS_H


#include <vector>
#include <Eigen/Eigen>
#include<opencv2/core/core.hpp>

namespace photo {

    // Image pyramid
    struct imgStr {
        float imageScale;
        std::vector< std::vector< float> > image;
        std::vector< std::vector< Eigen::Vector2f > > gradient;
    };

    // Gets the subpixel value using bilinear interpolation
    double getSubpixImageValue(double u, double v, std::vector <std::vector<float>> &image);
    double getSubpixImageValue(double u, double v, cv::Mat image);

    // Computes distance from point3D to plane defined by normal
    double getDistanceToPlane(const Eigen::Vector3d &point3D, const Eigen::Vector3d &normal);

    // Creates camera matrix from (fx, y, cx, cy)
    Eigen::Matrix3d getCameraMatrix(float fx, float fy, float cx, float cy);

    // Computes homography
    Eigen::Matrix3d computeHomography(Eigen::Matrix4d Tba, Eigen::Vector3d n, double d, Eigen::Matrix3d Ka, Eigen::Matrix3d Kb);

    // Computes the difference between patches with affine parameters (a,b) esitmation with (a*ref + b = cur)
    double computePatchDiffAffine(const std::vector< double> &refPatch,  const std::vector< double> &curPatch);

    // Computes the difference between patches with mean subtraction
    double computePatchDiffAvg(const std::vector< double> &refPatch,  const std::vector< double> &curPatch);
};

#endif // _PHOTOMETRIC_ERROR_FUNCTIONS_H