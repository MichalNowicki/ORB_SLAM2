//
// Created by mnowicki on 04.09.18.
//

#ifndef ORB_SLAM2_PHOTOTRACKER_H
#define ORB_SLAM2_PHOTOTRACKER_H


#include"MapPoint.h"
#include"KeyFrame.h"
#include"Frame.h"

#include "Thirdparty/g2o/g2o/types/types_six_dof_photo.h"
#include <opencv2/core/eigen.hpp>

namespace ORB_SLAM2 {

    class PhotoTracker {
    public:

        // Initializes photometric tracking of features
        PhotoTracker(double photoThreshold = 20);

        // Performs tracking of feature motion
        int SearchByPhoto(Frame &CurrentFrame, const Frame &LastFrame);

    private:

        // Computes distance from point3D to plane defined by normal
        double getDistanceToPlane(const cv::Mat &point3D, const cv::Mat &normal);

        // Normalizes to (u,v,1)
        cv::Mat normalize2D(cv::Mat p);

        // Computes the inverse pose
        cv::Mat getInversePose(cv::Mat Tcw);

        // Creates camera matrix from (fx, y, cx, cy)
        cv::Mat getCameraMatrix(float fx, float fy, float cx, float cy);

        // Computes homography
        cv::Mat computeHomography(cv::Mat Tba, cv::Mat n, double d, cv::Mat Ka, cv::Mat Kb);

        // Gets the subpixel value using bilinear interpolation
        double getSubpixImageValue(double u, double v, std::vector< std::vector< float> > &image);

        // Computes the difference between patches with affine parameters (a,b) esitmation with (a*ref + b = cur)
        double computePatchDiffAffine(const std::vector< double> &refPatch,  const std::vector< double> &curPatch);

        // Computes the difference between patches with mean subtraction
        double computePatchDiffAvg(const std::vector< double> &refPatch,  const std::vector< double> &curPatch);

        // Defines the neighbourhood used in comparisons
        std::vector< std::pair<double, double> > neighbours;


        ///

        // Photometric error threshold to rescue feature
        double photoThreshold;

    };
} // namespace ORB_SLAM
#endif //ORB_SLAM2_PHOTOTRACKER_H
