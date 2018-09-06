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

        // Performs tracking of mappoint found in neighbouring KFs
        int SearchByPhoto(Frame &CurrentFrame, const vector<MapPoint*> &vpMapPoints);

    private:

        // Computes distance from point3D to plane defined by normal
        double getDistanceToPlane(const Eigen::Vector3d &point3D, const Eigen::Vector3d &normal);

        // Computes the inverse pose
        cv::Mat getInversePose(cv::Mat Tcw);

        // Creates camera matrix from (fx, y, cx, cy)
        Eigen::Matrix3d getCameraMatrix(float fx, float fy, float cx, float cy);

        // Computes homography
        Eigen::Matrix3d computeHomography(Eigen::Matrix4d Tba, Eigen::Vector3d n, double d, Eigen::Matrix3d Ka, Eigen::Matrix3d Kb);

        // Gets the subpixel value using bilinear interpolation
        double getSubpixImageValue(double u, double v, std::vector< std::vector< float> > &image);

        // Computes the difference between patches with affine parameters (a,b) esitmation with (a*ref + b = cur)
        double computePatchDiffAffine(const std::vector< double> &refPatch,  const std::vector< double> &curPatch);

        // Computes the difference between patches with mean subtraction
        double computePatchDiffAvg(const std::vector< double> &refPatch,  const std::vector< double> &curPatch);

        // Add tracked map point to currently observed
        //      It simulates the correct matching be placing the keypoint in the projected location
        void addTrackedMapPoint(Frame &CurrentFrame, MapPoint *pMP, cv::KeyPoint kp, double currentU, double currentV);

        bool trackMapPoint(MapPoint *pMP, Frame &CurrentFrame,
                Eigen::Vector3d featureInLast, Eigen::Matrix4d Tba, Eigen::Matrix3d Ka,
                g2o::imgStr *lastImage, cv::KeyPoint kp);

        /// Variables

        // Defines the neighbourhood used in comparisons
        std::vector< std::pair<double, double> > neighbours;

        // Photometric error threshold to rescue feature
        double photoThreshold;

    };
} // namespace ORB_SLAM
#endif //ORB_SLAM2_PHOTOTRACKER_H
