//
// Created by mnowicki on 04.09.18.
//

#ifndef ORB_SLAM2_PHOTOTRACKER_H
#define ORB_SLAM2_PHOTOTRACKER_H


#include"MapPoint.h"
#include"KeyFrame.h"
#include"Frame.h"

#include "photometricErrorFunctions.h"
#include "g2oPhotoError.h"
#include <opencv2/core/eigen.hpp>
#include <utility>

#include <fstream>

namespace ORB_SLAM2 {

    class PhotoTracker {
    public:

        // Initializes photometric tracking of features
        PhotoTracker(double photoThreshold = 20, int kltMaxIterations=30, double kltEPS=0.01, double kltZNCCThr=0.9,
                int kltPatchSize=9, bool verbose = false);

        // Performs tracking of feature motion
        int SearchByPhoto(Frame &CurrentFrame, Frame &LastFrame);

        // Performs tracking of mappoint found in neighbouring KFs
        int SearchByPhoto(Frame &CurrentFrame, const vector<MapPoint*> &vpMapPoints);

        // KLT
        std::pair<int,int> SearchByKLT(Frame &CurrentFrame, Frame &LastFrame);
        int SearchByKLT(Frame &CurrentFrame, const vector<MapPoint*> &vpMapPoints);

    private:

        // Add tracked map point to currently observed
        //      It simulates the correct matching be placing the keypoint in the projected location
        void addTrackedMapPoint(Frame &CurrentFrame, MapPoint *pMP, cv::KeyPoint kp, double currentU, double currentV);

        bool trackMapPoint(MapPoint *pMP, Frame &CurrentFrame,
                Eigen::Vector3d featureInLast, Eigen::Matrix4d Tba, Eigen::Matrix3d Ka,
                photo::imgStr *lastImage, cv::KeyPoint kp);

        double computeZNCC(Frame &LastFrame, Frame &CurrentFrame, cv::Point2f &prevPts, cv::Point2f &nextPts, const int &patchSize);

        /// Variables

        // Defines the neighbourhood used in comparisons
        std::vector< std::pair<double, double> > neighbours;

        // Photometric error threshold to rescue feature
        double photoThreshold;

        // KLT threshold
        int kltMaxIterations;
        double kltEPS;
        double kltZNCCThr;
        int patchSize;

        bool verbose;


    };
} // namespace ORB_SLAM
#endif //ORB_SLAM2_PHOTOTRACKER_H
