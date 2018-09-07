#include"PhotoTracker.h"

namespace ORB_SLAM2 {


    PhotoTracker::PhotoTracker(double photoThreshold) : photoThreshold(photoThreshold) {

        //   x
        //  xxx
        // xxxxx - neighbours used in photo check
        //  xxx
        //   x
        neighbours.push_back(make_pair(0, -2));

        neighbours.push_back(make_pair(-1, -1));
        neighbours.push_back(make_pair(0, -1));
        neighbours.push_back(make_pair(1, -1));

        neighbours.push_back(make_pair(-2, 0));
        neighbours.push_back(make_pair(-1, 0));
        neighbours.push_back(make_pair(0, 0));
        neighbours.push_back(make_pair(1, 0));
        neighbours.push_back(make_pair(2, 0));

        neighbours.push_back(make_pair(-1, 1));
        neighbours.push_back(make_pair(0, 1));
        neighbours.push_back(make_pair(1, 1));

        neighbours.push_back(make_pair(0, 2));
    }

    double PhotoTracker::getDistanceToPlane(const Eigen::Vector3d &point3D, const Eigen::Vector3d &normal) {
        return -normal.transpose() * point3D;
    }

    Eigen::Matrix3d PhotoTracker::computeHomography(Eigen::Matrix4d Tba, Eigen::Vector3d n, double d, Eigen::Matrix3d Ka, Eigen::Matrix3d Kb) {
        // Getting R,t
        Eigen::Matrix3d R21 = Tba.block<3,3>(0,0);
        Eigen::Vector3d t21 = Tba.block<3,1>(3,0);

        Eigen::Matrix3d H = Ka * (R21 - t21 * n.transpose() / d) * Kb.inverse();

        // Homography
        return H;
    }

    Eigen::Matrix3d PhotoTracker::getCameraMatrix(float fx, float fy, float cx, float cy) {
        Eigen::Matrix3d cameraMatrix = Eigen::Matrix3d::Identity();
        cameraMatrix(0,0) = fx;
        cameraMatrix(0,2) = cx;
        cameraMatrix(1,1) = fy;
        cameraMatrix(1,2) = cy;
        cameraMatrix(1,2) = cy;

        return cameraMatrix;
    }

    double PhotoTracker::getSubpixImageValue(double u, double v, std::vector<std::vector<float> > &image) {

        const double xInt = int(u), yInt = int(v);
        const double xSub = u - xInt, ySub = v - yInt;

        const double topLeft = (1.0 - xSub) * (1.0 - ySub);
        const double topRight = xSub * (1.0 - ySub);
        const double bottomLeft = (1.0 - xSub) * ySub;
        const double bottomRight = xSub * ySub;


        if (yInt < 0 || xInt < 0 || yInt + 1 >= image.size() || xInt + 1 >= image[0].size()) {
            return -1;
        }

        return topLeft * image[yInt][xInt] +
               topRight * image[yInt][xInt + 1] +
               bottomLeft * image[yInt + 1][xInt] +
               bottomRight * image[yInt + 1][xInt + 1];
    }

    double PhotoTracker::computePatchDiffAffine(const std::vector< double> &refPatch,  const std::vector< double> &curPatch) {
        // Least squares solution to linear brightness change: (alfa * curPatch + beta = refPatch)
        Eigen::MatrixXf A = Eigen::MatrixXf::Ones(refPatch.size(), 2);
        Eigen::VectorXf b = Eigen::VectorXf::Zero(curPatch.size());
        for (int i=0;i<refPatch.size();i++) {
            A(i,0) = curPatch[i];
            b(i) = refPatch[i];
        }
        Eigen::VectorXf brightness = A.colPivHouseholderQr().solve(b);

        std::cout << "Affine: " << brightness(0) << " " << brightness(1) << std::endl;

        // Computing the difference
        double rmse = 0;
        for (int i=0;i<refPatch.size();i++) {
            float diff = brightness(0) * curPatch[i] + brightness(1) - refPatch[i];
            rmse += diff*diff;
        }
        return sqrt(rmse / refPatch.size());
    }

    double PhotoTracker::computePatchDiffAvg(const std::vector< double> &refPatch,  const std::vector< double> &curPatch) {
        double averageRefPatch = std::accumulate( refPatch.begin(), refPatch.end(), 0.0)/refPatch.size();
        double averageCurPatch = std::accumulate( curPatch.begin(), curPatch.end(), 0.0)/curPatch.size();

        // Computing the difference
        double rmse = 0;
        for (int i=0;i<refPatch.size();i++) {
            float diff = (curPatch[i]-averageCurPatch) - (refPatch[i] - averageRefPatch);
            rmse += diff*diff;
        }
        return sqrt(rmse / refPatch.size());
    }


    int PhotoTracker::SearchByPhoto(Frame &CurrentFrame, Frame &LastFrame) {
        // Number of tracked features
        int nmatches = 0;

        // Relative transformation between frames in Eigen
        Eigen::Matrix4d Taw = LastFrame.getPose();
        Eigen::Matrix4d Tbw = CurrentFrame.getPose();
        Eigen::Matrix4d Tba = Tbw * Taw.inverse();

        // For all map points observed in last frame
        for (int i = 0; i < LastFrame.N; i++) {
            MapPoint *pMP = LastFrame.mvpMapPoints[i];

            if (pMP) {
                if (!LastFrame.mvbOutlier[i]) {

                    // Lets check if this feature was matched using desciptors
                    bool featureMatched = false;
                    for (int j=0;j<CurrentFrame.N;j++) {
                        MapPoint *pMP2 = CurrentFrame.mvpMapPoints[j];
                        if (pMP == pMP2) {
                            featureMatched = true;
                            break;
                        }
                    }
                    if(featureMatched)
                        continue;

                    // Project it onto current and last frame to check if depth is positive
                    Eigen::Vector3d featureInGlobal = pMP->GetWorldPosEigen();
                    Eigen::Vector3d featureInLast = Taw.block<3,3>(0,0) * featureInGlobal + Taw.block<3,1>(0,3);
                    Eigen::Vector3d featureInCurrent = Tbw.block<3,3>(0,0) * featureInGlobal + Tbw.block<3,1>(0,3);
                    if (featureInCurrent(3) < 0 || featureInLast(3) < 0)
                        continue;

                    /// Information about last frame needed for tracking
                    // Camera matrix
                    Eigen::Matrix3d Ka = getCameraMatrix(LastFrame.fx, LastFrame.fy, LastFrame.cx, LastFrame.cy);

                    // Point to track
                    cv::KeyPoint kp = LastFrame.mvKeysUn[i];

                    // Getting the octave in last frame
                    int pyramidIndex = kp.octave;

                    // Getting the image pyramid for tracking
                    g2o::imgStr *lastImage = LastFrame.mpORBextractorLeft->photobaImagePyramid[pyramidIndex];

                    // Perform tracking
                    if ( trackMapPoint(pMP, CurrentFrame, featureInLast, Tba, Ka, lastImage, kp) ) {
                        nmatches++;

                        pMP->fForRescue = &LastFrame;
                        pMP->featureIndexForRescue = i;

                    }
                }
            }
        }



        return nmatches;
    }

    int PhotoTracker::SearchByPhoto(Frame &CurrentFrame, const vector<MapPoint*> &vpMapPoints) {
        // Number of tracked features
        int nmatches=0;

        // Current frame position
        Eigen::Matrix4d Taw, Tbw = CurrentFrame.getPose();

        // For all neighbouring features
        for(size_t iMP=0; iMP<vpMapPoints.size(); iMP++) {
            MapPoint *pMP = vpMapPoints[iMP];

            // Feature was not projected onto current frame
            if (!pMP->mbTrackInView)
                continue;

            if (pMP->isBad())
                continue;

            // Let's not consider already matched or tracked features in previous VO step
            bool found = false;
            for (int j = 0; j < CurrentFrame.N; j++) {
                if (CurrentFrame.mvpMapPoints[j] == pMP) {
                    found = true;
                    break;
                }
            }
            if (found) {
                continue;
            }


            // Project it onto current and last frame to check if depth is positive
            Eigen::Vector3d featureInGlobal = pMP->GetWorldPosEigen();
            Eigen::Vector3d featureInLast, featureInCurrent = Tbw.block<3,3>(0,0) * featureInGlobal + Tbw.block<3,1>(0,3);
            if (featureInCurrent(3) < 0)
                 continue;

            // Selecting the frame for tracking - one with the closest viewing angle that still contains image pyramid
            std::map<KeyFrame*,size_t> observations = pMP->GetObservations();

            KeyFrame* pKF;
            int pointInKFIndex = 0;
            double bestAngle = 2;
            for(map<KeyFrame*,size_t>::iterator mit=observations.begin(), mend=observations.end(); mit!=mend; mit++)
            {
                KeyFrame* tmpKF = mit->first;

                if (tmpKF->imagePyramidLeft.size() > 0)
                {
                    // Last frame position
                    Taw = tmpKF->GetPoseEigen();

                    // Feature in last frame
                    featureInLast = Taw.block<3,3>(0,0) * featureInGlobal + Taw.block<3,1>(0,3);

                    // Is depth positive
                    if (featureInLast(3) < 0)
                        continue;

                    // Computing the difference between observation angles
                    double diffCos = 1 - featureInLast.dot(featureInCurrent) / featureInLast.norm() / featureInCurrent.norm();

                    // We found a better KF for photo tracking
                    if (bestAngle > diffCos)
                    {
                        bestAngle = diffCos;
                        pKF = tmpKF;
                        pointInKFIndex = mit->second;
                    }
                }
            }
            if (!pKF)
                continue;

            // Relative transformation between poses
            Eigen::Matrix4d Tba = Tbw * Taw.inverse();

            /// Information about last frame needed for tracking
            // Camera matrix
            Eigen::Matrix3d Ka = getCameraMatrix(pKF->fx, pKF->fy, pKF->cx, pKF->cy);

            // Point to track
            cv::KeyPoint kp = pKF->mvKeysUn[pointInKFIndex];

            // Getting the octave in last frame
            int pyramidIndex = kp.octave;

            // Getting the image pyramid for tracking
            g2o::imgStr *lastImage = pKF->imagePyramidLeft[pyramidIndex];

            // Perform tracking
            if ( trackMapPoint(pMP, CurrentFrame, featureInLast, Tba, Ka, lastImage, kp) ) {
                nmatches++;
                pMP->kfForRescue = pKF;
                pMP->featureIndexForRescue = pointInKFIndex;
            }

        }

        return nmatches;
    }

    bool PhotoTracker::trackMapPoint(MapPoint *pMP, Frame &CurrentFrame,
            Eigen::Vector3d featureInLast, Eigen::Matrix4d Tba, Eigen::Matrix3d Ka,
            g2o::imgStr *lastImage, cv::KeyPoint kp) {

        // The patch normal in ref observation is assumed to be [0,0,-1] TODO: We could use estimated normal
        Eigen::Vector3d n(0,0,-1);

        // Computation of the distance to patch plane in image A
        double d = getDistanceToPlane(featureInLast, n);

        // Computation of the homography between A and B
        Eigen::Matrix3d Kb = getCameraMatrix(CurrentFrame.fx, CurrentFrame.fy, CurrentFrame.cx, CurrentFrame.cy);
        Eigen::Matrix3d H = computeHomography(Tba, n, d, Ka, Kb);

        //        const int pyramidIndex = pMP->mnTrackScaleLevel; // TODO: Predicted instead of octave of detection?
        int pyramidIndex = kp.octave;

        // Getting the image pyramids
        g2o::imgStr *currentImage = CurrentFrame.mpORBextractorLeft->photobaImagePyramid[pyramidIndex];

        // Getting the scale of the selected lvl
        const float pyramidScaleA = lastImage->imageScale;
        const float pyramidScaleB = currentImage->imageScale;

        // For all points that are considered neighbours (can be patch)
        std::vector<double> refPatch, curPatch;
        refPatch.reserve(neighbours.size());
        curPatch.reserve(neighbours.size());

        for (int i = 0; i < neighbours.size(); i++) {

            // Getting the patch value in last frame
            double refU = kp.pt.x / pyramidScaleA + neighbours[i].first;
            double refV = kp.pt.y / pyramidScaleA + neighbours[i].second;
            double refValue = getSubpixImageValue(refU, refV, lastImage->image);
            refPatch.push_back(refValue);

            // Projecting (x,y) from image A into image B with H
            Eigen::Vector3d pInB = H * Eigen::Vector3d(refU * pyramidScaleA, refV * pyramidScaleA, 1);
            pInB = pInB / pInB(2);

            // Getting the patch value in current frame
            double obsU = pInB(0) / pyramidScaleB;
            double obsV = pInB(1) / pyramidScaleB;
            double obsValue = getSubpixImageValue(obsU, obsV, currentImage->image);
            curPatch.push_back(obsValue);

            // Either of values is outside of the image
            if (refValue < 0 || obsValue < 0) {
                return false;
            }
        }

        //double errorWithAffine = computePatchDiffAffine(refPatch, curPatch);
        double errorAvg = computePatchDiffAvg(refPatch, curPatch);

        // If RMSE is lower than threshold then we are successful
        if (errorAvg < photoThreshold) {

            //std::cout << "[" << i << "] -> Matched: " << featureMatched << " | RMSE_Avg = " << errorAvg << std::endl;

            // Location in current frame
            Eigen::Vector3d current = H * Eigen::Vector3d(kp.pt.x, kp.pt.y, 1);
            current = current / current(2);

            addTrackedMapPoint(CurrentFrame, pMP, kp, current(0), current(1));
            return true;
        }
        return false;
    }


    void PhotoTracker::addTrackedMapPoint(Frame &CurrentFrame, MapPoint *pMP, cv::KeyPoint kp, double currentU, double currentV) {

        // Informing about the state
        pMP->rescuedLast = true;
        pMP->rescuedAtLeastOnce = true;

        // Artificially increasing the number of points detected in this case
        CurrentFrame.N++;

        // Adding map point in the end and saying it is an inlier
        CurrentFrame.mvpMapPoints.push_back(pMP);
        CurrentFrame.mvbOutlier.push_back(false);

        // Simulating a new position of the keypoint in the image
        kp.pt.x = currentU;
        kp.pt.y = currentV;
        CurrentFrame.mvKeys.push_back(kp);
        CurrentFrame.mvKeysUn.push_back(kp);

        CurrentFrame.mvuRight.push_back(-1);
        CurrentFrame.mvDepth.push_back(-1);

        // Empty descriptor
        cv::Mat emptyDescriptor = cv::Mat::zeros(1, CurrentFrame.mDescriptors.cols, CurrentFrame.mDescriptors.type());
        CurrentFrame.mDescriptors.push_back(emptyDescriptor);
    }
}