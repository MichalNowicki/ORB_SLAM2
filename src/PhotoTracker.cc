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

    cv::Mat PhotoTracker::getInversePose(cv::Mat Tcw) {
        cv::Mat Rcw = Tcw.rowRange(0, 3).colRange(0, 3);
        cv::Mat tcw = Tcw.rowRange(0, 3).col(3);
        cv::Mat Rwc = Rcw.t();
        cv::Mat Ow = -Rwc * tcw;

        cv::Mat Twc = cv::Mat::eye(4, 4, Tcw.type());
        Rwc.copyTo(Twc.rowRange(0, 3).colRange(0, 3));
        Ow.copyTo(Twc.rowRange(0, 3).col(3));
        return Twc;
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




    int PhotoTracker::SearchByPhoto(Frame &CurrentFrame, const Frame &LastFrame) {

        // Current frame position
        const cv::Mat Rcw = CurrentFrame.mTcw.rowRange(0, 3).colRange(0, 3);
        const cv::Mat tcw = CurrentFrame.mTcw.rowRange(0, 3).col(3);

        // Last frame position
        const cv::Mat Rlw = LastFrame.mTcw.rowRange(0, 3).colRange(0, 3);
        const cv::Mat tlw = LastFrame.mTcw.rowRange(0, 3).col(3);

        cv::Mat Twa = getInversePose(LastFrame.mTcw);
        cv::Mat Tbw = CurrentFrame.mTcw;
        Eigen::Matrix4d Twa_eig, Tbw_eig;
        cv::cv2eigen(Twa, Twa_eig);
        cv::cv2eigen(Tbw, Tbw_eig);

        Eigen::Matrix4d Tba = Tbw_eig * Twa_eig;

        int nmatches = 0;

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
                    cv::Mat x3Dw = pMP->GetWorldPos();
                    cv::Mat x3Dc = Rcw * x3Dw + tcw;
                    cv::Mat x3Dl = Rlw * x3Dw + tlw;

                    Eigen::Vector3d featureInCurrent, featureInLast;
                    cv::cv2eigen(x3Dc, featureInCurrent);
                    cv::cv2eigen(x3Dl, featureInLast);
                    if (featureInCurrent(3) < 0 || featureInLast(3) < 0)
                        continue;

                    // The patch normal in ref observation is assumed to be [0,0,-1]
                    Eigen::Vector3d n(0,0,-1);

                    // Computation of the distance to patch plane in image A
                    double d = getDistanceToPlane(featureInLast, n);

                    // Computation of the homography between A and B
                    Eigen::Matrix3d Ka = getCameraMatrix(LastFrame.fx, LastFrame.fy, LastFrame.cx, LastFrame.cy);
                    Eigen::Matrix3d Kb = getCameraMatrix(CurrentFrame.fx, CurrentFrame.fy, CurrentFrame.cx, CurrentFrame.cy);
                    Eigen::Matrix3d H = computeHomography(Tba, n, d, Ka, Kb);

                    // Location in last frame
                    float lastU = LastFrame.mvKeysUn[i].pt.x;
                    float lastV = LastFrame.mvKeysUn[i].pt.y;

                    // Location in current frame
                    Eigen::Vector3d last(lastU, lastV, 1);
                    Eigen::Vector3d current = H * last;
                    current = current / current(2);
                    float currentU = current(0), currentV = current(1);

                    // Getting the octave in last frame
                    int pyramidIndex = LastFrame.mvKeysUn[i].octave;

                    // Getting the image pyramids
                    g2o::imgStr *lastImage = LastFrame.mpORBextractorLeft->photobaImagePyramid[pyramidIndex];
                    g2o::imgStr *currentImage = CurrentFrame.mpORBextractorLeft->photobaImagePyramid[pyramidIndex];

                    // Getting the scale of the selected lvl
                    const float pyramidScale = currentImage->imageScale;

                    // For all points that are considered neighbours (can be patch)
                    double error = 0;
                    bool success = true;
                    std::vector<double> refPatch, curPatch;
                    refPatch.reserve(neighbours.size());
                    curPatch.reserve(neighbours.size());

                    for (int i = 0; i < neighbours.size(); i++) {

                        // Getting the patch value in last frame
                        double refU = lastU / pyramidScale + neighbours[i].first;
                        double refV = lastV / pyramidScale + neighbours[i].second;
                        double refValue = getSubpixImageValue(refU, refV, lastImage->image);
                        refPatch.push_back(refValue);

                        // Projecting (x,y) from image A into image B with  H
                        Eigen::Vector3d pInA(refU * pyramidScale, refV * pyramidScale, 1);
                        Eigen::Vector3d pInB = H * pInA;
                        pInB = pInB / pInB(2);

                        // Getting the patch value in current frame
                        double obsU = pInB(0) / pyramidScale;
                        double obsV = pInB(1) / pyramidScale;
                        double obsValue = getSubpixImageValue(obsU, obsV, currentImage->image);
                        curPatch.push_back(obsValue);

                        // std::cout << "(" << refU << "," << refV <<") = " << refValue << "\t (" << obsU << "," << obsV<<") = " << obsValue << std::endl;

                        // Either of values is outside of the image
                        if (refValue < 0 || obsValue < 0) {
                            success = false;
                            break;
                        }
                    }

                    if (success) {
                        //double errorWithAffine = computePatchDiffAffine(refPatch, curPatch);
                        double errorAvg = computePatchDiffAvg(refPatch, curPatch);


                        // If RMSE is lower than threshold then we are succesful
                        if (errorAvg < photoThreshold) {

                            //std::cout << "[" << i << "] -> Matched: " << featureMatched << " | RMSE_Avg = " << errorAvg << std::endl;

                            // TODO: How to only mark for future?
                            pMP->rescuedLast = true;
                            pMP->rescuedAtLeastOnce = true;
                            CurrentFrame.N++;
                            CurrentFrame.mvpMapPoints.push_back(pMP);
                            CurrentFrame.mvbOutlier.push_back(false);

                            cv::KeyPoint kp = LastFrame.mvKeys[i];
                            kp.pt.x = currentU;
                            kp.pt.y = currentV;
                            CurrentFrame.mvKeys.push_back(kp);

                            kp = LastFrame.mvKeysUn[i];
                            kp.pt.x = currentU;
                            kp.pt.y = currentV;
                            CurrentFrame.mvKeysUn.push_back(kp);


                            CurrentFrame.mvuRight.push_back(-1);
                            CurrentFrame.mvDepth.push_back(-1);

                            cv::Mat emptyDescriptor = cv::Mat::zeros(1, CurrentFrame.mDescriptors.cols,
                                                                     CurrentFrame.mDescriptors.type());
                            CurrentFrame.mDescriptors.push_back(emptyDescriptor);

                            // Simulate descriptor ?

                            nmatches++;


                        }
                    }
                }
            }
        }



        return nmatches;
    }

    int PhotoTracker::SearchByPhoto(Frame &CurrentFrame, const vector<MapPoint*> &vpMapPoints) {
        int nmatches=0;

        // Current frame position
        const cv::Mat Rcw = CurrentFrame.mTcw.rowRange(0, 3).colRange(0, 3);
        const cv::Mat tcw = CurrentFrame.mTcw.rowRange(0, 3).col(3);
        cv::Mat Tbw = CurrentFrame.mTcw;

        // For all neighbouring features
        for(size_t iMP=0; iMP<vpMapPoints.size(); iMP++) {
            MapPoint *pMP = vpMapPoints[iMP];
            if (!pMP->mbTrackInView)
                continue;

            if (pMP->isBad())
                continue;

            // Let's not consider already matched or tracked features
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

            // Selecting the frame for tracking
            // TODO: We probably should select the KF with the closest viewing angle that still contains necessary image for optimization
            std::map<KeyFrame*,size_t> observations = pMP->GetObservations();

            KeyFrame* pKF;
            int pKFIndex = 0;
            for(map<KeyFrame*,size_t>::iterator mit=observations.begin(), mend=observations.end(); mit!=mend; mit++)
            {
                pKF = mit->first;
                pKFIndex = mit->second;
                break;
            }

            // Last frame position
            const cv::Mat Rlw = pKF->GetPose().rowRange(0, 3).colRange(0, 3);
            const cv::Mat tlw = pKF->GetPose().rowRange(0, 3).col(3);

            cv::Mat Twa = getInversePose(pKF->GetPose());
            Eigen::Matrix4d Twa_eig, Tbw_eig;
            cv::cv2eigen(Twa, Twa_eig);
            cv::cv2eigen(Tbw, Tbw_eig);

            Eigen::Matrix4d Tba = Tbw_eig * Twa_eig;

            // Project it onto current and last frame to check if depth is positive
            cv::Mat x3Dw = pMP->GetWorldPos();
            cv::Mat x3Dc = Rcw * x3Dw + tcw;
            cv::Mat x3Dl = Rlw * x3Dw + tlw;

            // Getting the 3D feature location in last and current frame
            Eigen::Vector3d featureInCurrent, featureInLast;
            cv::cv2eigen(x3Dc, featureInCurrent);
            cv::cv2eigen(x3Dl, featureInLast);
            if (featureInCurrent(3) < 0 || featureInLast(3) < 0)
                continue;

            // The patch normal in ref observation is assumed to be [0,0,-1]
            Eigen::Vector3d n(0,0,-1);

            // Computation of the distance to patch plane in image A
            double d = getDistanceToPlane(featureInLast, n);

            // Computation of the homography between A and B
            Eigen::Matrix3d Ka = getCameraMatrix(pKF->fx, pKF->fy, pKF->cx, pKF->cy);
            Eigen::Matrix3d Kb = getCameraMatrix(CurrentFrame.fx, CurrentFrame.fy, CurrentFrame.cx, CurrentFrame.cy);
            Eigen::Matrix3d H = computeHomography(Tba, n, d, Ka, Kb);

            // Location in last frame
            int lastU = pKF->mvKeysUn[pKFIndex].pt.x;
            int lastV = pKF->mvKeysUn[pKFIndex].pt.y;

            // Location in current frame
            Eigen::Vector3d last(lastU, lastV, 1);
            Eigen::Vector3d current = H * last;
            current = current / current(2);
            float currentU = current(0), currentV = current(1);

            // Getting the octave in last frame
//            const int &nPredictedLevel = pMP->mnTrackScaleLevel;
            int pyramidIndex = pKF->mvKeysUn[pKFIndex].octave;

            // Getting the image pyramids
            g2o::imgStr *lastImage = pKF->imagePyramidLeft[pyramidIndex];
            g2o::imgStr *currentImage = CurrentFrame.mpORBextractorLeft->photobaImagePyramid[pyramidIndex];

            // Getting the scale of the selected lvl
            const float pyramidScale = currentImage->imageScale;

            // For all points that are considered neighbours (can be patch)
            bool success = true;
            std::vector<double> refPatch, curPatch;
            refPatch.reserve(neighbours.size());
            curPatch.reserve(neighbours.size());

            for (int i = 0; i < neighbours.size(); i++) {

                // Getting the patch value in last frame
                double refU = lastU / pyramidScale + neighbours[i].first;
                double refV = lastV / pyramidScale + neighbours[i].second;
                double refValue = getSubpixImageValue(refU, refV, lastImage->image);
                refPatch.push_back(refValue);

                // Projecting (x,y) from image A into image B with  H
                Eigen::Vector3d pInA(refU * pyramidScale, refV * pyramidScale, 1);
                Eigen::Vector3d pInB = H * pInA;
                pInB = pInB / pInB(2);

                // Getting the patch value in current frame
                double obsU = pInB(0) / pyramidScale;
                double obsV = pInB(1) / pyramidScale;
                double obsValue = getSubpixImageValue(obsU, obsV, currentImage->image);
                curPatch.push_back(obsValue);

                // std::cout << "(" << refU << "," << refV <<") = " << refValue << "\t (" << obsU << "," << obsV<<") = " << obsValue << std::endl;

                // Either of values is outside of the image
                if (refValue < 0 || obsValue < 0) {
                    success = false;
                    break;
                }
            }

            if (success) {
                //double errorWithAffine = computePatchDiffAffine(refPatch, curPatch);
                double errorAvg = computePatchDiffAvg(refPatch, curPatch);

                // If RMSE is lower than threshold then we are successful
                if (errorAvg < photoThreshold) {

                    //std::cout << "[" << i << "] -> Matched: " << featureMatched << " | RMSE_Avg = " << errorAvg << std::endl;

                    pMP->rescuedLast = true;
                    pMP->rescuedAtLeastOnce = true;
                    CurrentFrame.N++;
                    CurrentFrame.mvpMapPoints.push_back(pMP);
                    CurrentFrame.mvbOutlier.push_back(false);

                    cv::KeyPoint kp = pKF->mvKeys[pKFIndex];
                    kp.pt.x = currentU;
                    kp.pt.y = currentV;
                    CurrentFrame.mvKeys.push_back(kp);

                    kp = pKF->mvKeysUn[pKFIndex];
                    kp.pt.x = currentU;
                    kp.pt.y = currentV;
                    CurrentFrame.mvKeysUn.push_back(kp);


                    CurrentFrame.mvuRight.push_back(-1);
                    CurrentFrame.mvDepth.push_back(-1);

                    // Empty descriptor
                    cv::Mat emptyDescriptor = cv::Mat::zeros(1, CurrentFrame.mDescriptors.cols,
                                                             CurrentFrame.mDescriptors.type());
                    CurrentFrame.mDescriptors.push_back(emptyDescriptor);


                    nmatches++;
                }
            }
        }

        return nmatches;
    }

}