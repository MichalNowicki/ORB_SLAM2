#include"PhotoTracker.h"

namespace ORB_SLAM2 {


    PhotoTracker::PhotoTracker(double photoThreshold, int kltMaxIterations, double kltEPS, double kltZNCCThr, int kltPatchSize,
            bool verbose, double kltMaxMovement) :
        photoThreshold(photoThreshold), kltMaxIterations(kltMaxIterations), kltEPS(kltEPS), kltZNCCThr(kltZNCCThr), patchSize(kltPatchSize),
        verbose(verbose), kltMaxMovement(kltMaxMovement) {

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

                    // Project it onto current and last frame to check if depth is positive
                    Eigen::Vector3d featureInGlobal = pMP->GetWorldPosEigen();
                    Eigen::Vector3d featureInLast = Taw.block<3,3>(0,0) * featureInGlobal + Taw.block<3,1>(0,3);
                    Eigen::Vector3d featureInCurrent = Tbw.block<3,3>(0,0) * featureInGlobal + Tbw.block<3,1>(0,3);
                    if (featureInCurrent(3) < 0 || featureInLast(3) < 0)
                        continue;

                    /// Information about last frame needed for tracking
                    // Camera matrix
                    Eigen::Matrix3d Ka = photo::getCameraMatrix(LastFrame.fx, LastFrame.fy, LastFrame.cx, LastFrame.cy);

                    // Point to track
                    cv::KeyPoint kp = LastFrame.mvKeysUn[i];


                    // TODO: What lvl should I use?
                    // Getting the octave in last frame
//                    int pyramidIndex = kp.octave;
                    int pyramidIndex = 2;

                    // Getting the image pyramid for tracking
                    photo::imgStr *lastImage = LastFrame.imagePyramidLeft[pyramidIndex];

                    // Perform tracking
                    if ( trackMapPoint(pMP, CurrentFrame, featureInLast, Tba, Ka, lastImage, kp) ) {
                        nmatches++;

                        pMP->fForRescue = &LastFrame;
                        pMP->kfForRescue = static_cast<KeyFrame*>(NULL);
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

            // Let's not consider already matched & tracked features in previous VO step
//            if (pMP->matchedLast || pMP->rescuedLast)
//                break;

            // Project it onto current and last frame to check if depth is positive
            Eigen::Vector3d featureInGlobal = pMP->GetWorldPosEigen();
            Eigen::Vector3d featureInLast, featureInCurrent = Tbw.block<3,3>(0,0) * featureInGlobal + Tbw.block<3,1>(0,3);
            if (featureInCurrent(3) < 0)
                 continue;

            // Selecting the frame for tracking - one with the closest viewing angle that still contains image pyramid
            std::map<KeyFrame*,size_t> observations = pMP->GetObservations();

            KeyFrame* pKF = static_cast<KeyFrame*>(NULL);
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
            Eigen::Matrix3d Ka = photo::getCameraMatrix(pKF->fx, pKF->fy, pKF->cx, pKF->cy);

            // Point to track
            cv::KeyPoint kp = pKF->mvKeysUn[pointInKFIndex];

            // Getting the octave in last frame
            // TODO: What lvl should I use?
//            int pyramidIndex = kp.octave;
            int pyramidIndex = 2;

            // Getting the image pyramid for tracking
            photo::imgStr *lastImage = pKF->imagePyramidLeft[pyramidIndex];

            // Perform tracking
            if ( trackMapPoint(pMP, CurrentFrame, featureInLast, Tba, Ka, lastImage, kp) ) {
                nmatches++;
                pMP->fForRescue = static_cast<Frame*>(NULL);
                pMP->kfForRescue = pKF;
                pMP->featureIndexForRescue = pointInKFIndex;
            }

        }

        return nmatches;
    }

    photoTrackerResult PhotoTracker::SearchByKLT(Frame &CurrentFrame, Frame &LastFrame) {
        // Number of tracked features
        int nmatches = 0;

        // Relative transformation between frames in Eigen
        Eigen::Matrix4d Taw = LastFrame.getPose();
        Eigen::Matrix4d Tbw = CurrentFrame.getPose();
        Eigen::Matrix4d Tba = Tbw * Taw.inverse();

        // For KLT
        std::vector<int> indices;
        std::vector<cv::Point2f> prevPts, nextPts;

        // For all map points observed in last frame
        for (int i = 0; i < LastFrame.N; i++) {
            MapPoint *pMP = LastFrame.mvpMapPoints[i];

            if (pMP) {
                if (!LastFrame.mvbOutlier[i]) {

                    // Do only tracking if matching failed
                    // TODO: Let's give preference for tracking
                    if (!pMP->matchedLast) {

                        // Project it onto current and last frame to check if depth is positive
                        Eigen::Vector3d featureInGlobal = pMP->GetWorldPosEigen();
                        Eigen::Vector3d featureInLast = Taw.block<3, 3>(0, 0) * featureInGlobal + Taw.block<3, 1>(0, 3);
                        Eigen::Vector3d featureInCurrent =
                                Tbw.block<3, 3>(0, 0) * featureInGlobal + Tbw.block<3, 1>(0, 3);
                        if (featureInCurrent(3) < 0 || featureInLast(3) < 0)
                            continue;

                        /// Information about last frame needed for tracking
                        // Camera matrix
                        Eigen::Matrix3d Ka = photo::getCameraMatrix(LastFrame.fx, LastFrame.fy, LastFrame.cx,
                                                                    LastFrame.cy);

                        Eigen::Vector3d projectedInCurrent = Ka * featureInCurrent;
                        projectedInCurrent = projectedInCurrent / projectedInCurrent(2);

                        // Projection onto the current image
                        cv::Point2f point;
                        point.x = projectedInCurrent(0);
                        point.y = projectedInCurrent(1);

                        if (point.x < CurrentFrame.mnMinX || point.x > CurrentFrame.mnMaxX)
                            continue;
                        if (point.y < CurrentFrame.mnMinY || point.y > CurrentFrame.mnMaxY)
                            continue;

                        // Point to track in previous and in current
                        cv::KeyPoint kp = LastFrame.mvKeysUn[i];
                        //std::cout << "Comparison : " << kp.pt.x << " " << kp.pt.y << " --- " << projectedInCurrent(0) << " " << projectedInCurrent(1) << std::endl;

                        prevPts.push_back(kp.pt);
                        nextPts.push_back(point);
                        indices.push_back(i);
                    }
                }
            }
        }

        std::vector<cv::Point2f> motionGuessPts = nextPts;

        /*  KLT TRACKING
         *      Size 	winSize = Size(21, 21),
                int 	maxLevel = 3,
                TermCriteria 	criteria = TermCriteria(TermCriteria::COUNT+TermCriteria::EPS, 30, 0.01),
                int 	flags = 0,
                double 	minEigThreshold = 1e-4
         */
        vector<uchar> status;
        vector<float> err;
        cv::calcOpticalFlowPyrLK(LastFrame.origImgPyramid, CurrentFrame.origImgPyramid,
                                 prevPts, nextPts, status, err, cv::Size(patchSize,patchSize), 3,
                                 cv::TermCriteria(cv::TermCriteria::COUNT+cv::TermCriteria::EPS, kltMaxIterations, kltEPS), cv::OPTFLOW_USE_INITIAL_FLOW);

//        cv::imwrite("logs/last.png", LastFrame.origImgPyramid[0]);
//        cv::imwrite("logs/current.png", CurrentFrame.origImgPyramid[0]);

        int belowThCount = 0, extra = 0;
        double avgTravel = 0;
        for(int i=0;i<status.size();i++) {
            if (status[i]) {

                nmatches++;

                double zncc = computeZNCC(LastFrame, CurrentFrame, prevPts[i], nextPts[i], patchSize);
                double travelDistance = cv::norm(motionGuessPts[i]-nextPts[i]);

                if (zncc > kltZNCCThr && travelDistance < kltMaxMovement) {
                    belowThCount++;
                    avgTravel += travelDistance;


                    int index = indices[i];

                    MapPoint *pMP = LastFrame.mvpMapPoints[index];

                    // Add artificial feature if it was not matched with descriptors
//                     TODO: Let's also track if matched - either one or the other due to map for observations
                    if (!pMP->matchedLast) {

                        // Informing about the state
                        pMP->rescuedLast = true;
                        pMP->rescuedAtLeastOnce = true;
                        pMP->lastZNCC = zncc;

                        pMP->fForRescue = &LastFrame;
                        pMP->kfForRescue = static_cast<KeyFrame *>(NULL);
                        pMP->featureIndexForRescue = index;

                        cv::KeyPoint kp = LastFrame.mvKeysUn[index];
                        addTrackedMapPoint(CurrentFrame, pMP, kp, nextPts[i].x, nextPts[i].y);

                        // TODO: For now duplicating the last descriptor
                        CurrentFrame.mDescriptors.row(CurrentFrame.mDescriptors.rows - 1) = LastFrame.mDescriptors.row(
                                index);

                        extra++;
                    }
                }
            }
        }

        if (verbose)
            std::cout << "Tracked: " << nmatches << " below thr: " << belowThCount << " out of " << status.size() <<
                " | Extra: " << extra << " Avg travel dist: " << avgTravel/belowThCount << std::endl;

        photoTrackerResult result;
        result.all = status.size();
        result.tracked = nmatches;
        result.trackedBelowTh = belowThCount;
        result.extraOverMatchings = extra;
        result.avgTravelDistForInliers = avgTravel/belowThCount;

//        exit(0);
        return result;
    }

    int PhotoTracker::SearchByKLT(Frame &CurrentFrame, const vector<MapPoint*> &vpMapPoints) {
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

            // Let's not consider already matched & tracked features in previous VO step
            if (pMP->matchedLast || pMP->rescuedLast)
                break;

            // Project it onto current and last frame to check if depth is positive
            Eigen::Vector3d featureInGlobal = pMP->GetWorldPosEigen();
            Eigen::Vector3d featureInLast, featureInCurrent = Tbw.block<3,3>(0,0) * featureInGlobal + Tbw.block<3,1>(0,3);
            if (featureInCurrent(3) < 0)
                continue;

            // Selecting the frame for tracking - one with the closest viewing angle that still contains image
            std::map<KeyFrame*,size_t> observations = pMP->GetObservations();

            KeyFrame* pKF = static_cast<KeyFrame*>(NULL);
            int pointInKFIndex = 0;
            double bestAngle = 2;
            for(map<KeyFrame*,size_t>::iterator mit=observations.begin(), mend=observations.end(); mit!=mend; mit++)
            {
                KeyFrame* tmpKF = mit->first;

                if (!tmpKF->origImgPyramid.empty())
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

            /// Information about last frame needed for tracking
            // Camera matrix
            Eigen::Matrix3d Ka = photo::getCameraMatrix(CurrentFrame.fx, CurrentFrame.fy, CurrentFrame.cx, CurrentFrame.cy);

            Eigen::Vector3d projectedInCurrent = Ka * featureInCurrent;
            projectedInCurrent = projectedInCurrent / projectedInCurrent(2);

            // Projection onto the current image
            cv::Point2f point;
            point.x = projectedInCurrent(0);
            point.y = projectedInCurrent(1);

            if(point.x<CurrentFrame.mnMinX || point.x>CurrentFrame.mnMaxX)
                continue;
            if(point.y<CurrentFrame.mnMinY || point.y>CurrentFrame.mnMaxY)
                continue;


            // Point to track
            cv::KeyPoint kp = pKF->mvKeysUn[pointInKFIndex];


            /*  KLT TRACKING
                    Size 	winSize = Size(21, 21),
                    int 	maxLevel = 3,
                    TermCriteria 	criteria = TermCriteria(TermCriteria::COUNT+TermCriteria::EPS, 30, 0.01),
                    int 	flags = 0,
                    double 	minEigThreshold = 1e-4
            */
            vector<uchar> status;
            vector<float> err;
            vector<cv::Point2f> prevPts, nextPts;
            prevPts.push_back(kp.pt);
            nextPts.push_back(point);

            cv::calcOpticalFlowPyrLK(pKF->origImgPyramid, CurrentFrame.origImgPyramid,
                                     prevPts, nextPts, status, err, cv::Size(9,9), 3,
                                     cv::TermCriteria(cv::TermCriteria::COUNT+cv::TermCriteria::EPS, kltMaxIterations, kltEPS), cv::OPTFLOW_USE_INITIAL_FLOW);


            // TODO: error is not properly defined with ZNCC
            if (status[0]) // & ZNCC
            {
                nmatches++;

                // Informing about the state
                pMP->rescuedLast = true;
                pMP->rescuedAtLeastOnce = true;
                pMP->fForRescue = static_cast<Frame*>(NULL);
                pMP->kfForRescue = pKF;
                pMP->featureIndexForRescue = pointInKFIndex;

                addTrackedMapPoint(CurrentFrame, pMP, kp, nextPts[0].x, nextPts[0].y);

                // TODO: For now duplicating the last descriptor
                CurrentFrame.mDescriptors.row(CurrentFrame.mDescriptors.rows-1) = pKF->mDescriptors.row(pointInKFIndex);
            }
        }

        return nmatches;
    }

    bool PhotoTracker::trackMapPoint(MapPoint *pMP, Frame &CurrentFrame,
            Eigen::Vector3d featureInLast, Eigen::Matrix4d Tba, Eigen::Matrix3d Ka,
            photo::imgStr *lastImage, cv::KeyPoint kp) {

        // The patch normal in ref observation is assumed to be [0,0,-1] TODO: We could use estimated normal
        Eigen::Vector3d n(0,0,-1);

        // Computation of the distance to patch plane in image A
        double d = photo::getDistanceToPlane(featureInLast, n);

        // Computation of the homography between A and B
        Eigen::Matrix3d Kb = photo::getCameraMatrix(CurrentFrame.fx, CurrentFrame.fy, CurrentFrame.cx, CurrentFrame.cy);
        Eigen::Matrix3d H = photo::computeHomography(Tba, n, d, Ka, Kb);

        // TODO: What lvl should I use?
        //        const int pyramidIndex = pMP->mnTrackScaleLevel; // TODO: Predicted instead of octave of detection?
//        int pyramidIndex = kp.octave;
        int pyramidIndex = 2;

        // Getting the image pyramids
        photo::imgStr *currentImage = CurrentFrame.mpORBextractorLeft->photobaImagePyramid[pyramidIndex];

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
            double refValue = photo::getSubpixImageValue(refU, refV, lastImage->image);
            refPatch.push_back(refValue);

            // Projecting (x,y) from image A into image B with H
            Eigen::Vector3d pInB = H * Eigen::Vector3d(refU * pyramidScaleA, refV * pyramidScaleA, 1);
            pInB = pInB / pInB(2);

            // Getting the patch value in current frame
            double obsU = pInB(0) / pyramidScaleB;
            double obsV = pInB(1) / pyramidScaleB;
            double obsValue = photo::getSubpixImageValue(obsU, obsV, currentImage->image);
            curPatch.push_back(obsValue);

            // Either of values is outside of the image
            if (refValue < 0 || obsValue < 0) {
                return false;
            }
        }

        //double errorWithAffine = computePatchDiffAffine(refPatch, curPatch);
        double errorAvg = photo::computePatchDiffAvg(refPatch, curPatch);

        // If RMSE is lower than threshold then we are successful
        if (errorAvg < photoThreshold) {

            //std::cout << "RMSE_Avg = " << errorAvg << std::endl;

            // Location in current frame
            Eigen::Vector3d current = H * Eigen::Vector3d(kp.pt.x, kp.pt.y, 1);
            current = current / current(2);

            // Informing about the state
            pMP->rescuedLast = true;
            pMP->rescuedAtLeastOnce = true;

            // Add artificial feature if it was not matched with descriptors
            if (!pMP->matchedLast)
                addTrackedMapPoint(CurrentFrame, pMP, kp, current(0), current(1));

            return true;
        }
        return false;
    }


    void PhotoTracker::addTrackedMapPoint(Frame &CurrentFrame, MapPoint *pMP, cv::KeyPoint kp, double currentU, double currentV) {
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

    double PhotoTracker::computeZNCC(Frame &LastFrame, Frame &CurrentFrame, cv::Point2f &lastPt, cv::Point2f &curPt, const int &patchSize) {
        const int halfPatchSize = patchSize/2;
        const int sqPatchSize = patchSize * patchSize;

        // ZNNC
        double avgOrig = 0, avgCur = 0;
        for (int xShift=-halfPatchSize;xShift<=halfPatchSize;xShift++) {
            for (int yShift=-halfPatchSize;yShift<=halfPatchSize;yShift++) {

                double origVal = photo::getSubpixImageValue(lastPt.x+xShift,lastPt.y+yShift,LastFrame.origImgPyramid[0]);
                double curVal = photo::getSubpixImageValue(curPt.x+xShift,curPt.y+yShift,CurrentFrame.origImgPyramid[0]);

                avgOrig += origVal;
                avgCur += curVal;
            }
        }

        avgOrig /= sqPatchSize;
        avgCur /= sqPatchSize;

        double nom = 0, sigmaOrig = 0, sigmaCur = 0;
        for (int xShift=-halfPatchSize;xShift<=halfPatchSize;xShift++) {
            for (int yShift=-halfPatchSize;yShift<=halfPatchSize;yShift++) {

                double origVal = photo::getSubpixImageValue(lastPt.x+xShift,lastPt.y+yShift,LastFrame.origImgPyramid[0]);
                double curVal = photo::getSubpixImageValue(curPt.x+xShift,curPt.y+yShift,CurrentFrame.origImgPyramid[0]);

                nom += (origVal - avgOrig)*(curVal - avgCur);
                sigmaOrig += pow(origVal - avgOrig, 2);
                sigmaCur += pow(curVal - avgCur, 2);
            }
        }
        return nom / sqrt(sigmaOrig * sigmaCur);
    }

}