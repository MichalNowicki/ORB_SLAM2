/**
* This file is part of ORB-SLAM2.
*
* Copyright (C) 2014-2016 Raúl Mur-Artal <raulmur at unizar dot es> (University of Zaragoza)
* For more information see <https://github.com/raulmur/ORB_SLAM2>
*
* ORB-SLAM2 is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* ORB-SLAM2 is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with ORB-SLAM2. If not, see <http://www.gnu.org/licenses/>.
*/

#include "MapPoint.h"
#include "ORBmatcher.h"
#include <Eigen/Eigen>

#include "Modifications/PatchRefinement.h"
#include<mutex>

namespace ORB_SLAM2
{

long unsigned int MapPoint::nNextId=0;
mutex MapPoint::mGlobalMutex;

MapPoint::MapPoint(const cv::Mat &Pos, KeyFrame *pRefKF, Map* pMap):
    mnFirstKFid(pRefKF->mnId), mnFirstFrame(pRefKF->mnFrameId), nObs(0), mnTrackReferenceForFrame(0),
    mnLastFrameSeen(0), mnBALocalForKF(0), mnFuseCandidateForKF(0), mnLoopPointForKF(0), mnCorrectedByKF(0),
    mnCorrectedReference(0), mnBAGlobalForKF(0), mpRefKF(pRefKF), mnVisible(1), mnFound(1), mbBad(false),
    mpReplaced(static_cast<MapPoint*>(NULL)), mfMinDistance(0), mfMaxDistance(0), mpMap(pMap)
{
    Pos.copyTo(mWorldPos);
    mNormalVector = cv::Mat::zeros(3,1,CV_32F);

    // MapPoints can be created from Tracking and Local Mapping. This mutex avoid conflicts with id.
    unique_lock<mutex> lock(mpMap->mMutexPointCreation);
    mnId=nNextId++;
}

MapPoint::MapPoint(const cv::Mat &Pos, Map* pMap, Frame* pFrame, const int &idxF):
    mnFirstKFid(-1), mnFirstFrame(pFrame->mnId), nObs(0), mnTrackReferenceForFrame(0), mnLastFrameSeen(0),
    mnBALocalForKF(0), mnFuseCandidateForKF(0),mnLoopPointForKF(0), mnCorrectedByKF(0),
    mnCorrectedReference(0), mnBAGlobalForKF(0), mpRefKF(static_cast<KeyFrame*>(NULL)), mnVisible(1),
    mnFound(1), mbBad(false), mpReplaced(NULL), mpMap(pMap)
{
    Pos.copyTo(mWorldPos);
    cv::Mat Ow = pFrame->GetCameraCenter();
    mNormalVector = mWorldPos - Ow;
    mNormalVector = mNormalVector/cv::norm(mNormalVector);

    cv::Mat PC = Pos - Ow;
    const float dist = cv::norm(PC);
    const int level = pFrame->mvKeysUn[idxF].octave;
    const float levelScaleFactor =  pFrame->mvScaleFactors[level];
    const int nLevels = pFrame->mnScaleLevels;

    mfMaxDistance = dist*levelScaleFactor;
    mfMinDistance = mfMaxDistance/pFrame->mvScaleFactors[nLevels-1];

    pFrame->mDescriptors.row(idxF).copyTo(mDescriptor);
    pFrame->mPatches[idxF].copyTo(mPatch);

    // MapPoints can be created from Tracking and Local Mapping. This mutex avoid conflicts with id.
    unique_lock<mutex> lock(mpMap->mMutexPointCreation);
    mnId=nNextId++;
}

void MapPoint::SetWorldPos(const cv::Mat &Pos)
{
    unique_lock<mutex> lock2(mGlobalMutex);
    unique_lock<mutex> lock(mMutexPos);
    Pos.copyTo(mWorldPos);
}

cv::Mat MapPoint::GetWorldPos()
{
    unique_lock<mutex> lock(mMutexPos);
    return mWorldPos.clone();
}

cv::Mat MapPoint::GetNormal()
{
    unique_lock<mutex> lock(mMutexPos);
    return mNormalVector.clone();
}

KeyFrame* MapPoint::GetReferenceKeyFrame()
{
    unique_lock<mutex> lock(mMutexFeatures);
    return mpRefKF;
}

void MapPoint::AddObservation(KeyFrame* pKF, size_t idx)
{
    unique_lock<mutex> lock(mMutexFeatures);
    if(mObservations.count(pKF))
        return;
    mObservations[pKF]=idx;

    if(pKF->mvuRight[idx]>=0)
        nObs+=2;
    else
        nObs++;
}

void MapPoint::EraseObservation(KeyFrame* pKF)
{
    bool bBad=false;
    {
        unique_lock<mutex> lock(mMutexFeatures);
        if(mObservations.count(pKF))
        {
            int idx = mObservations[pKF];
            if(pKF->mvuRight[idx]>=0)
                nObs-=2;
            else
                nObs--;

            mObservations.erase(pKF);

            if(mpRefKF==pKF) // TODO: When we change ref keyframe, it is necessary to recompute patch matchings
                mpRefKF=mObservations.begin()->first;

            // If only 2 observations or less, discard point
            if(nObs<=2)
                bBad=true;
        }
    }

    if(bBad)
        SetBadFlag();
}

map<KeyFrame*, size_t> MapPoint::GetObservations()
{
    unique_lock<mutex> lock(mMutexFeatures);
    return mObservations;
}

int MapPoint::Observations()
{
    unique_lock<mutex> lock(mMutexFeatures);
    return nObs;
}

void MapPoint::SetBadFlag()
{
    map<KeyFrame*,size_t> obs;
    {
        unique_lock<mutex> lock1(mMutexFeatures);
        unique_lock<mutex> lock2(mMutexPos);
        mbBad=true;
        obs = mObservations;
        mObservations.clear();
    }
    for(map<KeyFrame*,size_t>::iterator mit=obs.begin(), mend=obs.end(); mit!=mend; mit++)
    {
        KeyFrame* pKF = mit->first;
        pKF->EraseMapPointMatch(mit->second);
    }

    mpMap->EraseMapPoint(this);
}

MapPoint* MapPoint::GetReplaced()
{
    unique_lock<mutex> lock1(mMutexFeatures);
    unique_lock<mutex> lock2(mMutexPos);
    return mpReplaced;
}

void MapPoint::Replace(MapPoint* pMP)
{
    if(pMP->mnId==this->mnId)
        return;

    int nvisible, nfound;
    map<KeyFrame*,size_t> obs;
    {
        unique_lock<mutex> lock1(mMutexFeatures);
        unique_lock<mutex> lock2(mMutexPos);
        obs=mObservations;
        mObservations.clear();
        mbBad=true;
        nvisible = mnVisible;
        nfound = mnFound;
        mpReplaced = pMP;
    }

    for(map<KeyFrame*,size_t>::iterator mit=obs.begin(), mend=obs.end(); mit!=mend; mit++)
    {
        // Replace measurement in keyframe
        KeyFrame* pKF = mit->first;

        if(!pMP->IsInKeyFrame(pKF))
        {
            pKF->ReplaceMapPointMatch(mit->second, pMP);
            pMP->AddObservation(pKF,mit->second);
        }
        else
        {
            pKF->EraseMapPointMatch(mit->second);
        }
    }
    pMP->IncreaseFound(nfound);
    pMP->IncreaseVisible(nvisible);
    pMP->ComputeDistinctiveDescriptors();

    mpMap->EraseMapPoint(this);
}

bool MapPoint::isBad()
{
    unique_lock<mutex> lock(mMutexFeatures);
    unique_lock<mutex> lock2(mMutexPos);
    return mbBad;
}

void MapPoint::IncreaseVisible(int n)
{
    unique_lock<mutex> lock(mMutexFeatures);
    mnVisible+=n;
}

void MapPoint::IncreaseFound(int n)
{
    unique_lock<mutex> lock(mMutexFeatures);
    mnFound+=n;
}

float MapPoint::GetFoundRatio()
{
    unique_lock<mutex> lock(mMutexFeatures);
    return static_cast<float>(mnFound)/mnVisible;
}

void MapPoint::ComputeDistinctiveDescriptors()
{
    // Retrieve all observed descriptors
    vector<cv::Mat> vDescriptors;

    map<KeyFrame*,size_t> observations;

    {
        unique_lock<mutex> lock1(mMutexFeatures);
        if(mbBad)
            return;
        observations=mObservations;
    }

    if(observations.empty())
        return;

    vDescriptors.reserve(observations.size());

    for(map<KeyFrame*,size_t>::iterator mit=observations.begin(), mend=observations.end(); mit!=mend; mit++)
    {
        KeyFrame* pKF = mit->first;

        if(!pKF->isBad())
            vDescriptors.push_back(pKF->mDescriptors.row(mit->second));
    }

    if(vDescriptors.empty())
        return;

    // Compute distances between them
    const size_t N = vDescriptors.size();

    float Distances[N][N];
    for(size_t i=0;i<N;i++)
    {
        Distances[i][i]=0;
        for(size_t j=i+1;j<N;j++)
        {
            int distij = ORBmatcher::DescriptorDistance(vDescriptors[i],vDescriptors[j]);
            Distances[i][j]=distij;
            Distances[j][i]=distij;
        }
    }

    // Take the descriptor with least median distance to the rest
    int BestMedian = INT_MAX;
    int BestIdx = 0;
    for(size_t i=0;i<N;i++)
    {
        vector<int> vDists(Distances[i],Distances[i]+N);
        sort(vDists.begin(),vDists.end());
        int median = vDists[0.5*(N-1)];

        if(median<BestMedian)
        {
            BestMedian = median;
            BestIdx = i;
        }
    }

    {
        unique_lock<mutex> lock(mMutexFeatures);
        mDescriptor = vDescriptors[BestIdx].clone();
    }
}

cv::Mat MapPoint::GetDescriptor()
{
    unique_lock<mutex> lock(mMutexFeatures);
    return mDescriptor.clone();
}

int MapPoint::GetIndexInKeyFrame(KeyFrame *pKF)
{
    unique_lock<mutex> lock(mMutexFeatures);
    if(mObservations.count(pKF))
        return mObservations[pKF];
    else
        return -1;
}

bool MapPoint::IsInKeyFrame(KeyFrame *pKF)
{
    unique_lock<mutex> lock(mMutexFeatures);
    return (mObservations.count(pKF));
}

void MapPoint::UpdateNormalAndDepth()
{
    map<KeyFrame*,size_t> observations;
    KeyFrame* pRefKF;
    cv::Mat Pos;
    {
        unique_lock<mutex> lock1(mMutexFeatures);
        unique_lock<mutex> lock2(mMutexPos);
        if(mbBad)
            return;
        observations=mObservations;
        pRefKF=mpRefKF;
        Pos = mWorldPos.clone();
    }

    if(observations.empty())
        return;

    cv::Mat normal = cv::Mat::zeros(3,1,CV_32F);
    int n=0;
    for(map<KeyFrame*,size_t>::iterator mit=observations.begin(), mend=observations.end(); mit!=mend; mit++)
    {
        KeyFrame* pKF = mit->first;
        cv::Mat Owi = pKF->GetCameraCenter();
        cv::Mat normali = mWorldPos - Owi;
        normal = normal + normali/cv::norm(normali);
        n++;
    }

    cv::Mat PC = Pos - pRefKF->GetCameraCenter();
    const float dist = cv::norm(PC);
    const int level = pRefKF->mvKeysUn[observations[pRefKF]].octave;
    const float levelScaleFactor =  pRefKF->mvScaleFactors[level];
    const int nLevels = pRefKF->mnScaleLevels;

    {
        unique_lock<mutex> lock3(mMutexPos);
        mfMaxDistance = dist*levelScaleFactor;
        mfMinDistance = mfMaxDistance/pRefKF->mvScaleFactors[nLevels-1];
        mNormalVector = normal/n;
    }
}

bool MapPoint::RefineSubPix(KeyFrame* currentKF, size_t idx, int patchSize)
{
    static const bool verbose = 0;

    int halfPatchSize = (patchSize - 1)/ 2;

    // Retrieve the current pose of feature observations
    map<KeyFrame*,size_t> observations;
    {
        unique_lock<mutex> lock1(mMutexFeatures);
        if(mbBad)
            return false;
        observations=mObservations;
    }

    // We end if no observations are available
    if(observations.empty())
        return false;

    // Current patch
    mPatch = currentKF->mPatches[idx];

    // Finding the refPatch of the original detection
    int refIdx = 0;

    cv::Mat curObsNormal = mWorldPos - currentKF->GetCameraCenter();
    curObsNormal = curObsNormal / cv::norm(curObsNormal);

    double minAngle = 3.1415265;
    KeyFrame *patchRefKF;
    bool success = false;
    for(map<KeyFrame*,size_t>::iterator mit=observations.begin(), mend=observations.end(); mit!=mend; mit++) {
        KeyFrame *pKF = mit->first;


        // We compute the observation angle just for tests
        cv::Mat kfObsNormal = mWorldPos - pKF->GetCameraCenter();
        kfObsNormal = kfObsNormal / cv::norm(kfObsNormal);
        cv::Mat dotProduct = curObsNormal.t() * kfObsNormal;
        double angle = acos(dotProduct.at<float>(0,0)) ;

        // We do it agains the frame we first observed the feature
        if (pKF == mpRefKF) {
            patchRefKF = pKF;
            refIdx = mit->second;
            minAngle = angle;
            success = true;
            break;
        }


//        // We look for the keyframe with most similar observation angle
//        if ( angle < minAngle ) {
//            patchRefKF = pKF;
//            refIdx = mit->second;
//            minAngle = angle;
//        }
    }

    if (!success)
    {
        std::cout << "Ref keyframe was not found " << std::endl;
        exit(0);
    }


//    std::cout << "Chosen angle: " << minAngle * 180/3.1415265 << std::endl;

    cv::Mat refPatch = patchRefKF->mPatches[refIdx];
    cv::KeyPoint refKp = patchRefKF->mvKeys[refIdx];
    const int refOctave = patchRefKF->mvKeys[refIdx].octave;
    const float refKpScale =  patchRefKF->mvScaleFactors[refOctave];

    cv::Point2f currentKp = currentKF->mvKeys[idx].pt;

    if (!mPatch.empty() && !refPatch.empty()) {
        // Data we need to use
        //      mPatch      - patch in the current frame A
        //      refPatch    - patch in the previous frame B
        //                  - pose of B in A
        //                  - normal of the patch
        //                  - d of the plane
        //                  - keypoint positions and corresponding scale

        if (verbose) {
            std::cout << "-------------------------------" << std::endl;
            std::cout << "Subpix refinement " << std::endl;
        }

        float dataNormal[3] = { 0, 0, -1 };
        cv::Mat n = cv::Mat(3, 1, CV_32F, dataNormal);
//        std::cout <<"normal: " << n.t() << std::endl;

        // Frame a and b in world coordinates
        cv::Mat Twa = patchRefKF->GetPoseInverse();
        cv::Mat Twb = currentKF->GetPoseInverse();

        if(verbose)
            std::cout << "Keypoints: " << refKp.pt.x << " " << refKp.pt.y << " vs " << currentKF->mvKeys[idx].pt.x << " " << currentKF->mvKeys[idx].pt.y << std::endl;

        const int currentOctave = currentKF->mvKeys[idx].octave;
        const float currentKpScale =  currentKF->mvScaleFactors[currentOctave];

        if(verbose) {
            std::cout << "currentKpScale: " << currentKpScale << std::endl;
            std::cout << "refKpScale: " << refKpScale << std::endl;

            std::cout << "Keypoints in octave: " << refKp.pt.x / refKpScale << " " << refKp.pt.y / refKpScale << " vs "
                      << currentKF->mvKeys[idx].pt.x / currentKpScale << " "
                      << currentKF->mvKeys[idx].pt.y / currentKpScale << std::endl;
        }

        cv::Mat Ka = patchRefKF->getCameraMatrix();
//        std::cout << "Ka " << std::endl << Ka << std::endl;

        cv::Mat Kb = currentKF->getCameraMatrix();
//        std::cout << "Kb " << std::endl << Kb << std::endl;

        // Code for patch refinement
        PatchRefinement patchRefinement;


        //
//        std::cout << "mWorldPos; " << mWorldPos << std::endl;
        cv::Mat projection1 =  Ka * (patchRefKF->GetRotation() * mWorldPos + patchRefKF->GetTranslation());
        projection1 = patchRefinement.normalize2D(projection1);
        cv::Mat projection2 =  Kb * (currentKF->GetRotation() * mWorldPos + currentKF->GetTranslation());
        projection2 = patchRefinement.normalize2D(projection2);
//        std::cout << "Projected point 1: " << projection1 << std::endl;
//        std::cout << "Projected point 2: " << projection2 << std::endl;

        float img1ReprojError = std::sqrt(pow(projection1.at<float>(0,0) - refKp.pt.x, 2) + pow(projection1.at<float>(1,0) - refKp.pt.y, 2));
        float img2ReprojError = std::sqrt(pow(projection2.at<float>(0,0) - currentKp.x, 2) + pow(projection2.at<float>(1,0) - currentKp.y, 2));


        cv::Mat pointInA = patchRefKF->GetRotation() * mWorldPos + patchRefKF->GetTranslation();
        double d = patchRefinement.getDistanceToPlane(pointInA, n);
        cv::Mat H = patchRefinement.ComputeHomography(Twa, Twb, n, d, Ka, Kb);

        if(verbose)
            std::cout << "Homography" << std::endl << H << std::endl;

        // Lets compute the warped patch
        Eigen::Matrix3d Heig = patchRefinement.cv2eigen(H);

        // Homography test
        Eigen::Vector3d kp2 = Eigen::Vector3d(currentKp.x,currentKp.y, 1.0);
        Eigen::Vector3d kp2In1 = Heig.inverse()*kp2;
        kp2In1 = kp2In1 / kp2In1[2];

        if(verbose)
            std::cout << "Kp2 in img1: " <<  kp2In1[0]  <<" " << kp2In1[1] << std::endl;

        double dist = sqrt((refKp.pt.x - kp2In1[0]) * (refKp.pt.x - kp2In1[0]) + (refKp.pt.y - kp2In1[1])*(refKp.pt.y - kp2In1[1])) / refKpScale;

        if(verbose)
            std::cout << "DIST: " << dist << std::endl;


        if(verbose)
            std::cout << "minAngle: " << minAngle * 180/3.1415265 << "\t Dist: " << dist << " Reproj1: " << img1ReprojError<< " Reproj2: " << img2ReprojError << std::endl;

//        if ( dist > 15) {
//            std::cout << "Kp1: " << refKp.pt.x << " " << refKp.pt.y << " Projected: " << projection1.at<float>(0,0) << " " << projection1.at<float>(1,0) << std::endl;
//            std::cout << "Kp2 " << currentKp.x << " " << currentKp.y << " Projected: " << projection2.at<float>(0,0) << " " << projection2.at<float>(1,0) << std::endl;
//            std::cout << "World: " << mWorldPos << std::endl;
//            std::cout << "H: " << std::endl << H << std::endl;
//            std::cout << "Kp2in2 " << kp2In1[0] << " " << kp2In1[1] << std::endl;
//            std::cout << "Tab: " << std::endl << Twa.inv() * Twb << std::endl;
//            std::cout << "point in A: " << std::endl << patchRefKF->GetRotation() * mWorldPos + patchRefKF->GetTranslation() << std::endl;
//            std::cout << "Ka: " << std::endl << Ka << " " << Kb << std::endl;
//
//        }

        // TODO: Threshold to test
        if ( dist < 4) {

            if(verbose) {
                std::cout << "refPatch " << std::endl << refPatch.cols << " " << refPatch.rows << std::endl;
                std::cout << "curPatch " << std::endl << mPatch.cols << " " << mPatch.rows << std::endl;
            }

            // Current patch
            std::vector<Eigen::Vector2d> gradient;
            Eigen::Matrix2d HessianInv;
            patchRefinement.computeImageGradient(refPatch, gradient, HessianInv);


            std::vector<Eigen::Vector2d> subPosGradient(patchSize * patchSize, Eigen::Vector2d::Zero());
            double centerPos = (refPatch.rows - 1) / 2;

            std::vector<double> currentPatch = patchRefinement.computePatch(mPatch, centerPos, centerPos,
                                                                             mPatch.rows, patchSize, gradient,
                                                                             subPosGradient);

            if(verbose) {
                std::cout << "Current patch" << std::endl;
                patchRefinement.printPatch(currentPatch, patchSize);
            }


            std::vector<double> refPatchDouble = patchRefinement.computePatchOnSubImage(refPatch, patchSize, Heig.inverse(),
                                                                                        currentKp, currentKpScale,
                                                                                        refKp.pt, refKpScale);
            if (refPatchDouble.size()  == 0) {
                if (verbose)
                    std::cout << "\tRefPatch is empty -> warp was outside saved subImg" << std::endl;
                return false;
            }



            if(verbose) {
                std::cout << "RefPatch patch" << std::endl;
                patchRefinement.printPatch(refPatchDouble, patchSize);
            }


            // Lets start the optimization from provided position
            float optX = centerPos, optY = centerPos;

            // We perform 100 iterations or until stop condition is met
            for (int i = 0; i < 100; i++) {
                if(verbose)
                    std::cout <<"Iteration: " << i << " Position: " << optX << ", " << optY << std::endl;

                // Might happen if hessian is not invertible
                if(isnan(optX) || isnan(optY)) {
                    if (verbose)
                        std::cout << "\tNaN positions. Probably Hessian is not invertible" << std::endl;
                    return false;
                }

                // Compute new patch
                std::vector<double> currentPatch = patchRefinement.computePatch(mPatch, optX, optY,
                                                                                mPatch.rows, patchSize, gradient,
                                                                                subPosGradient);

                // Estimate residuals
                Eigen::Vector2d res = patchRefinement.computePatchDifference(refPatchDouble, currentPatch, subPosGradient);

                // Obtained error
                if(verbose)
                    std::cout << "Error: " << std::endl << res.transpose() << std::endl;

                // Step based on Hessian
                Eigen::Vector2d step = HessianInv * res;

                // Obtained error
                if(verbose)
                    std::cout << "Proposed step: " << std::endl << step.transpose() << std::endl;

                // Stop condition
                if (step.squaredNorm() < 0.0001)
                    break;

                // We do the step
                optX = optX + step[0];
                optY = optY + step[1];

                // dx, dy is too big - we might go to inf
//                double dx = (optX - centerPos)*currentKpScale;
//                double dy = (optY - centerPos)*currentKpScale;
//
//                if (sqrt(dx*dx+dy*dy) > halfPatchSize) {
                if(optX > (centerPos + halfPatchSize / 2) || optX < (centerPos - halfPatchSize / 2) || optY > (centerPos + halfPatchSize / 2) || optY < (centerPos - halfPatchSize / 2)) {
                    if (verbose)
                        std::cout << "\tWe moved too much! optX = " << optX << " optY = " << optY << std::endl;
                    return false;
                }

            }

            if(verbose) {
                std::cout << "Final position diff in keypoints octave: " << optX - centerPos << ", " << optY - centerPos
                          << std::endl;
                std::cout << "Base image changed: " << (optX - centerPos) * currentKpScale << ", "
                          << (optY - centerPos) * currentKpScale << std::endl;
            }

            currentKF->mvKeys[idx].pt.x = currentKF->mvKeys[idx].pt.x + (optX - centerPos)*currentKpScale;
            currentKF->mvKeys[idx].pt.y = currentKF->mvKeys[idx].pt.y + (optY - centerPos)*currentKpScale;

            // TODO: We should remove distortion !!!
            currentKF->mvKeysUn[idx].pt = currentKF->mvKeys[idx].pt;

            if (verbose)
                std::cout << "\tSuccess !" << std::endl;
            return true;
        }
        else if (verbose)
            std::cout <<"\tDistance above threshold" << std::endl;

    }
    else if (verbose)
        std::cout << "\tPatch missing" << std::endl;

    return false;
//    double minNormalError = 0;
//
//    cv::Mat currentNormal = mWorldPos - currentKF->GetCameraCenter();
//     Let's find the keyframe with the closest normal to current
//    for(map<KeyFrame*,size_t>::iterator mit=observations.begin(), mend=observations.end(); mit!=mend; mit++)
//    {
//        KeyFrame* pKF = mit->first;
//        int idx = mit->second;
//
//        if(!pKF->isBad()) {
//            cv::Mat Owi = pKF->GetCameraCenter();
//            cv::Mat normali = mWorldPos - Owi;
//
//            double error = cv::norm(currentNormal, normali);
//            if (minIter == mit || minNormalError > error ) {
//                minNormalError = error;
//            }
//        }
//    }
//
//    KeyFrame* relativeKF = minIter->first;
//    int relativeKF_idx = minIter->second;


//        cv::Mat Tpc = currentKF->GetPoseInverse() * mWorldPos;
//        float depth = Tpc.at<float>(2);
//        ComputeHomography(patchRefKF, currentKF, mNormalVector, depth);



}

float MapPoint::GetMinDistanceInvariance()
{
    unique_lock<mutex> lock(mMutexPos);
    return 0.8f*mfMinDistance;
}

float MapPoint::GetMaxDistanceInvariance()
{
    unique_lock<mutex> lock(mMutexPos);
    return 1.2f*mfMaxDistance;
}

int MapPoint::PredictScale(const float &currentDist, KeyFrame* pKF)
{
    float ratio;
    {
        unique_lock<mutex> lock(mMutexPos);
        ratio = mfMaxDistance/currentDist;
    }

    int nScale = ceil(log(ratio)/pKF->mfLogScaleFactor);
    if(nScale<0)
        nScale = 0;
    else if(nScale>=pKF->mnScaleLevels)
        nScale = pKF->mnScaleLevels-1;

    return nScale;
}

int MapPoint::PredictScale(const float &currentDist, Frame* pF)
{
    float ratio;
    {
        unique_lock<mutex> lock(mMutexPos);
        ratio = mfMaxDistance/currentDist;
    }

    int nScale = ceil(log(ratio)/pF->mfLogScaleFactor);
    if(nScale<0)
        nScale = 0;
    else if(nScale>=pF->mnScaleLevels)
        nScale = pF->mnScaleLevels-1;

    return nScale;
}


} //namespace ORB_SLAM
