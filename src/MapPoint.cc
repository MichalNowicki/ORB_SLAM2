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
    mpReplaced(static_cast<MapPoint*>(NULL)), mfMinDistance(0), mfMaxDistance(0), mpMap(pMap), refKFChanged(false)
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
    mnFound(1), mbBad(false), mpReplaced(NULL), mpMap(pMap), refKFChanged(false)
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

            if(mpRefKF==pKF) {// TODO: When we change ref keyframe, it is necessary to recompute patch matchings
//                std::cout << "RefKF was removed!" << std::endl;
                refKFChanged = true;
                mpRefKF = mObservations.begin()->first;
            }

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

bool MapPoint::RefineSubPix(KeyFrame* currentKF, size_t idx, int patchSize, int &numberOfIterations, double & errBefore, double &errAfter)
{
    const double reprojectionThr = 3;
    static const bool verbose = 0;

    enum PATCHVERSION { AGAINSTFIRST, AGAINSTCLOSESTANGLE };
    PATCHVERSION selectedPatchVersion = PATCHVERSION::AGAINSTFIRST;

     // We get the copy of the observations and world position of a point to work on
    map<KeyFrame*,size_t> observations;
    cv::Mat worldPos;
    {
        unique_lock<mutex> lock1(mMutexFeatures);
        unique_lock<mutex> lock2(mMutexPos);
        if(mbBad)
            return false;
        observations=mObservations;
        worldPos = mWorldPos;
    }

    // We end if no observations are available
    if(observations.empty())
        return false;

    // Current observation normal
    cv::Mat curObsNormal = worldPos - currentKF->GetCameraCenter();
    curObsNormal = curObsNormal / cv::norm(curObsNormal);

    // Finding the refPatch of the original detection or the one with closest observation angle
    double minAngle = 3.1415265;
    KeyFrame *patchRefKF;
    int refIdx = 0;
    bool success = false;

    // For all observations of the feature
    for(map<KeyFrame*,size_t>::iterator mit=observations.begin(), mend=observations.end(); mit!=mend; mit++) {
        KeyFrame *pKF = mit->first;

        // We compute the observation angle just for tests
        cv::Mat kfObsNormal = worldPos - pKF->GetCameraCenter();
        kfObsNormal = kfObsNormal / cv::norm(kfObsNormal);
        cv::Mat dotProduct = curObsNormal.t() * kfObsNormal;
        double angle = acos(dotProduct.at<float>(0,0)) ;

        // We look for the first frame and restore kp id
        if (selectedPatchVersion == PATCHVERSION::AGAINSTFIRST)
        {
            if (pKF == mpRefKF) {
                patchRefKF = pKF;
                refIdx = mit->second;
                minAngle = angle;
                success = true;
                break;
            }
        }
        else if(selectedPatchVersion == PATCHVERSION::AGAINSTCLOSESTANGLE)
        {
            // We look for the keyframe with most similar observation angle
            if ( angle < minAngle ) {
                patchRefKF = pKF;
                refIdx = mit->second;
                minAngle = angle;
                success = true;
            }
        }

    }

    // If no patchRefKF was found then something is completly wrong
    assert(success);

    // Ref patch in image A
    cv::Mat refPatch = patchRefKF->mPatches[refIdx];
    cv::Point2f refKp = patchRefKF->mvKeys[refIdx].pt;
    cv::Point2f refKpPatchLoc = patchRefKF->mPatchLocations[refIdx];
    const int refOctave = patchRefKF->mvKeys[refIdx].octave;
    const float refKpScale =  patchRefKF->mvScaleFactors[refOctave];

    // Current patch in image B
    mPatch = currentKF->mPatches[idx];
    cv::Point2f currentKp = currentKF->mvKeys[idx].pt;
    cv::Point2f curKpPatchLoc = currentKF->mPatchLocations[idx];
    const int currentOctave = currentKF->mvKeys[idx].octave;
    const float currentKpScale =  currentKF->mvScaleFactors[currentOctave];

    // TODO: Once it works, it is possible to just exchange the reference frame and recompute everything
    if ( refKFChanged ) {
        if (verbose)
            std::cout << "\tFeature rejected due to keyframe changed!" << std::endl;
        return false;
    }

    // Check if both patches exist - it may not exist if the feature is located close to image border
    if (mPatch.empty() || refPatch.empty()) {
        if (verbose)
            std::cout << "\tPatch missing" << std::endl;
        return false;
    }

    // Retrieve the camera matrices
    cv::Mat Ka = patchRefKF->getCameraMatrix();
    cv::Mat Kb = currentKF->getCameraMatrix();

    // Computation of the reprojection error based on 3D position
    cv::Mat pointInA = patchRefKF->GetRotation() * mWorldPos + patchRefKF->GetTranslation();
    cv::Mat projInA = Ka * pointInA;
    projInA = PatchRefinement::normalize2D(projInA);
    cv::Mat pointInB = currentKF->GetRotation() * mWorldPos + currentKF->GetTranslation();
    cv::Mat projInB = Kb * pointInB;
    projInB = PatchRefinement::normalize2D(projInB);
    float imgAReprojErrScaled = std::sqrt(
            pow(projInA.at<float>(0, 0) - refKp.x, 2) + pow(projInA.at<float>(1, 0) - refKp.y, 2)) /
                            refKpScale;
    float imgBReprojErrScaled = std::sqrt(
            pow(projInB.at<float>(0, 0) - currentKp.x, 2) + pow(projInB.at<float>(1, 0) - currentKp.y, 2)) /
                            currentKpScale;

    if (imgAReprojErrScaled > reprojectionThr || imgBReprojErrScaled > reprojectionThr) {
        if (verbose)
            std::cout << "\tThe reprojection error is too large to even try subpix refinement: " << imgAReprojErrScaled << " " << imgBReprojErrScaled << std::endl;
        return false;
    }


    // The patch normal in ref observation is assumed to be [0,0,-1]
    // TODO: Possible improvements are available when normal is estimated in a better way
    float dataNormal[3] = {0, 0, -1};
    cv::Mat n = cv::Mat(3, 1, CV_32F, dataNormal);

    // Computation of the distance to patch plane in image A
    double d = PatchRefinement::getDistanceToPlane(pointInA, n);

    // Computation of the homography between A and B
    cv::Mat Twa = patchRefKF->GetPoseInverse();
    cv::Mat Tbw = currentKF->GetPose();
    cv::Mat Tba = Tbw * Twa;
    cv::Mat H = PatchRefinement::computeHomography(Tba, n, d, Ka, Kb);


    // Homography error verification
    Eigen::Matrix3d Heig = PatchRefinement::cv2eigen(H);
    Eigen::Vector3d kp1 = Eigen::Vector3d(refKp.x, refKp.y, 1.0), kp2 = Eigen::Vector3d(currentKp.x, currentKp.y, 1.0);
    Eigen::Vector3d kp2In1 = Heig.inverse() * kp2, kp1In2 = Heig * kp1;
    kp1In2 = kp1In2 / kp1In2[2];
    kp2In1 = kp2In1 / kp2In1[2];
    double distA = sqrt((refKp.x - kp2In1[0]) * (refKp.x - kp2In1[0]) + (refKp.y - kp2In1[1]) * (refKp.y - kp2In1[1])) /
                   refKpScale;
    double distB = sqrt((currentKp.x - kp1In2[0]) * (currentKp.x - kp1In2[0]) +
                        (currentKp.y - kp1In2[1]) * (currentKp.y - kp1In2[1])) / currentKpScale;

    if (distA > reprojectionThr || distB > reprojectionThr) {
        if (verbose)
            std::cout << "\tThe homography error is too large to even try subpix refinement" << std::endl;
        return false;
    }


    // Let's optimize the final position
    cv::Point2f correction = cv::Point2f(0, 0);

    PatchRefinement patchRefinement(patchSize);
    success = patchRefinement.optimizePosition(refPatch, refKpPatchLoc, refKp, refKpScale,
                                               mPatch, curKpPatchLoc, currentKpScale,
                                               Heig, correction);

    // Success -> update the positon of the feature
    if (success) {
        if (verbose)
            std::cout << "\tSuccess !" << std::endl;

        // We add those increments to the mvKeys positions
        currentKF->mvKeys[idx].pt = currentKp + correction;

        // Removing the distortion
        currentKF->mvKeysUn[idx].pt = UndistortPoint(currentKF->mvKeys[idx].pt, currentKF->mK,
                                                     currentKF->mDistCoef);

        // Saving some statistics
        numberOfIterations = patchRefinement.numberOfPerformedIterations;
        errBefore = patchRefinement.errBefore;
        errAfter = patchRefinement.errAfter;
    }

    return true;

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


//cv::Point2f MapPoint::distortPoint(const cv::Point2f p, const float fx, const float fy, const float cx,
//                                   const float cy, const float invfx, const float invfy, const cv::Mat mDistCoef) {
//    float k1 = mDistCoef.at<float>(0);
//    float k2 = mDistCoef.at<float>(1);
//    float p1 = mDistCoef.at<float>(2);
//    float p2 = mDistCoef.at<float>(3);
//    float k3 = mDistCoef.at<float>(4);
//
//    // Normalized image coordinates
//    float x = (p.x - cx) / invfx;
//    float y = (p.y - cy) / invfy;
//
//    // r2 as it is easier to write those eq
//    float r2 = x*x + y*y;
//
//    // Radial
//    float xDis = x * (1 + k1 * r2 + k2 * r2 * r2 + k3 * r2 * r2 * r2);
//    float yDis = y * (1 + k1 * r2 + k2 * r2 * r2 + k3 * r2 * r2 * r2);
//
//    // Tangential
//    xDis = xDis + (2 * p1 * x * y + p2 * (r2 + 2 * x * x));
//    yDis = yDis + (p1 * (r2 + 2 * y * y) + 2 * p2 * x * y);
//
//    // Unnormalize coordinates
//    cv::Point2f ret;
//    ret.x = xDis * fx + cx;
//    ret.y = xDis * fy + cy;
//
//    return ret;
//}

cv::Point2f MapPoint::UndistortPoint(cv::Point2f p, const cv::Mat K, const cv::Mat mDistCoef) {
    if (mDistCoef.at<float>(0) == 0.0)
        return p;

    // Fill matrix with point
    cv::Mat mat(1, 2, CV_32F);

    mat.at<float>(0) = p.x;
    mat.at<float>(1) = p.y;

    // Undistort points
    mat=mat.reshape(2);
    cv::undistortPoints(mat, mat, K, mDistCoef, cv::Mat(), K);
    mat=mat.reshape(1);


    p.x = mat.at<float>(0);
    p.y = mat.at<float>(1);
    return p;
}

} //namespace ORB_SLAM
