//
// Created by michalnowicki on 08.09.17.
//

#include "OptimizerBA.h"

#include "Thirdparty/g2o/g2o/core/block_solver.h"
#include "Thirdparty/g2o/g2o/core/optimization_algorithm_levenberg.h"
#include "Thirdparty/g2o/g2o/core/optimization_algorithm_gauss_newton.h"
#include "Thirdparty/g2o/g2o/solvers/linear_solver_eigen.h"
#include "Thirdparty/g2o/g2o/types/types_six_dof_expmap.h"
#include "Thirdparty/g2o/g2o/core/robust_kernel_impl.h"
#include "Thirdparty/g2o/g2o/solvers/linear_solver_dense.h"
#include "Thirdparty/g2o/g2o/types/types_seven_dof_expmap.h"


#include<Eigen/StdVector>

#include "Converter.h"

#include<mutex>

namespace ORB_SLAM2 {
    void OptimizerBA::saveBAProblem(KeyFrame *pKF, bool* pbStopFlag, Map* pMap, float sigma, std::string name)
    {
        list<KeyFrame*> lLocalKeyFrames;

        lLocalKeyFrames.push_back(pKF);
        pKF->mnBALocalForKF = pKF->mnId;

        const vector<KeyFrame*> vNeighKFs = pKF->GetVectorCovisibleKeyFrames();
        std::cout << "pKF->GetVectorCovisibleKeyFrames().size() = " << vNeighKFs.size() << std::endl;
        for(int i=0, iend=vNeighKFs.size(); i<iend; i++)
        {
            KeyFrame* pKFi = vNeighKFs[i];
            pKFi->mnBALocalForKF = pKF->mnId;
            if(!pKFi->isBad())
                lLocalKeyFrames.push_back(pKFi);
        }
        std::cout << "lLocalKeyFrames.size() = " << lLocalKeyFrames.size() << std::endl;

        // Local MapPoints seen in Local KeyFrames
        list<MapPoint*> lLocalMapPoints;
        for(list<KeyFrame*>::iterator lit=lLocalKeyFrames.begin() , lend=lLocalKeyFrames.end(); lit!=lend; lit++)
        {
            vector<MapPoint*> vpMPs = (*lit)->GetMapPointMatches();
            for(vector<MapPoint*>::iterator vit=vpMPs.begin(), vend=vpMPs.end(); vit!=vend; vit++)
            {
                MapPoint* pMP = *vit;
                if(pMP)
                    if(!pMP->isBad())
                        if(pMP->mnBALocalForKF!=pKF->mnId)
                        {
                            lLocalMapPoints.push_back(pMP);
                            pMP->mnBALocalForKF=pKF->mnId;
                        }
            }
        }

        // TODO: Stupid code to repeat BA
        for(list<KeyFrame*>::iterator lit=lLocalKeyFrames.begin() , lend=lLocalKeyFrames.end(); lit!=lend; lit++)
        {
            vector<MapPoint*> vpMPs = (*lit)->GetMapPointMatches();
            for(vector<MapPoint*>::iterator vit=vpMPs.begin(), vend=vpMPs.end(); vit!=vend; vit++) {
                MapPoint *pMP = *vit;
                if (pMP)
                    if (!pMP->isBad())
                        pMP->mnBALocalForKF = pKF->mnId - 1;
            }
        }

        int obs = 0;
        for(list<MapPoint*>::iterator lit=lLocalMapPoints.begin(), lend=lLocalMapPoints.end(); lit!=lend; lit++) {
            MapPoint *pMP = *lit;

            const map<KeyFrame *, size_t> observations = pMP->GetObservations();
            //Set edges
            for (map<KeyFrame *, size_t>::const_iterator mit = observations.begin(), mend = observations.end();
                 mit != mend; mit++) {
                KeyFrame *pKFi = mit->first;

                if (!pKFi->isBad()) {
                    obs++;
                }
            }
        }


        std::ofstream saveStream(name.c_str());

        saveStream << lLocalKeyFrames.size() << " "<< lLocalMapPoints.size() << " " << obs << std::endl;

        // Set Local KeyFrame vertices
        unsigned long maxKFid = 0;
        for(list<KeyFrame*>::iterator lit=lLocalKeyFrames.begin(), lend=lLocalKeyFrames.end(); lit!=lend; lit++) {
            KeyFrame *pKFi = *lit;

            cv::Mat R = pKFi->GetRotation(), Rexp, t = pKFi->GetTranslation();
            cv::Rodrigues(R, Rexp);

            saveStream << pKFi->mnId << " " << Rexp.at<float>(0) << " " << Rexp.at<float>(1) << " " << Rexp.at<float>(2) << " "
                       << t.at<float>(0) << " " << t.at<float>(1) << " " << t.at<float>(2) << " " << pKFi->fx << " "
                       << pKFi->fy << " "
                       << pKFi->cx << " " << pKFi->cy << std::endl;

            // Set Local KeyFrame vertices
            for (list<KeyFrame *>::iterator lit = lLocalKeyFrames.begin(), lend = lLocalKeyFrames.end();
                 lit != lend; lit++) {
                KeyFrame *pKFi = *lit;
                if (pKFi->mnId > maxKFid)
                    maxKFid = pKFi->mnId;
            }
        }

        for(list<MapPoint*>::iterator lit=lLocalMapPoints.begin(), lend=lLocalMapPoints.end(); lit!=lend; lit++) {
            MapPoint *pMP = *lit;

            int id = pMP->mnId + maxKFid + 1;
            Eigen::Vector3d worldPos = Converter::toVector3d(pMP->GetWorldPos());
            saveStream << id << " " << worldPos[0] << " " << worldPos[1] << " " << worldPos[2] << std::endl;
        }


        for(list<MapPoint*>::iterator lit=lLocalMapPoints.begin(), lend=lLocalMapPoints.end(); lit!=lend; lit++) {
            MapPoint *pMP = *lit;
            int id = pMP->mnId + maxKFid + 1;

            const map<KeyFrame *, size_t> observations = pMP->GetObservations();

            //Set edges
            for (map<KeyFrame *, size_t>::const_iterator mit = observations.begin(), mend = observations.end();
                 mit != mend; mit++) {
                KeyFrame *pKFi = mit->first;

                if (!pKFi->isBad()) {
                    const cv::KeyPoint &kpUn = pKFi->mvKeysUn[mit->second];
                    const float &invSigma2 = pKFi->mvInvLevelSigma2[kpUn.octave];

                    saveStream << pKFi->mnId << " " << id << " " << kpUn.pt.x << " " << kpUn.pt.y << " " << invSigma2 << std::endl;
                }
            }
        }

        saveStream.close();
    }

    void OptimizerBA::VerifyReprojectionError(list<KeyFrame *> lLocalKeyFrames, set<MapPoint *> lLocalMapPoints, list<KeyFrame *> lFixedCameras) {

        std::vector<double> reprojectionErr;
        for (set<MapPoint *>::iterator lit = lLocalMapPoints.begin(), lend = lLocalMapPoints.end(); lit != lend; lit++) {
            MapPoint *pMP = *lit;
            cv::Mat mWorldPos = pMP->GetWorldPos();

            const map<KeyFrame *, size_t> observations = pMP->GetObservations();

//Set edges
            for (map<KeyFrame *, size_t>::const_iterator mit = observations.begin(), mend = observations.end();
                 mit != mend; mit++) {
                KeyFrame *pKFi = mit->first;

                if (!pKFi->isBad()) {
                    const cv::KeyPoint &kpUn = pKFi->mvKeysUn[mit->second];

                    cv::Mat pointInA = pKFi->GetRotation() * mWorldPos + pKFi->GetTranslation();
                    cv::Mat projectionInA = pKFi->mK * pointInA;
                    double z = projectionInA.at<float>(2, 0);
                    projectionInA.col(0) = projectionInA.col(0) / z;


                    const float refKpScale = pKFi->mvScaleFactors[kpUn.octave];

                    float img1ReprojError = std::sqrt(pow(projectionInA.at<float>(0, 0) - kpUn.pt.x, 2) +
                                                      pow(projectionInA.at<float>(1, 0) - kpUn.pt.y, 2)) / refKpScale;

//                std::cout << projectionInA.at<float>(0, 0) << " " << projectionInA.at<float>(1, 0) << " " << kpUn.pt.x << " " << kpUn.pt.y << std::endl;

                    if (std::isnan(img1ReprojError))
                        std::cout << "NAN" << std::endl;
                    reprojectionErr.push_back(img1ReprojError);
                }
            }


        }

        std::cout << "\tAvg reprojection error computed in each scale: " << accumulate( reprojectionErr.begin(), reprojectionErr.end(), 0.0)/reprojectionErr.size() << std::endl;
    }


    void OptimizerBA::prepareKFsAndMapPoints(KeyFrame *pKF, Map *pMap, list<KeyFrame *> &lLocalKeyFrames, set<MapPoint *> &lLocalMapPoints, list<KeyFrame *> &lFixedCameras) {
        // Local KeyFrames: First Breath Search from Current Keyframe
        lLocalKeyFrames.clear();
        lLocalKeyFrames.push_back(pKF);
        pKF->mnBALocalForKF = pKF->mnId;

        const vector<KeyFrame *> vNeighKFs = pKF->GetVectorCovisibleKeyFrames();
        for (int i = 0, iend = vNeighKFs.size(); i < iend; i++) {
            KeyFrame *pKFi = vNeighKFs[i];
            pKFi->mnBALocalForKF = pKF->mnId;
            if (!pKFi->isBad())
                lLocalKeyFrames.push_back(pKFi);
        }

        // Local MapPoints seen in Local KeyFrames
        lLocalMapPoints.clear();
        for (list<KeyFrame *>::iterator lit = lLocalKeyFrames.begin(), lend = lLocalKeyFrames.end(); lit != lend; lit++) {
            vector<MapPoint *> vpMPs = (*lit)->GetMapPointMatches();
            for (vector<MapPoint *>::iterator vit = vpMPs.begin(), vend = vpMPs.end(); vit != vend; vit++) {
                MapPoint *pMP = *vit;
                if (pMP)
                    if (!pMP->isBad()) {
                            lLocalMapPoints.insert(pMP);
                            pMP->mnBALocalForKF = pKF->mnId;
                        }
            }
        }

        // Fixed Keyframes. Keyframes that see Local MapPoints but that are not Local Keyframes
        lFixedCameras.clear();
        for (set<MapPoint *>::iterator lit = lLocalMapPoints.begin(), lend = lLocalMapPoints.end(); lit != lend; lit++) {
            map<KeyFrame *, size_t> observations = (*lit)->GetObservations();
            for (map<KeyFrame *, size_t>::iterator mit = observations.begin(), mend = observations.end();
                 mit != mend; mit++) {
                KeyFrame *pKFi = mit->first;

                if (pKFi->mnBALocalForKF != pKF->mnId && pKFi->mnBAFixedForKF != pKF->mnId) {
                    pKFi->mnBAFixedForKF = pKF->mnId;
                    if (!pKFi->isBad())
                        lFixedCameras.push_back(pKFi);
                }
            }
        }
    }


    std::vector<double>  OptimizerBA::LocalBundleAdjustment(KeyFrame *pKF, bool *pbStopFlag, Map *pMap, OptimizerBA::TYPE type, float sigma)
    {
        list<KeyFrame *> lLocalKeyFrames, lFixedCameras;
        set<MapPoint *> lLocalMapPoints;

        // Find the keyframes that will take part in optimization
        OptimizerBA::prepareKFsAndMapPoints(pKF, pMap, lLocalKeyFrames, lLocalMapPoints, lFixedCameras);

        if (type == OptimizerBA::TYPE::INVERSE_DEPTH)
            return OptimizerBA::BundleAdjustmentInvDepth(lLocalKeyFrames, lFixedCameras, lLocalMapPoints, pbStopFlag, pMap, sigma);
        if (type == OptimizerBA::TYPE::INVERSE_DEPTH_SINGLE_PARAM)
            return OptimizerBA::BundleAdjustmentInvDepthSingleParam(lLocalKeyFrames, lFixedCameras, lLocalMapPoints, pbStopFlag, pMap, sigma);
        if (type == OptimizerBA::TYPE::INVERSE_DEPTH_SINGLE_PARAM_PATCH)
            return OptimizerBA::BundleAdjustmentInvDepthSingleParamPatch(lLocalKeyFrames, lFixedCameras, lLocalMapPoints, pbStopFlag, pMap, sigma);
        if (type == OptimizerBA::TYPE::INVERSE_DEPTH_SINGLE_PARAM_PATCH_BRIGHTNESS)
            return OptimizerBA::BundleAdjustmentInvDepthSingleParamPatchBright(lLocalKeyFrames, lFixedCameras, lLocalMapPoints, pbStopFlag, pMap, sigma);
    }


    std::vector<double>  OptimizerBA::GlobalBundleAdjustment(Map* pMap, OptimizerBA::TYPE type, float sigma)
    {
        vector<KeyFrame*> vpKFs = pMap->GetAllKeyFrames();
        vector<MapPoint*> vpMP = pMap->GetAllMapPoints();

        bool pbStopFlag = false; // TODO
        list<KeyFrame *> lLocalKeyFrames(vpKFs.begin(), vpKFs.end()), lFixedCameras;
        set<MapPoint *> lLocalMapPoints(vpMP.begin(), vpMP.end());

        if (type == OptimizerBA::TYPE::INVERSE_DEPTH)
            return OptimizerBA::BundleAdjustmentInvDepth(lLocalKeyFrames, lFixedCameras, lLocalMapPoints, &pbStopFlag, pMap, sigma);
        if (type == OptimizerBA::TYPE::INVERSE_DEPTH_SINGLE_PARAM)
            return OptimizerBA::BundleAdjustmentInvDepthSingleParam(lLocalKeyFrames, lFixedCameras, lLocalMapPoints, &pbStopFlag, pMap, sigma);
        if (type == OptimizerBA::TYPE::INVERSE_DEPTH_SINGLE_PARAM_PATCH)
            return OptimizerBA::BundleAdjustmentInvDepthSingleParamPatch(lLocalKeyFrames, lFixedCameras, lLocalMapPoints, &pbStopFlag, pMap, sigma);
        if (type == OptimizerBA::TYPE::INVERSE_DEPTH_SINGLE_PARAM_PATCH_BRIGHTNESS)
            return OptimizerBA::BundleAdjustmentInvDepthSingleParamPatchBright(lLocalKeyFrames, lFixedCameras, lLocalMapPoints, &pbStopFlag, pMap, sigma);
    }


    std::vector<double> OptimizerBA::BundleAdjustmentInvDepth(list<KeyFrame *> lLocalKeyFrames, list<KeyFrame *> lFixedCameras, set<MapPoint *> lLocalMapPoints, bool *pbStopFlag, Map *pMap, float sigma)
    {
        const float thHuberMono = 5.991;
        const int blockSolverCameras = 6;
        const int blockSolverPoses = 1;

        // Setup optimizer
        g2o::SparseOptimizer optimizer;

        typedef g2o::BlockSolver< g2o::BlockSolverTraits<blockSolverCameras, blockSolverPoses> > OurBlockSolver;
        OurBlockSolver::LinearSolverType *linearSolver;
        linearSolver = new g2o::LinearSolverEigen<OurBlockSolver::PoseMatrixType>();
        OurBlockSolver *solver_ptr = new OurBlockSolver(linearSolver);

        g2o::OptimizationAlgorithmLevenberg *solver = new g2o::OptimizationAlgorithmLevenberg(solver_ptr);
        optimizer.setAlgorithm(solver);

        if (pbStopFlag)
            optimizer.setForceStopFlag(pbStopFlag);


        KeyFrame *pKF = *lLocalKeyFrames.begin(); // TODO: Now setting it as camera params of 1st camera
        Eigen::Vector2d principal_point(pKF->cx, pKF->cy);
        g2o::CameraParameters * cam_params
                = new g2o::CameraParameters (pKF->fx, pKF->fy, principal_point, 0.);
        cam_params->setId(0);
        if (!optimizer.addParameter(cam_params)){
            assert(false);
        }


        unsigned long maxKFid = 0;

        // Set Local KeyFrame vertices
        for (list<KeyFrame *>::iterator lit = lLocalKeyFrames.begin(), lend = lLocalKeyFrames.end(); lit != lend; lit++) {
            KeyFrame *pKFi = *lit;

            g2o::VertexSE3Expmap *v = new g2o::VertexSE3Expmap();
            v->setEstimate(Converter::toSE3Quat(pKFi->GetPose()));
            v->setId(pKFi->mnId);
            v->setFixed(pKFi->mnId == 0);
            optimizer.addVertex(v);
            if (pKFi->mnId > maxKFid)
                maxKFid = pKFi->mnId;
        }

        // Set Fixed KeyFrame vertices
        for (list<KeyFrame *>::iterator lit = lFixedCameras.begin(), lend = lFixedCameras.end(); lit != lend; lit++) {
            KeyFrame *pKFi = *lit;
            g2o::VertexSE3Expmap *v = new g2o::VertexSE3Expmap();
            v->setEstimate(Converter::toSE3Quat(pKFi->GetPose()));
            v->setId(pKFi->mnId);
            v->setFixed(true);
            optimizer.addVertex(v);
            if (pKFi->mnId > maxKFid)
                maxKFid = pKFi->mnId;
        }

        // Set MapPoint vertices
        const int nExpectedSize = (lLocalKeyFrames.size() + lFixedCameras.size()) * lLocalMapPoints.size();

        vector<g2o::EdgeProjectPSI2UV *> vpEdgesMono;
        vpEdgesMono.reserve(nExpectedSize);

        vector<KeyFrame *> vpEdgeKFMono;
        vpEdgeKFMono.reserve(nExpectedSize);

        vector<MapPoint *> vpMapPointEdgeMono;
        vpMapPointEdgeMono.reserve(nExpectedSize);


        // For all points
        for (set<MapPoint *>::iterator lit = lLocalMapPoints.begin(), lend = lLocalMapPoints.end(); lit != lend; lit++) {
            MapPoint *pMP = *lit;
            KeyFrame *refKF = pMP->GetReferenceKeyFrame();
            const map<KeyFrame *, size_t> observations = pMP->GetObservations();


            g2o::VertexSBAPointXYZ *vPoint = new g2o::VertexSBAPointXYZ();

            // Moving from world coordinate to local coordinate in ref kf and then saving depth as inverse Z
            cv::Mat worldPos = pMP->GetWorldPos();
            cv::Mat pointInRef = refKF->GetRotation() * worldPos + refKF->GetTranslation();

            Eigen::Vector3d pointInverted;
            pointInverted[0] = pointInRef.at<float>(0) / pointInRef.at<float>(2);
            pointInverted[1] = pointInRef.at<float>(1) / pointInRef.at<float>(2);
            pointInverted[2] = 1. / pointInRef.at<float>(2);
            vPoint->setEstimate(pointInverted );


            int id = pMP->mnId + maxKFid + 1;
            vPoint->setId(id);
            vPoint->setMarginalized(true);
            optimizer.addVertex(vPoint);

//Set edges
            for (map<KeyFrame *, size_t>::const_iterator mit = observations.begin(), mend = observations.end();
                 mit != mend; mit++) {
                KeyFrame *pKFi = mit->first;

                if (!pKFi->isBad()) {
                    const cv::KeyPoint &kpUn = pKFi->mvKeysUn[mit->second];

// Monocular observation
                    if (pKFi->mvuRight[mit->second] < 0) {

                        Eigen::Matrix<double, 2, 1> obs;
                        obs << kpUn.pt.x, kpUn.pt.y;

                        g2o::EdgeProjectPSI2UV *e = new g2o::EdgeProjectPSI2UV();
                        e->resize(3);

                        // RefKp might not have been added yet
                        g2o::OptimizableGraph::Vertex * v = dynamic_cast<g2o::OptimizableGraph::Vertex *>(optimizer.vertex(refKF->mnId));
                        if (!v)
                        {
                            g2o::VertexSE3Expmap *vSE3 = new g2o::VertexSE3Expmap();
                            vSE3->setEstimate(Converter::toSE3Quat(refKF->GetPose()));
                            vSE3->setId(refKF->mnId);
                            vSE3->setFixed(true);
                            optimizer.addVertex(vSE3);
                        }

                        e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex *>(optimizer.vertex(id)));
                        e->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex *>(optimizer.vertex(pKFi->mnId)));
                        e->setVertex(2, dynamic_cast<g2o::OptimizableGraph::Vertex *>(optimizer.vertex(refKF->mnId)));

                        e->setMeasurement(obs);
                        const float &invSigma2 = pKFi->mvInvLevelSigma2[kpUn.octave];
                        const float assumedDetInvSigma2 = 1.0 / (sigma * sigma);
                        e->setInformation(Eigen::Matrix2d::Identity() * invSigma2 * assumedDetInvSigma2);

                        g2o::RobustKernelHuber *rk = new g2o::RobustKernelHuber;
                        e->setRobustKernel(rk);
                        rk->setDelta(sqrt(thHuberMono));

                        e->setParameterId(0, 0);
                        optimizer.addEdge(e);
                        vpEdgesMono.push_back(e);
                        vpEdgeKFMono.push_back(pKFi);
                        vpMapPointEdgeMono.push_back(pMP);

                    }
                }
            }
        }

        if (pbStopFlag)
            if (*pbStopFlag)
                return std::vector<double>(); // Modified


        // Optimization for 5 iterations
        optimizer.initializeOptimization();
        optimizer.optimize(5);

        bool bDoMore = true;

        if (pbStopFlag)
            if (*pbStopFlag)
                bDoMore = false;

        if (bDoMore) {

            // Check inlier observations
            for (size_t i = 0, iend = vpEdgesMono.size(); i < iend; i++) {
                g2o::EdgeProjectPSI2UV *e = vpEdgesMono[i];
                MapPoint *pMP = vpMapPointEdgeMono[i];

                if (pMP->isBad())
                    continue;

                if (e->chi2() > thHuberMono || !e->isDepthPositive()) {
                    e->setLevel(1);
                }

                e->setRobustKernel(0);
            }


            // Optimize again without the outliers
            optimizer.initializeOptimization(0);
            optimizer.optimize(10);
        }

        vector<pair<KeyFrame *, MapPoint *> > vToErase;
        vToErase.reserve(vpEdgesMono.size());

        // Check inlier observations
        std::vector<double> chi2Statistics;
        for (size_t i = 0, iend = vpEdgesMono.size(); i < iend; i++) {
            g2o::EdgeProjectPSI2UV *e = vpEdgesMono[i];
            MapPoint *pMP = vpMapPointEdgeMono[i];

            if (pMP->isBad())
                continue;

            if (e->chi2() > thHuberMono || !e->isDepthPositive()) {
                KeyFrame *pKFi = vpEdgeKFMono[i];
                vToErase.push_back(make_pair(pKFi, pMP));
            } else
                chi2Statistics.push_back(e->chi2());

        }
        std::cout << "\tLocalBA: avg chi2() = "
                  << accumulate(chi2Statistics.begin(), chi2Statistics.end(), 0.0) / chi2Statistics.size() << " over "
                  << chi2Statistics.size() << " measurements" << std::endl;


        // Get Map Mutex
        unique_lock<mutex> lock(pMap->mMutexMapUpdate);

        if (!vToErase.empty()) {
            for (size_t i = 0; i < vToErase.size(); i++) {
                KeyFrame *pKFi = vToErase[i].first;
                MapPoint *pMPi = vToErase[i].second;
                pKFi->EraseMapPointMatch(pMPi);
                pMPi->EraseObservation(pKFi);
            }
        }

// Recover optimized data

//Keyframes
        for (list<KeyFrame *>::iterator lit = lLocalKeyFrames.begin(), lend = lLocalKeyFrames.end(); lit != lend; lit++) {
            KeyFrame *pKF = *lit;
            g2o::VertexSE3Expmap *vSE3 = static_cast<g2o::VertexSE3Expmap *>(optimizer.vertex(pKF->mnId));
            g2o::SE3Quat SE3quat = vSE3->estimate();
            pKF->SetPose(Converter::toCvMat(SE3quat));
        }

//Points
        for (set<MapPoint *>::iterator lit = lLocalMapPoints.begin(), lend = lLocalMapPoints.end(); lit != lend; lit++) {
            MapPoint *pMP = *lit;

            g2o::VertexSBAPointXYZ *vPoint = static_cast<g2o::VertexSBAPointXYZ *>(optimizer.vertex(
                    pMP->mnId + maxKFid + 1));
            KeyFrame *refKF = pMP->GetReferenceKeyFrame();

            Eigen::Vector3d pointEstimate = vPoint->estimate();

            cv::Mat pointInFirst = cv::Mat(3,1, CV_32F);
            pointInFirst.at<float>(0) = pointEstimate[0] / pointEstimate[2];
            pointInFirst.at<float>(1) = pointEstimate[1] / pointEstimate[2];
            pointInFirst.at<float>(2) = 1. / pointEstimate[2];


            cv::Mat worldPointPos = refKF->GetRotation().t() * pointInFirst - refKF->GetRotation().t() * refKF->GetTranslation();

            pMP->SetWorldPos(worldPointPos);
            pMP->UpdateNormalAndDepth();
        }



        // We verify the reprojection error
        OptimizerBA::VerifyReprojectionError(lLocalKeyFrames, lLocalMapPoints, lFixedCameras);

        return chi2Statistics;
    }




    std::vector<double> OptimizerBA::BundleAdjustmentInvDepthSingleParam(list<KeyFrame *> lLocalKeyFrames, list<KeyFrame *> lFixedCameras, set<MapPoint *> lLocalMapPoints, bool *pbStopFlag, Map *pMap, float sigma)
    {
        const float thHuberMono = sqrt(5.991);
        const int blockSolverCameras = 6;
        const int blockSolverPoses = 1;

        // Setup optimizer
        g2o::SparseOptimizer optimizer;

        typedef g2o::BlockSolver< g2o::BlockSolverTraits<blockSolverCameras, blockSolverPoses> > OurBlockSolver;
        OurBlockSolver::LinearSolverType *linearSolver;
        linearSolver = new g2o::LinearSolverEigen<OurBlockSolver::PoseMatrixType>();
        OurBlockSolver *solver_ptr = new OurBlockSolver(linearSolver);

        g2o::OptimizationAlgorithmLevenberg *solver = new g2o::OptimizationAlgorithmLevenberg(solver_ptr);
        optimizer.setAlgorithm(solver);

        if (pbStopFlag)
            optimizer.setForceStopFlag(pbStopFlag);

        KeyFrame *pKF = *lLocalKeyFrames.begin(); // TODO: Now setting it as camera params of 1st camera
        Eigen::Vector2d principal_point(pKF->cx, pKF->cy);
        g2o::CameraParameters * cam_params
                = new g2o::CameraParameters (pKF->fx, pKF->fy, principal_point, 0.);
        cam_params->setId(0);
        if (!optimizer.addParameter(cam_params)){
            assert(false);
        }


        unsigned long maxKFid = 0;

        // Set Local KeyFrame vertices
        for (list<KeyFrame *>::iterator lit = lLocalKeyFrames.begin(), lend = lLocalKeyFrames.end(); lit != lend; lit++) {
            KeyFrame *pKFi = *lit;

            g2o::VertexSE3Expmap *v = new g2o::VertexSE3Expmap();
            v->setEstimate(Converter::toSE3Quat(pKFi->GetPose()));
            v->setId(pKFi->mnId);
            v->setFixed(pKFi->mnId == 0);
            optimizer.addVertex(v);
            if (pKFi->mnId > maxKFid)
                maxKFid = pKFi->mnId;
        }

        // Set Fixed KeyFrame vertices
        for (list<KeyFrame *>::iterator lit = lFixedCameras.begin(), lend = lFixedCameras.end(); lit != lend; lit++) {
            KeyFrame *pKFi = *lit;
            g2o::VertexSE3Expmap *v = new g2o::VertexSE3Expmap();
            v->setEstimate(Converter::toSE3Quat(pKFi->GetPose()));
            v->setId(pKFi->mnId);
            v->setFixed(true);
            optimizer.addVertex(v);
            if (pKFi->mnId > maxKFid)
                maxKFid = pKFi->mnId;
        }

        // Set MapPoint vertices
        const int nExpectedSize = (lLocalKeyFrames.size() + lFixedCameras.size()) * lLocalMapPoints.size();


        vector<g2o::EdgeProjectPSI2UVSingleParam *> vpEdgesMono;
        vpEdgesMono.reserve(nExpectedSize);

        vector<KeyFrame *> vpEdgeKFMono;
        vpEdgeKFMono.reserve(nExpectedSize);

        vector<MapPoint *> vpMapPointEdgeMono;
        vpMapPointEdgeMono.reserve(nExpectedSize);


        for (set<MapPoint *>::iterator lit = lLocalMapPoints.begin(), lend = lLocalMapPoints.end(); lit != lend; lit++) {
            MapPoint *pMP = *lit;
            KeyFrame *refKF = pMP->GetReferenceKeyFrame();
            const map<KeyFrame *, size_t> observations = pMP->GetObservations();


            g2o::VertexSBAPointInvD *vPoint = new g2o::VertexSBAPointInvD();

            // Moving from world coordinate to local coordinate in ref kf and then saving depth as inverse Z
            cv::Mat worldPos = pMP->GetWorldPos();
            cv::Mat pointInRef = refKF->GetRotation() * worldPos + refKF->GetTranslation();
            vPoint->setEstimate(1. / pointInRef.at<float>(2));

            //  original observation
            int indexOfFirstObs = observations.find(refKF)->second;
            vPoint->u0 = refKF->mvKeysUn[indexOfFirstObs].pt.x;
            vPoint->v0 = refKF->mvKeysUn[indexOfFirstObs].pt.y;

            int id = pMP->mnId + maxKFid + 1;
            vPoint->setId(id);
            vPoint->setMarginalized(true);
            optimizer.addVertex(vPoint);

//Set edges
            for (map<KeyFrame *, size_t>::const_iterator mit = observations.begin(), mend = observations.end();
                 mit != mend; mit++) {
                KeyFrame *pKFi = mit->first;

// First observation is treated differently -> it is not added at all
                if (refKF == pKFi)
                    continue;

                if (!pKFi->isBad()) {
                    const cv::KeyPoint &kpUn = pKFi->mvKeysUn[mit->second];

// Monocular observation
                    if (pKFi->mvuRight[mit->second] < 0) {

                        Eigen::Matrix<double, 2, 1> obs;
                        obs << kpUn.pt.x, kpUn.pt.y;

                        g2o::EdgeProjectPSI2UVSingleParam *e = new g2o::EdgeProjectPSI2UVSingleParam();
                        e->resize(3);

// RefKp might not have been added yet
                        g2o::OptimizableGraph::Vertex *v = dynamic_cast<g2o::OptimizableGraph::Vertex *>(optimizer.vertex(
                                refKF->mnId));
                        if (!v) {
                            g2o::VertexSE3Expmap *vSE3 = new g2o::VertexSE3Expmap();
                            vSE3->setEstimate(Converter::toSE3Quat(refKF->GetPose()));
                            vSE3->setId(refKF->mnId);
                            vSE3->setFixed(true);
                            optimizer.addVertex(vSE3);
                        }

                        e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex *>(optimizer.vertex(id)));
                        e->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex *>(optimizer.vertex(pKFi->mnId)));
                        e->setVertex(2, dynamic_cast<g2o::OptimizableGraph::Vertex *>(optimizer.vertex(refKF->mnId)));

                        e->setMeasurement(obs);
                        const float &invSigma2 = pKFi->mvInvLevelSigma2[kpUn.octave];
                        const float assumedDetInvSigma2 = 1.0 / (sigma * sigma);
                        e->setInformation(Eigen::Matrix2d::Identity() * invSigma2 * assumedDetInvSigma2);

                        g2o::RobustKernelHuber *rk = new g2o::RobustKernelHuber;
                        e->setRobustKernel(rk);
                        rk->setDelta(thHuberMono);

                        e->setParameterId(0, 0);

                        g2o::OptimizableGraph::Vertex *v1 = dynamic_cast<g2o::OptimizableGraph::Vertex *>(optimizer.vertex(id));
                        g2o::OptimizableGraph::Vertex *v2 = dynamic_cast<g2o::OptimizableGraph::Vertex *>(optimizer.vertex(pKFi->mnId));
                        g2o::OptimizableGraph::Vertex *v3 = dynamic_cast<g2o::OptimizableGraph::Vertex *>(optimizer.vertex(refKF->mnId));
                        if (!v1 || !v2 || !v3)
                            std::cout << "InvDepth: Vertices should exists but at least one doesn't. Why? Don't know yet" << std::endl;
                        else
                            optimizer.addEdge(e);
                        vpEdgesMono.push_back(e);
                        vpEdgeKFMono.push_back(pKFi);
                        vpMapPointEdgeMono.push_back(pMP);

                    }
                }
            }
        }

        if (pbStopFlag)
            if (*pbStopFlag)
                return std::vector<double>(); // Modified


        double initchi2 = 0;
        int count = 0;
        for (size_t i = 0, iend = vpEdgesMono.size(); i < iend; i++) {
            g2o::EdgeProjectPSI2UVSingleParam *e = vpEdgesMono[i];
            MapPoint *pMP = vpMapPointEdgeMono[i];

            e->computeError();
            if (pMP->isBad())
                continue;

            if (e->chi2() <= 5.991 && e->isDepthPositive()) {
                initchi2 += e->chi2();
                count++;
            }
        }
        std::cout << "\tLocalBA: avg initial chi2 = " << initchi2 / count << " over " << count << " measurements "
                  << std::endl;


        auto start = chrono::steady_clock::now();
        optimizer.initializeOptimization(0);
        auto end = chrono::steady_clock::now();
        auto diff = end - start;

//    std::cout << "!!!!!!!!! ---------- !!!!!!!!!! " << std::endl;
        auto start2 = chrono::steady_clock::now();
        optimizer.optimize(5);
        auto end2 = chrono::steady_clock::now();
        auto diff2 = end2 - start2;
//    std::cout << "!!!!!!!!! ---------- !!!!!!!!!! " << std::endl;

        auto diff3 = diff2, diff4 = diff2, diff5 = diff2;

        bool bDoMore = true;

        if (pbStopFlag)
            if (*pbStopFlag)
                bDoMore = false;

        if (bDoMore) {

// Check inlier observations
            for (size_t i = 0, iend = vpEdgesMono.size(); i < iend; i++) {
                g2o::EdgeProjectPSI2UVSingleParam *e = vpEdgesMono[i];
                MapPoint *pMP = vpMapPointEdgeMono[i];

                if (pMP->isBad())
                    continue;

                if (e->chi2() > 5.991 || !e->isDepthPositive()) {
                    e->setLevel(1);
                }

                e->setRobustKernel(0);
            }

// Optimize again without the outliers

//        std::cout << " ------- " << std::endl;

            auto start3 = chrono::steady_clock::now();
            optimizer.initializeOptimization(0);
            auto end3 = chrono::steady_clock::now();
            diff3 = end3 - start3;

//        std::cout << " ------- " << std::endl;

            auto start4 = chrono::steady_clock::now();
            optimizer.optimize(5);
            auto end4 = chrono::steady_clock::now();
            diff4 = end4 - start4;

//        std::cout << " ------- " << std::endl;

            auto start5 = chrono::steady_clock::now();
            optimizer.optimize(5);
            auto end5 = chrono::steady_clock::now();
            diff5 = end5 - start5;

//        std::cout << " ------- " << std::endl;

        }

        vector<pair<KeyFrame *, MapPoint *> > vToErase;
        vToErase.reserve(vpEdgesMono.size());

// Check inlier observations
        std::vector<double> chi2Statistics;
        for (size_t i = 0, iend = vpEdgesMono.size(); i < iend; i++) {
            g2o::EdgeProjectPSI2UVSingleParam *e = vpEdgesMono[i];
            MapPoint *pMP = vpMapPointEdgeMono[i];

            if (pMP->isBad())
                continue;

            if (e->chi2() > 5.991 || !e->isDepthPositive()) {
                KeyFrame *pKFi = vpEdgeKFMono[i];
                vToErase.push_back(make_pair(pKFi, pMP));
            } else {
//            double * err = e->errorData();
//            double * inf = e->informationData();
//            double eDist = std::sqrt(err[0] * err[0] + err[1] * err[1]);
//            std::cout << "ERR: " <<  err[0] << " " <<  err[1] << " ReprojError: " << eDist <<" Chi2: " << e->chi2() << std::endl;
                chi2Statistics.push_back(e->chi2());
            }
        }
        std::cout << "\tLocalBA: avg chi2() = "
                  << accumulate(chi2Statistics.begin(), chi2Statistics.end(), 0.0) / chi2Statistics.size() << " over "
                  << chi2Statistics.size() << " measurements" << std::endl;


// Get Map Mutex
        unique_lock<mutex> lock(pMap->mMutexMapUpdate);

        if (!vToErase.empty()) {
            for (size_t i = 0; i < vToErase.size(); i++) {
                KeyFrame *pKFi = vToErase[i].first;
                MapPoint *pMPi = vToErase[i].second;
                pKFi->EraseMapPointMatch(pMPi);
                pMPi->EraseObservation(pKFi);
            }
        }

// Recover optimized data

//Keyframes
        for (list<KeyFrame *>::iterator lit = lLocalKeyFrames.begin(), lend = lLocalKeyFrames.end(); lit != lend; lit++) {
            KeyFrame *pKF = *lit;
            g2o::VertexSE3Expmap *vSE3 = static_cast<g2o::VertexSE3Expmap *>(optimizer.vertex(pKF->mnId));
            g2o::SE3Quat SE3quat = vSE3->estimate();
            pKF->SetPose(Converter::toCvMat(SE3quat));
        }

//Points
        for (set<MapPoint *>::iterator lit = lLocalMapPoints.begin(), lend = lLocalMapPoints.end(); lit != lend; lit++) {
            MapPoint *pMP = *lit;

            g2o::VertexSBAPointInvD *vPoint = static_cast<g2o::VertexSBAPointInvD *>(optimizer.vertex(
                    pMP->mnId + maxKFid + 1));
            KeyFrame *refKF = pMP->GetReferenceKeyFrame();

            cv::Mat pointInFirst = cv::Mat(3,1, CV_32F), worldPos = cv::Mat(3,1, CV_32F);
            pointInFirst.at<float>(2) = 1. / vPoint->estimate();
            pointInFirst.at<float>(0) = (vPoint->u0 - refKF->cx) * pointInFirst.at<float>(2) / refKF->fx;
            pointInFirst.at<float>(1) = (vPoint->v0 - refKF->cy) * pointInFirst.at<float>(2) / refKF->fy;


            cv::Mat worldPointPos = refKF->GetRotation().t() * pointInFirst - refKF->GetRotation().t() * refKF->GetTranslation();

            pMP->SetWorldPos(worldPointPos);
            pMP->UpdateNormalAndDepth();
        }


        // We verify the reprojection error
        OptimizerBA::VerifyReprojectionError(lLocalKeyFrames, lLocalMapPoints, lFixedCameras);


        cout << "\t LBA InvD Single: Initialization time: " << chrono::duration <double, milli> (diff).count() << " "<<  chrono::duration <double, milli> (diff3).count() << " ms" << endl;
        cout << "\t LBA InvD Single: Pure optimization time: " << chrono::duration <double, milli> (diff2).count() << " " << chrono::duration <double, milli> (diff4).count()
             << " " << chrono::duration <double, milli> (diff5).count() << " ms" << endl;


        return chi2Statistics;
    }




    std::vector<double> OptimizerBA::BundleAdjustmentInvDepthSingleParamPatch(list<KeyFrame *> lLocalKeyFrames, list<KeyFrame *> lFixedCameras, set<MapPoint *> lLocalMapPoints, bool *pbStopFlag, Map *pMap, float sigma)
    {
        const float thHuber = 9;
        const int blockSolverCameras = 6;
        const int blockSolverPoses = 1;

        // Setup optimizer
        g2o::SparseOptimizer optimizer;

        typedef g2o::BlockSolver< g2o::BlockSolverTraits<blockSolverCameras, blockSolverPoses> > OurBlockSolver;
        OurBlockSolver::LinearSolverType *linearSolver;
        linearSolver = new g2o::LinearSolverEigen<OurBlockSolver::PoseMatrixType>();
        OurBlockSolver *solver_ptr = new OurBlockSolver(linearSolver);

        g2o::OptimizationAlgorithmLevenberg *solver = new g2o::OptimizationAlgorithmLevenberg(solver_ptr);
        optimizer.setAlgorithm(solver);

        if (pbStopFlag)
            optimizer.setForceStopFlag(pbStopFlag);

        KeyFrame *pKF = *lLocalKeyFrames.begin(); // TODO: Now setting it as camera params of 1st camera
        Eigen::Vector2d principal_point(pKF->cx, pKF->cy);
        g2o::CameraParameters *cam_params
                = new g2o::CameraParameters(pKF->fx, pKF->fy, principal_point, 0.);
        cam_params->setId(0);

        if (!optimizer.addParameter(cam_params)) {
            assert(false);
        }

        unsigned long maxKFid = 0;

// Set Local KeyFrame vertices
        for (list<KeyFrame *>::iterator lit = lLocalKeyFrames.begin(), lend = lLocalKeyFrames.end(); lit != lend; lit++) {
            KeyFrame *pKFi = *lit;
            g2o::VertexSE3Expmap *vSE3 = new g2o::VertexSE3Expmap();
            vSE3->setEstimate(Converter::toSE3Quat(pKFi->GetPose()));
            vSE3->setId(pKFi->mnId);
            vSE3->setFixed(pKFi->mnId == 0);
            optimizer.addVertex(vSE3);
            if (pKFi->mnId > maxKFid)
                maxKFid = pKFi->mnId;
        }

// Set Fixed KeyFrame vertices
        for (list<KeyFrame *>::iterator lit = lFixedCameras.begin(), lend = lFixedCameras.end(); lit != lend; lit++) {
            KeyFrame *pKFi = *lit;
            g2o::VertexSE3Expmap *vSE3 = new g2o::VertexSE3Expmap();
            vSE3->setEstimate(Converter::toSE3Quat(pKFi->GetPose()));
            vSE3->setId(pKFi->mnId);
            vSE3->setFixed(true);
            optimizer.addVertex(vSE3);
            if (pKFi->mnId > maxKFid)
                maxKFid = pKFi->mnId;
        }

// Set MapPoint vertices
        const int nExpectedSize = (lLocalKeyFrames.size() + lFixedCameras.size()) * lLocalMapPoints.size();

        vector<g2o::EdgeProjectPSI2UVSingleParamPatch *> vpEdgesMono;
        vpEdgesMono.reserve(nExpectedSize);

        vector<KeyFrame *> vpEdgeKFMono;
        vpEdgeKFMono.reserve(nExpectedSize);

        vector<MapPoint *> vpMapPointEdgeMono;
        vpMapPointEdgeMono.reserve(nExpectedSize);


        for (set<MapPoint *>::iterator lit = lLocalMapPoints.begin(), lend = lLocalMapPoints.end(); lit != lend; lit++) {
            MapPoint *pMP = *lit;
            KeyFrame *refKF = pMP->GetReferenceKeyFrame();
            const map<KeyFrame *, size_t> observations = pMP->GetObservations();


            g2o::VertexSBAPointInvD *vPoint = new g2o::VertexSBAPointInvD();

// Moving from world coordinate to local coordinate in ref kf and then saving depth as inverse Z
            cv::Mat worldPos = pMP->GetWorldPos();
            cv::Mat pointInRef = refKF->GetRotation() * worldPos + refKF->GetTranslation();
            vPoint->setEstimate(1. / pointInRef.at<float>(2));

//  original observation
            int indexOfFirstObs = observations.find(refKF)->second;
            vPoint->u0 = refKF->mvKeysUn[indexOfFirstObs].pt.x;
            vPoint->v0 = refKF->mvKeysUn[indexOfFirstObs].pt.y;
            const int level = refKF->mvKeysUn[indexOfFirstObs].octave;
            const double refScaleFactor = refKF->mvScaleFactors[level];
            cv::Mat refPatch = refKF->mPatches[indexOfFirstObs];

            std::vector<double> refPatchVec;
            for (int y = 0; y < refPatch.rows; y++)
                for (int x = 0; x < refPatch.cols; x++)
                    refPatchVec.push_back((double) refPatch.at<uchar>(y, x));


            int id = pMP->mnId + maxKFid + 1;
            vPoint->setId(id);
            vPoint->setMarginalized(true);
            optimizer.addVertex(vPoint);

//Set edges
            for (map<KeyFrame *, size_t>::const_iterator mit = observations.begin(), mend = observations.end();
                 mit != mend; mit++) {
                KeyFrame *pKFi = mit->first;

// First observation is treated differently -> it is not added at all
                if (refKF == pKFi)
                    continue;

                if (!pKFi->isBad()) {
                    const cv::KeyPoint &kpUn = pKFi->mvKeysUn[mit->second];
                    cv::Mat curPatch = pKFi->mPatches[mit->second];

                    std::vector<double> curPatchVec;
                    for (int y = 0; y < curPatch.rows; y++)
                        for (int x = 0; x < curPatch.cols; x++)
                            curPatchVec.push_back((double) curPatch.at<uchar>(y, x));

// Monocular observation
                    if (pKFi->mvuRight[mit->second] < 0) {

                        Eigen::Matrix<double, 9, 1> obs;
                        obs << kpUn.pt.x, kpUn.pt.y, 0, 0, 0, 0, 0, 0, 0;

                        g2o::EdgeProjectPSI2UVSingleParamPatch *e = new g2o::EdgeProjectPSI2UVSingleParamPatch();
                        e->resize(3);

// RefKp might not have been added yet
                        g2o::OptimizableGraph::Vertex *v = dynamic_cast<g2o::OptimizableGraph::Vertex *>(optimizer.vertex(
                                refKF->mnId));
                        if (!v) {
                            g2o::VertexSE3Expmap *vSE3 = new g2o::VertexSE3Expmap();
                            vSE3->setEstimate(Converter::toSE3Quat(refKF->GetPose()));
                            vSE3->setId(refKF->mnId);
                            vSE3->setFixed(true);
                            optimizer.addVertex(vSE3);
                        }

                        e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex *>(optimizer.vertex(id)));
                        e->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex *>(optimizer.vertex(pKFi->mnId)));
                        e->setVertex(2, dynamic_cast<g2o::OptimizableGraph::Vertex *>(optimizer.vertex(refKF->mnId)));


// TODO: Setting the necessary values
                        const double curScaleFactor = pKFi->mvScaleFactors[kpUn.octave];
                        e->setAdditionalData(refPatchVec, curPatchVec, refScaleFactor, curScaleFactor);


                        e->setMeasurement(obs);
                        const float &invSigma2 = pKFi->mvInvLevelSigma2[kpUn.octave];
                        const float assumedDetInvSigma2 = 1.0 / (sigma * sigma);
//                    e->setInformation(Eigen::Matrix2d::Identity() * invSigma2 * assumedDetInvSigma2);
                        e->setInformation(Eigen::Matrix<double, 9, 9>::Identity() * 0.01);

                        g2o::RobustKernelHuber *rk = new g2o::RobustKernelHuber;
                        e->setRobustKernel(rk);
                        rk->setDelta(thHuber);

                        e->setParameterId(0, 0);


                        g2o::OptimizableGraph::Vertex *v1 = dynamic_cast<g2o::OptimizableGraph::Vertex *>(optimizer.vertex(id));
                        g2o::OptimizableGraph::Vertex *v2 = dynamic_cast<g2o::OptimizableGraph::Vertex *>(optimizer.vertex(pKFi->mnId));
                        g2o::OptimizableGraph::Vertex *v3 = dynamic_cast<g2o::OptimizableGraph::Vertex *>(optimizer.vertex(refKF->mnId));
                        if (!v1 || !v2 || !v3)
                            std::cout << "Vertices should exists but at least one doesn't. Why? Don't know yet" << std::endl;
                        else
                            optimizer.addEdge(e);
                        vpEdgesMono.push_back(e);
                        vpEdgeKFMono.push_back(pKFi);
                        vpMapPointEdgeMono.push_back(pMP);

                    }
                }
            }
        }

        if (pbStopFlag)
            if (*pbStopFlag)
                return std::vector<double>(); // Modified

        auto start = chrono::steady_clock::now();
        optimizer.initializeOptimization(0);
        auto end = chrono::steady_clock::now();
        auto diff = end - start;


        double initchi2 = 0;
        int count = 0;
        for (size_t i = 0, iend = vpEdgesMono.size(); i < iend; i++) {
            g2o::EdgeProjectPSI2UVSingleParamPatch *e = vpEdgesMono[i];
            MapPoint *pMP = vpMapPointEdgeMono[i];

            e->computeError();
            if (pMP->isBad())
                continue;

            if (e->chi2() <= thHuber && e->isDepthPositive()) {
//            std::cout << " ??? " << e->chi2() << std::endl;
                initchi2 += e->chi2();
                count ++;
            }
        }
        std::cout << "\tLocalBA Patches: avg initial chi2 = " << initchi2 / count<< " over " << count << " measurements " << std::endl;



//    std::cout << "!!!!!!!!! ---------- !!!!!!!!!! " << std::endl;
        auto start2 = chrono::steady_clock::now();
        optimizer.optimize(5);
        auto end2 = chrono::steady_clock::now();
        auto diff2 = end2 - start2;
//    std::cout << "!!!!!!!!! ---------- !!!!!!!!!! " << std::endl;

        auto diff3 = diff2, diff4 = diff2, diff5 = diff2;

        bool bDoMore = true;

        if (pbStopFlag)
            if (*pbStopFlag)
                bDoMore = false;

        if (bDoMore) {

// Check inlier observations
            for (size_t i = 0, iend = vpEdgesMono.size(); i < iend; i++) {
                g2o::EdgeProjectPSI2UVSingleParamPatch *e = vpEdgesMono[i];
                MapPoint *pMP = vpMapPointEdgeMono[i];

                if (pMP->isBad())
                    continue;

                if (e->chi2() > thHuber || !e->isDepthPositive()) {
                    e->setLevel(1);
                }

                e->setRobustKernel(0);
            }

// Optimize again without the outliers
//        std::cout << " ------- " << std::endl;

            auto start3 = chrono::steady_clock::now();
            optimizer.initializeOptimization(0);
            auto end3 = chrono::steady_clock::now();
            diff3 = end3 - start3;

//        std::cout << " ------- " << std::endl;

            auto start4 = chrono::steady_clock::now();
            optimizer.optimize(5);
            auto end4 = chrono::steady_clock::now();
            diff4 = end4 - start4;

//        std::cout << " ------- " << std::endl;

            auto start5 = chrono::steady_clock::now();
            optimizer.optimize(5);
            auto end5 = chrono::steady_clock::now();
            diff5 = end5 - start5;

//        std::cout << " ------- " << std::endl;

        }

        vector<pair<KeyFrame *, MapPoint *> > vToErase;
        vToErase.reserve(vpEdgesMono.size()); //+ vpEdgesStereo.size());

// Check inlier observations
        std::vector<double> chi2Statistics;
        for (size_t i = 0, iend = vpEdgesMono.size(); i < iend; i++) {
            g2o::EdgeProjectPSI2UVSingleParamPatch *e = vpEdgesMono[i];
            MapPoint *pMP = vpMapPointEdgeMono[i];

            if (pMP->isBad())
                continue;

            if (e->chi2() > thHuber || !e->isDepthPositive()) {
                KeyFrame *pKFi = vpEdgeKFMono[i];
                vToErase.push_back(make_pair(pKFi, pMP));
            } else {
//            double * err = e->errorData();
//            double * inf = e->informationData();
//            double eDist = std::sqrt(err[0] * err[0] + err[1] * err[1]);
//            std::cout << "ERR: " <<  err[0] << " " <<  err[1] << " ReprojError: " << eDist <<" Chi2: " << e->chi2() << std::endl;
                chi2Statistics.push_back(e->chi2());
            }
        }
        std::cout << "\tLocalBA Patches: avg chi2() = "
                  << accumulate(chi2Statistics.begin(), chi2Statistics.end(), 0.0) / chi2Statistics.size() << " over "
                  << chi2Statistics.size() << " measurements" << std::endl;



// Get Map Mutex
        unique_lock<mutex> lock(pMap->mMutexMapUpdate);

        if (!vToErase.empty()) {
            for (size_t i = 0; i < vToErase.size(); i++) {
                KeyFrame *pKFi = vToErase[i].first;
                MapPoint *pMPi = vToErase[i].second;
                pKFi->EraseMapPointMatch(pMPi);
                pMPi->EraseObservation(pKFi);
            }
        }

// Recover optimized data

//Keyframes
        for (list<KeyFrame *>::iterator lit = lLocalKeyFrames.begin(), lend = lLocalKeyFrames.end(); lit != lend; lit++) {
            KeyFrame *pKF = *lit;
            g2o::VertexSE3Expmap *vSE3 = static_cast<g2o::VertexSE3Expmap *>(optimizer.vertex(pKF->mnId));
            g2o::SE3Quat SE3quat = vSE3->estimate();
            pKF->SetPose(Converter::toCvMat(SE3quat));
        }

//Points
        for (set<MapPoint *>::iterator lit = lLocalMapPoints.begin(), lend = lLocalMapPoints.end(); lit != lend; lit++) {
            MapPoint *pMP = *lit;

            g2o::VertexSBAPointInvD *vPoint = static_cast<g2o::VertexSBAPointInvD *>(optimizer.vertex(
                    pMP->mnId + maxKFid + 1));
            KeyFrame *refKF = pMP->GetReferenceKeyFrame();

            cv::Mat pointInFirst = cv::Mat(3, 1, CV_32F), worldPos = cv::Mat(3, 1, CV_32F);
            pointInFirst.at<float>(2) = 1. / vPoint->estimate();
            pointInFirst.at<float>(0) = (vPoint->u0 - refKF->cx) * pointInFirst.at<float>(2) / refKF->fx;
            pointInFirst.at<float>(1) = (vPoint->v0 - refKF->cy) * pointInFirst.at<float>(2) / refKF->fy;


            cv::Mat worldPointPos =
                    refKF->GetRotation().t() * pointInFirst - refKF->GetRotation().t() * refKF->GetTranslation();

            pMP->SetWorldPos(worldPointPos);
            pMP->UpdateNormalAndDepth();
        }


        // We verify the reprojection error
        OptimizerBA::VerifyReprojectionError(lLocalKeyFrames, lLocalMapPoints, lFixedCameras);

        cout << "\t LBA InvD Patches: Initialization time: " << chrono::duration <double, milli> (diff).count() << " "<<  chrono::duration <double, milli> (diff3).count() << " ms" << endl;
        cout << "\t LBA InvD Patches: Pure optimization time: " << chrono::duration <double, milli> (diff2).count() << " " << chrono::duration <double, milli> (diff4).count()
             << " " << chrono::duration <double, milli> (diff5).count() << " ms" << endl;
        return chi2Statistics;
    }


    std::vector<double> OptimizerBA::BundleAdjustmentInvDepthSingleParamPatchBright(list<KeyFrame *> lLocalKeyFrames, list<KeyFrame *> lFixedCameras, set<MapPoint *> lLocalMapPoints,  bool *pbStopFlag, Map *pMap, float sigma) {

        const float thHuber = 81; // as in the DSO
        const float sqrtThHuber = sqrt(thHuber);
        const int blockSolverCameras = 8;
        const int blockSolverPoses = 1;

        // Setup optimizer
        g2o::SparseOptimizer optimizer;

        typedef g2o::BlockSolver< g2o::BlockSolverTraits<blockSolverCameras, blockSolverPoses> > OurBlockSolver;
        OurBlockSolver::LinearSolverType *linearSolver;
        linearSolver = new g2o::LinearSolverEigen<OurBlockSolver::PoseMatrixType>();
        OurBlockSolver *solver_ptr = new OurBlockSolver(linearSolver);

        g2o::OptimizationAlgorithmLevenberg *solver = new g2o::OptimizationAlgorithmLevenberg(solver_ptr);
        optimizer.setAlgorithm(solver);

        if (pbStopFlag)
            optimizer.setForceStopFlag(pbStopFlag);

        KeyFrame *pKF = *lLocalKeyFrames.begin(); // TODO: Now setting it as camera params of 1st camera
        Eigen::Vector2d principal_point(pKF->cx, pKF->cy);
        g2o::CameraParameters *cam_params
                = new g2o::CameraParameters(pKF->fx, pKF->fy, principal_point, 0.);
        cam_params->setId(0);

        if (!optimizer.addParameter(cam_params)) {
            assert(false);
        }

        unsigned long maxKFid = 0;

// Set Local KeyFrame vertices
        for (list<KeyFrame *>::iterator lit = lLocalKeyFrames.begin(), lend = lLocalKeyFrames.end(); lit != lend; lit++) {
            KeyFrame *pKFi = *lit;
            g2o::VertexSE3ExpmapBright *vSE3 = new g2o::VertexSE3ExpmapBright();
            g2o::SE3QuatBright estimate;
            estimate.se3quat = Converter::toSE3Quat(pKFi->GetPose());
            estimate.a = pKFi->affineA;
            estimate.b = pKFi->affineB;

            vSE3->setEstimate(estimate);
            vSE3->setId(pKFi->mnId);
            vSE3->setFixed(pKFi->mnId == 0);


            optimizer.addVertex(vSE3);
            if (pKFi->mnId > maxKFid)
                maxKFid = pKFi->mnId;




//        std::cout << "Starting poses: " << vSE3->id() << " a=" << estimate.a << " b=" <<estimate.b<< std::endl << estimate.se3quat << std::endl;
        }

// Set Fixed KeyFrame vertices
        for (list<KeyFrame *>::iterator lit = lFixedCameras.begin(), lend = lFixedCameras.end(); lit != lend; lit++) {
            KeyFrame *pKFi = *lit;
            g2o::VertexSE3ExpmapBright *vSE3 = new g2o::VertexSE3ExpmapBright();
            g2o::SE3QuatBright estimate;
            estimate.se3quat = Converter::toSE3Quat(pKFi->GetPose());
            estimate.a = pKFi->affineA;
            estimate.b = pKFi->affineB;

            vSE3->setEstimate(estimate);
            vSE3->setId(pKFi->mnId);
            vSE3->setFixed(true);
            optimizer.addVertex(vSE3);
            if (pKFi->mnId > maxKFid)
                maxKFid = pKFi->mnId;
        }

// Set MapPoint vertices
        const int nExpectedSize = (lLocalKeyFrames.size() + lFixedCameras.size()) * lLocalMapPoints.size();

        vector<g2o::EdgeProjectPSI2UVSingleParamPatchBright *> vpEdgesMono;
        vpEdgesMono.reserve(nExpectedSize);

        vector<KeyFrame *> vpEdgeKFMono;
        vpEdgeKFMono.reserve(nExpectedSize);

        vector<MapPoint *> vpMapPointEdgeMono;
        vpMapPointEdgeMono.reserve(nExpectedSize);



        for (set<MapPoint *>::iterator lit = lLocalMapPoints.begin(), lend = lLocalMapPoints.end(); lit != lend; lit++) {
            MapPoint *pMP = *lit;
            KeyFrame *refKF = pMP->GetReferenceKeyFrame();
            const map<KeyFrame *, size_t> observations = pMP->GetObservations();


            g2o::VertexSBAPointInvD *vPoint = new g2o::VertexSBAPointInvD();

// Moving from world coordinate to local coordinate in ref kf and then saving depth as inverse Z
            cv::Mat worldPos = pMP->GetWorldPos();
            cv::Mat pointInRef = refKF->GetRotation() * worldPos + refKF->GetTranslation();
            vPoint->setEstimate(1. / pointInRef.at<float>(2));

//  original observation
            int indexOfFirstObs = observations.find(refKF)->second;
            vPoint->u0 = refKF->mvKeysUn[indexOfFirstObs].pt.x;
            vPoint->v0 = refKF->mvKeysUn[indexOfFirstObs].pt.y;
            const int level = refKF->mvKeysUn[indexOfFirstObs].octave;
            const double refScaleFactor = refKF->mvScaleFactors[level];
            cv::Mat refPatch = refKF->mPatches[indexOfFirstObs];

            std::vector<double> refPatchVec;
            for (int y = 0; y < refPatch.rows; y++)
                for (int x = 0; x < refPatch.cols; x++)
                    refPatchVec.push_back((double) refPatch.at<uchar>(y, x));


            int id = pMP->mnId + maxKFid + 1;
            vPoint->setId(id);
            vPoint->setMarginalized(true);
            optimizer.addVertex(vPoint);

//Set edges
            for (map<KeyFrame *, size_t>::const_iterator mit = observations.begin(), mend = observations.end();
                 mit != mend; mit++) {
                KeyFrame *pKFi = mit->first;

// First observation is treated differently -> it is not added at all
                if (refKF == pKFi)
                    continue;

                if (!pKFi->isBad()) {
                    const cv::KeyPoint &kpUn = pKFi->mvKeysUn[mit->second];
                    cv::Mat curPatch = pKFi->mPatches[mit->second];

                    if ( curPatch.empty() )
                        continue;

                    std::vector<double> curPatchVec;
                    for (int y = 0; y < curPatch.rows; y++)
                        for (int x = 0; x < curPatch.cols; x++)
                            curPatchVec.push_back((double) curPatch.at<uchar>(y, x));

// Monocular observation
                    if (pKFi->mvuRight[mit->second] < 0) {

                        Eigen::Matrix<double, 9, 1> obs;
                        obs << kpUn.pt.x, kpUn.pt.y, 0, 0, 0, 0, 0, 0, 0;

                        g2o::EdgeProjectPSI2UVSingleParamPatchBright *e = new g2o::EdgeProjectPSI2UVSingleParamPatchBright();
                        e->resize(3);

// RefKp might not have been added yet
                        g2o::OptimizableGraph::Vertex *v = dynamic_cast<g2o::OptimizableGraph::Vertex *>(optimizer.vertex(
                                refKF->mnId));
                        if (!v) {
                            g2o::VertexSE3ExpmapBright *vSE3 = new g2o::VertexSE3ExpmapBright();

                            g2o::SE3QuatBright estimate;
                            estimate.se3quat = Converter::toSE3Quat(refKF->GetPose());
                            estimate.a = refKF->affineA;
                            estimate.b = refKF->affineB;
                            vSE3->setEstimate(estimate);

                            vSE3->setId(refKF->mnId);
                            vSE3->setFixed(true);
                            optimizer.addVertex(vSE3);
                        }

                        e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex *>(optimizer.vertex(id)));
                        e->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex *>(optimizer.vertex(pKFi->mnId)));
                        e->setVertex(2, dynamic_cast<g2o::OptimizableGraph::Vertex *>(optimizer.vertex(refKF->mnId)));


                    //  Setting the necessary values
                        const double curScaleFactor = pKFi->mvScaleFactors[kpUn.octave];
                        e->setAdditionalData(refPatchVec, curPatchVec, refScaleFactor, curScaleFactor);


                        e->setMeasurement(obs);

                        // TODO: It is possible to weight the error by the pyramind lvl of point but does it make sense?
                        const float &invSigma2 = pKFi->mvInvLevelSigma2[kpUn.octave];
                        const float assumedDetInvSigma2 = 1.0 / (sigma * sigma);

                        e->setInformation(Eigen::Matrix<double, 9, 9>::Identity());

                        g2o::RobustKernelHuber *rk = new g2o::RobustKernelHuber;
                        e->setRobustKernel(rk);
                        rk->setDelta(sqrtThHuber);

                        e->setParameterId(0, 0);


                        g2o::OptimizableGraph::Vertex *v1 = dynamic_cast<g2o::OptimizableGraph::Vertex *>(optimizer.vertex(id));
                        g2o::OptimizableGraph::Vertex *v2 = dynamic_cast<g2o::OptimizableGraph::Vertex *>(optimizer.vertex(pKFi->mnId));
                        g2o::OptimizableGraph::Vertex *v3 = dynamic_cast<g2o::OptimizableGraph::Vertex *>(optimizer.vertex(refKF->mnId));

                        if (!v1 || !v2 || !v3)
                            std::cout << "Vertices should exists but at least one doesn't. Why? Don't know yet" << std::endl;
                        else
                            optimizer.addEdge(e);
                        vpEdgesMono.push_back(e);
                        vpEdgeKFMono.push_back(pKFi);
                        vpMapPointEdgeMono.push_back(pMP);

                    }
                }
            }
        }

        if (pbStopFlag)
            if (*pbStopFlag)
                return std::vector<double>(); // Modified

        auto start = chrono::steady_clock::now();
        optimizer.initializeOptimization(0);
        auto end = chrono::steady_clock::now();
        auto diff = end - start;


        double initchi2 = 0;
        int count = 0;
        for (size_t i = 0, iend = vpEdgesMono.size(); i < iend; i++) {
            g2o::EdgeProjectPSI2UVSingleParamPatchBright *e = vpEdgesMono[i];
            MapPoint *pMP = vpMapPointEdgeMono[i];

            e->computeError();
            if (pMP->isBad())
                continue;

            if (e->chi2() <= thHuber && e->isDepthPositive()) {
//            std::cout << " ??? " << e->chi2() << std::endl;
                initchi2 += e->chi2();
                count++;
            }
        }
        std::cout << "\tLBA InvD Patches: avg initial chi2 = " << initchi2 / count << " over " << count << " measurements "
                  << std::endl;



//    std::cout << "!!!!!!!!! ---------- !!!!!!!!!! " << std::endl;
        auto start2 = chrono::steady_clock::now();
        optimizer.optimize(10);
        auto end2 = chrono::steady_clock::now();
        auto diff2 = end2 - start2;
//    std::cout << "!!!!!!!!! ---------- !!!!!!!!!! " << std::endl;

        auto diff3 = diff2, diff4 = diff2, diff5 = diff2;

        bool bDoMore = false; // We won't do more in that case

        if (pbStopFlag)
            if (*pbStopFlag)
                bDoMore = false;

        if (bDoMore) {

// Check inlier observations
            for (size_t i = 0, iend = vpEdgesMono.size(); i < iend; i++) {
                g2o::EdgeProjectPSI2UVSingleParamPatchBright *e = vpEdgesMono[i];
                MapPoint *pMP = vpMapPointEdgeMono[i];

                if (pMP->isBad())
                    continue;

                if (e->chi2() > thHuber || !e->isDepthPositive()) {
                    e->setLevel(1);
                }

                e->setRobustKernel(0);
            }



// Optimize again without the outliers
//        std::cout << " ------- " << std::endl;

            auto start3 = chrono::steady_clock::now();
            optimizer.initializeOptimization(0);
            auto end3 = chrono::steady_clock::now();
            diff3 = end3 - start3;

//        std::cout << " ------- " << std::endl;

            auto start4 = chrono::steady_clock::now();
            optimizer.optimize(5);
            auto end4 = chrono::steady_clock::now();
            diff4 = end4 - start4;

//        std::cout << " ------- " << std::endl;

            auto start5 = chrono::steady_clock::now();
            optimizer.optimize(5);
            auto end5 = chrono::steady_clock::now();
            diff5 = end5 - start5;

//        std::cout << " ------- " << std::endl;

        }

        vector<pair<KeyFrame *, MapPoint *> > vToErase;
        vToErase.reserve(vpEdgesMono.size()); //+ vpEdgesStereo.size());

// Check inlier observations
        std::vector<double> chi2Statistics;
        for (size_t i = 0, iend = vpEdgesMono.size(); i < iend; i++) {
            g2o::EdgeProjectPSI2UVSingleParamPatchBright *e = vpEdgesMono[i];
            MapPoint *pMP = vpMapPointEdgeMono[i];

            if (pMP->isBad())
                continue;

            if (e->chi2() > thHuber || !e->isDepthPositive()) {
                KeyFrame *pKFi = vpEdgeKFMono[i];
                vToErase.push_back(make_pair(pKFi, pMP));
            } else {
//            double * err = e->errorData();
//            double * inf = e->informationData();
//            double eDist = std::sqrt(err[0] * err[0] + err[1] * err[1]);
//            std::cout << "ERR: " <<  err[0] << " " <<  err[1] << " ReprojError: " << eDist <<" Chi2: " << e->chi2() << std::endl;
                chi2Statistics.push_back(e->chi2());
            }
        }
        std::cout << "\tLBA InvD Patches:: avg chi2() = "
                  << accumulate(chi2Statistics.begin(), chi2Statistics.end(), 0.0) / chi2Statistics.size() << " over "
                  << chi2Statistics.size() << " measurements" << std::endl;


// Get Map Mutex
        unique_lock<mutex> lock(pMap->mMutexMapUpdate);

        if (!vToErase.empty()) {
            for (size_t i = 0; i < vToErase.size(); i++) {
                KeyFrame *pKFi = vToErase[i].first;
                MapPoint *pMPi = vToErase[i].second;
                pKFi->EraseMapPointMatch(pMPi);
                pMPi->EraseObservation(pKFi);
            }
        }

// Recover optimized data

//Keyframes
        for (list<KeyFrame *>::iterator lit = lLocalKeyFrames.begin(), lend = lLocalKeyFrames.end(); lit != lend; lit++) {
            KeyFrame *pKF = *lit;
            g2o::VertexSE3ExpmapBright *vSE3 = static_cast<g2o::VertexSE3ExpmapBright *>(optimizer.vertex(pKF->mnId));
            g2o::SE3QuatBright est = vSE3->estimate();
            g2o::SE3Quat SE3quat = est.se3quat;

//            std::cout << "Retrieved poses: " << vSE3->id() << " a=" << est.a << " b=" <<est.b<< std::endl << est.se3quat << std::endl;
            pKF->SetPose(Converter::toCvMat(SE3quat));
            pKF->affineA = est.a;
            pKF->affineB = est.b;
        }

//Points
        for (set<MapPoint *>::iterator lit = lLocalMapPoints.begin(), lend = lLocalMapPoints.end(); lit != lend; lit++) {
            MapPoint *pMP = *lit;

            g2o::VertexSBAPointInvD *vPoint = static_cast<g2o::VertexSBAPointInvD *>(optimizer.vertex(
                    pMP->mnId + maxKFid + 1));
            KeyFrame *refKF = pMP->GetReferenceKeyFrame();

            cv::Mat pointInFirst = cv::Mat(3, 1, CV_32F), worldPos = cv::Mat(3, 1, CV_32F);
            pointInFirst.at<float>(2) = 1. / vPoint->estimate();
            pointInFirst.at<float>(0) = (vPoint->u0 - refKF->cx) * pointInFirst.at<float>(2) / refKF->fx;
            pointInFirst.at<float>(1) = (vPoint->v0 - refKF->cy) * pointInFirst.at<float>(2) / refKF->fy;


            cv::Mat worldPointPos =
                    refKF->GetRotation().t() * pointInFirst - refKF->GetRotation().t() * refKF->GetTranslation();

            pMP->SetWorldPos(worldPointPos);
            pMP->UpdateNormalAndDepth();
        }


        // We verify the reprojection error
        OptimizerBA::VerifyReprojectionError(lLocalKeyFrames, lLocalMapPoints, lFixedCameras);


        cout << "\t LBA InvD Patches: Pure optimization time: " << chrono::duration<double, milli>(diff2).count() << " "
             << chrono::duration<double, milli>(diff4).count()
             << " " << chrono::duration<double, milli>(diff5).count() << " ms" << endl;

        return chi2Statistics;
    }




}