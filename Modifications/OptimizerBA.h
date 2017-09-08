//
// Created by michalnowicki on 08.09.17.
//

#ifndef ORB_SLAM2_OPTIMIZER_H
#define ORB_SLAM2_OPTIMIZER_H


#include "Map.h"
#include "MapPoint.h"
#include "KeyFrame.h"
#include "LoopClosing.h"
#include "Frame.h"

#include "Thirdparty/g2o/g2o/types/types_seven_dof_expmap.h"
#include "Thirdparty/g2o/g2o/types/EdgeProjectPSI2UV.h"
#include "Thirdparty/g2o/g2o/types/EdgeProjectPSI2UVSingleParam.h"
#include "Thirdparty/g2o/g2o/types/EdgeProjectPSI2UVPatch.h"
#include "Thirdparty/g2o/g2o/types/EdgeProjectPSI2UVSingleParamPatch.h"
#include "Thirdparty/g2o/g2o/types/EdgeProjectPSI2UVSingleParamPatchBright.h"

namespace ORB_SLAM2
{

    class LoopClosing;

    class OptimizerBA
    {
    public:
        enum TYPE {
            INVERSE_DEPTH,
            INVERSE_DEPTH_SINGLE_PARAM,
            INVERSE_DEPTH_SINGLE_PARAM_PATCH,
            INVERSE_DEPTH_SINGLE_PARAM_PATCH_BRIGHTNESS,
        };


        std::vector<double> static LocalBundleAdjustment(KeyFrame *pKF, bool *pbStopFlag, Map *pMap, OptimizerBA::TYPE type, float sigma = 1.0);

        std::vector<double> static GlobalBundleAdjustment(Map* pMap, OptimizerBA::TYPE type, float sigma = 1.0);


        // Possible Bundle Adjustments
        std::vector<double> static BundleAdjustmentInvDepth(list<KeyFrame *> lLocalKeyFrames, list<KeyFrame *> lFixedCameras,
                                                            set<MapPoint *> lLocalMapPoints, bool *pbStopFlag, Map *pMap,
                                                            float sigma = 1.0);
        std::vector<double> static BundleAdjustmentInvDepthSingleParam(list<KeyFrame *> lLocalKeyFrames, list<KeyFrame *> lFixedCameras,
                                                                       set<MapPoint *> lLocalMapPoints, bool *pbStopFlag, Map *pMap,
                                                                       float sigma = 1.0);
        std::vector<double> static BundleAdjustmentInvDepthSingleParamPatch(list<KeyFrame *> lLocalKeyFrames, list<KeyFrame *> lFixedCameras,
                                                                            set<MapPoint *> lLocalMapPoints, bool *pbStopFlag, Map *pMap,
                                                                            float sigma = 1.0);
        std::vector<double> static BundleAdjustmentInvDepthSingleParamPatchBright(list<KeyFrame *> lLocalKeyFrames, list<KeyFrame *> lFixedCameras,
                                                                                  set<MapPoint *> lLocalMapPoints, bool *pbStopFlag, Map *pMap,
                                                                                  float sigma = 1.0);

        void static saveBAProblem(KeyFrame *pKF, bool* pbStopFlag, Map* pMap, float sigma, std::string name);


        void static prepareKFsAndMapPoints(KeyFrame *pKF, Map *pMap, list<KeyFrame *> &lLocalKeyFrames, set<MapPoint *> &lLocalMapPoints, list<KeyFrame *> &lFixedCameras);

        void static VerifyReprojectionError(list<KeyFrame *> lLocalKeyFrames, set<MapPoint *> lLocalMapPoints, list<KeyFrame *> lFixedCameras);

    };

} //namespace ORB_SLAM

#endif //ORB_SLAM2_OPTIMIZER_H
