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

    double PhotoTracker::getDistanceToPlane(const cv::Mat &point3D, const cv::Mat &normal) {
        cv::Mat nTimesPoint = normal.t() * point3D;
        return -nTimesPoint.at<float>(0, 0);
    }

    cv::Mat PhotoTracker::normalize2D(cv::Mat p) {
        double z = p.at<float>(2, 0);
        p.col(0) = p.col(0) / z;
        return p;
    }

    cv::Mat PhotoTracker::computeHomography(cv::Mat Tba, cv::Mat n, double d, cv::Mat Ka, cv::Mat Kb) {
        // Getting R,t
        cv::Mat R21 = Tba.rowRange(0, 3).colRange(0, 3);
        cv::Mat t21 = Tba.rowRange(0, 3).col(3);

        cv::Mat H = Ka * (R21 - t21 * n.t() / d) * Kb.inv();

        // Homography
        return H.clone();
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

    cv::Mat PhotoTracker::getCameraMatrix(float fx, float fy, float cx, float cy) {
        float dataK[9] = {fx, 0, cx, 0, fy, cy, 0, 0, 1};
        return cv::Mat(3, 3, CV_32F, dataK).clone();
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

        int nmatches = 0, nrescued =0;

        std::vector<double> errorMatched;

        // For all map points observed in last frame
        for (int i = 0; i < LastFrame.N; i++) {
            MapPoint *pMP = LastFrame.mvpMapPoints[i];



            if (pMP) { //  if (pMP && pMP2) { // TODO: For all not matched or all?
                if (!LastFrame.mvbOutlier[i]) {

                    bool featureMatched = false;
                    for (int j=0;j<CurrentFrame.N;j++) {
                        MapPoint *pMP2 = CurrentFrame.mvpMapPoints[j];
                        if (pMP == pMP2) {
                            featureMatched = true;
                            break;
                        }
                    }

                    // Project it onto current frame to check if depth is positive
                    cv::Mat x3Dw = pMP->GetWorldPos();
                    cv::Mat x3Dc = Rcw * x3Dw + tcw;
                    cv::Mat x3Dl = Rlw * x3Dw + tlw;
                    if (x3Dc.at<float>(3) < 0)
                        continue;


                    // The patch normal in ref observation is assumed to be [0,0,-1]
                    float dataNormal[3] = {0, 0, -1};
                    cv::Mat n = cv::Mat(3, 1, CV_32F, dataNormal);

                    // Computation of the distance to patch plane in image A
                    double d = getDistanceToPlane(x3Dl, n);

                    // Computation of the homography between A and B
                    cv::Mat Twa = getInversePose(LastFrame.mTcw);
                    cv::Mat Tbw = CurrentFrame.mTcw;
                    cv::Mat Tba = Tbw * Twa;
                    cv::Mat Ka = getCameraMatrix(LastFrame.fx, LastFrame.fy, LastFrame.cx, LastFrame.cy);
                    cv::Mat Kb = getCameraMatrix(CurrentFrame.fx, CurrentFrame.fy, CurrentFrame.cx, CurrentFrame.cy);
                    cv::Mat H = computeHomography(Tba, n, d, Ka, Kb);

                    Eigen::Matrix3d Heigen;
                    cv::cv2eigen(H, Heigen);

                    // Location in last frame
                    float lastU = LastFrame.mvKeysUn[i].pt.x;
                    float lastV = LastFrame.mvKeysUn[i].pt.y;

                    // Location in current frame
                    Eigen::Vector3d last(lastU, lastV, 1);
                    Eigen::Vector3d current = Heigen * last;
                    current = current / current(2);
                    float currentU = current(1), currentV = current(2);

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
                        Eigen::Vector3d pInB = Heigen * pInA;
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

                        if (featureMatched)
                            errorMatched.push_back(errorAvg);

                        // If RMSE is lower than threshold then we are succesful
                        if (errorAvg < photoThreshold) {

                            //std::cout << "[" << i << "] -> Matched: " << featureMatched << " | RMSE_Avg = " << errorAvg << std::endl;

                            if (!featureMatched) {
                                // TODO: How to do it?
                                pMP->rescued = true;
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

                                // Simulate descriptor ?

                                nmatches++;
                            }



                        }
                    }
                }
            }
        }

        double sum = 0;
        for (int j=0;j<errorMatched.size();j++)
            sum += errorMatched[j];
        std::cout <<"Average error of matched features: " << sum / errorMatched.size() << std::endl;


        return nmatches;
    }

}