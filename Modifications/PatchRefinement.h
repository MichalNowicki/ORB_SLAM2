//
// Created by michalnowicki on 03.07.17.
//

#ifndef ORB_SLAM2_PATCHREFINEMENT_H
#define ORB_SLAM2_PATCHREFINEMENT_H

#include<opencv2/core/core.hpp>
#include <Eigen/Eigen>
#include <iostream>
#include <vector>

class PatchRefinement {

public:
    PatchRefinement(cv::Mat img) {
         int size = img.rows;

        // Each point saves larger patch but for refinement below value is used. The rest of the patch is used for warping.
        static int patchSize = 5;

        // Steps:
        // 1. Compute image gradient for mPatch
        // 2. Warp refPatch by R,t
        // 3. Prepare
        std::vector<Eigen::Vector2f> gradient(size*size, Eigen::Vector2f::Zero());

        // compute gradient
        Eigen::Matrix2f Hessian = Eigen::Matrix2f::Zero(), HessianInv;
        for (int i = 1; i<size-1;i++)
        {
            for (int j = 1;j<size-1;j++) {
                Eigen::Vector2f Jacobian;
                Jacobian[0] = 0.5 * (img.at<uchar>(i,j+1) - img.at<uchar>(i,j-1));
                Jacobian[1] = 0.5 * (img.at<uchar>(i+1,j) - img.at<uchar>(i-1,j));
                gradient[j+i*size] = Jacobian;

                Hessian += Jacobian * Jacobian.transpose();
            }
        }
        HessianInv = Hessian.inverse();


//        std::cout << "Gradient X:" << std::endl;
//        for (int i=0;i<gradient.size();i++) {
//            std::cout << gradient[i][0] << " " ;
//            if (i%size == size - 1)
//                std::cout << std::endl;
//        }
//        std::cout << "Gradient Y:" << std::endl;
//        for (int i=0;i<gradient.size();i++) {
//            std::cout << gradient[i][1] << " ";
//            if (i % size == size - 1)
//                std::cout << std::endl;
//        }


        // Create some patches to test
        double origX = 100.0, origY = 45.0;
        std::vector<Eigen::Vector2f> gd(patchSize*patchSize, Eigen::Vector2f::Zero());
        std::vector<double> patch = computePatch(img, origX, origY, size, patchSize, gradient, gd);

        std::cout <<"Original patch for " << origX << " " << origY << std::endl;
        printPatch(patch, patchSize);


        // Lets optimize position
        float optX = 100.0, optY = 46;

        for (int i = 0; i < 1000; i++) {
            std::cout <<"----------- Simulated patch for " << optX << " " << optY << std::endl;
            std::vector<double> simulatedPatch = computePatch(img, optX, optY, size, patchSize, gradient, gd);
            printPatch(simulatedPatch, patchSize);


            Eigen::Vector2f res = patchDiff(patch, simulatedPatch, gd);
            std::cout << "Error: " << std::endl << res.transpose() << std::endl;

            Eigen::Vector2f step = HessianInv * res;
            std::cout << "Proposed step: "<< std::endl << step.transpose() << std::endl;

            optX = optX + step[0];
            optY = optY + step[1];
        }

    };

    void printPatch(std::vector<double> patch, int patchSize) {
        for (int i=0;i<patch.size();i++)
        {
            std::cout << patch[i] <<" " ;
            if (i%patchSize == patchSize - 1)
                std::cout<< std::endl;
        }
    }




    Eigen::Vector2f patchDiff(std::vector<double> patch, std::vector<double> optimizedPatch, std::vector<Eigen::Vector2f> imageGradient) {
        float averagePatch = std::accumulate( patch.begin(), patch.end(), 0.0)/patch.size();
        float averageOptPatch = std::accumulate( optimizedPatch.begin(), optimizedPatch.end(), 0.0)/optimizedPatch.size();


        Eigen::Vector2f res = Eigen::Vector2f::Zero();
        for (int i=0;i<patch.size();i++) {
                float d = (optimizedPatch[i]-averageOptPatch) - (patch[i] - averagePatch);

                res += imageGradient[i] * (-d);
        }
        return res;
    }

    //
    std::vector<double> computePatch(cv::Mat img, double x, double y, int size, int patchSize, std::vector<Eigen::Vector2f> gradient, std::vector<Eigen::Vector2f> &gd) {

        const double xInt = int(x), yInt = int(y);
        const double xSub = x - xInt, ySub = y - yInt;

        // From wiki: http://upload.wikimedia.org/math/9/b/4/9b4e1064436ecccd069ea238b656c063.png
        const double topLeft = (1.0 - xSub) * (1.0 - ySub);
        const double topRight = xSub * (1.0 - ySub);
        const double bottomLeft = (1.0 - xSub) * ySub;
        const double bottomRight = xSub * ySub;

        std::vector<double> patch (patchSize*patchSize, 0.0);

        int halfPatchSize = (patchSize - 1)/2;
        int index = 0;
        for (int j = yInt - halfPatchSize; j < yInt + halfPatchSize + 1; j++) {
            for (int i = xInt - halfPatchSize; i < xInt + halfPatchSize + 1; i++) {

                float value = topLeft * img.at<uchar>(j, i) + topRight * img.at<uchar>(j, i + 1) +
                              bottomLeft * img.at<uchar>(j + 1, i) + bottomRight * img.at<uchar>(j + 1, i + 1);
                patch[index] = value;

                gd[index] = topLeft * gradient[j*size+i] + topRight * gradient[j*size+i+1] +
                                    bottomLeft * gradient[(j+1)*size+i] + bottomRight * gradient[(j+1)*size+i+1];

                index++;
            }
        }

        return patch;
    }





    // TODO: Not sure about the coordinate systems! Something is wrong!
//        cv::Mat Tpc = currentKF->GetPoseInverse() * mWorldPos;
//        float depth = Tpc.at<float>(2);
//        ComputeHomography(mpRefKF, currentKF, mNormalVector, depth);

//    cv::Mat ComputeHomography(KeyFrame *&pKF1, KeyFrame *&pKF2, cv::Mat n, double d)
//    {
//        cv::Mat R1w = pKF1->GetRotation();
//        cv::Mat t1w = pKF1->GetTranslation();
//        cv::Mat R2w = pKF2->GetRotation();
//        cv::Mat t2w = pKF2->GetTranslation();
//
//        cv::Mat R12 = R1w*R2w.t();
//        cv::Mat t12 = -R1w*R2w.t()*t2w+t1w;
//
//        const cv::Mat &K1 = pKF1->mK;
//        const cv::Mat &K2 = pKF2->mK;
//
//        // K1'* (R - tn' / d)*K^-1
//        return K1.t()* (R12 -  t12 * n / d )*K2.inv();
//    }

};


#endif //ORB_SLAM2_PATCHREFINEMENT_H
