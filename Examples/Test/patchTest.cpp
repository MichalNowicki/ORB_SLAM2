//
// Created by michalnowicki on 03.07.17.
//

#include <iostream>
#include "Modifications/PatchRefinement.h"
#include<opencv2/core/core.hpp>
#include <opencv2/opencv.hpp>

int main(int argc, char **argv)
{

    //    uchar val[64] = {    1, 0, 0, 1, 1,
//                        1, 0, 0, 0, 1,
//                        1,1 ,1 , 1, 1,
//                         0 ,0 ,0 ,0, 0};
//    cv::Mat image = cv::Mat(8,8, CV_8UC1, val);


    // We test on real image
    cv::Mat imagecolor = cv::imread("lena.bmp"), image;
    cv::cvtColor(imagecolor,image,CV_BGR2GRAY);

    // Class to verify
    PatchRefinement patchRefinement;

    // Verify the optimization
  //  patchRefinement.testOptimization(image, 5, Eigen::Vector2d(50,50), Eigen::Vector2d(49.5,49.5));

    // Let's check the homography
    float data[3] = { 0, 0, 10};
    cv::Mat point3DInImg1 = cv::Mat(3, 1, CV_32F, data);

    float dataK[9] = { 500, 0, 320, 0, 500, 240, 0, 0, 1};
    cv::Mat K1 = cv::Mat(3, 3, CV_32F, dataK);
    cv::Mat K2 = K1;

//        float dataRotation[9] = {   0.866, 0.5, 0.0,
//                                    -0.5, 0.866, 0,
//                                    0.0, 0.0, 1.0};
    float dataRotation[9] = {   1.0, 0.0, 0.0,
                                0,  1.0, 0,
                                0.0, 0.0, 1.0};
    cv::Mat R = cv::Mat(3, 3, CV_32F, dataRotation);
    float dataTranslation[3] = { 3.2, 0, 0 };
    cv::Mat t = cv::Mat(3, 1, CV_32F, dataTranslation);

    int patchSize = 5;
    patchRefinement.testHomography(image, point3DInImg1, K1, K2, R, t, patchSize);

}