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
//    float data[3] = { 0, 0, 10};
//    cv::Mat point3DInImg1 = cv::Mat(3, 1, CV_32F, data);
//
//    float dataK[9] = { 500, 0, 320, 0, 500, 240, 0, 0, 1};
//    cv::Mat K1 = cv::Mat(3, 3, CV_32F, dataK);
//    cv::Mat K2 = K1;
//
////        float dataRotation[9] = {   0.866, 0.5, 0.0,
////                                    -0.5, 0.866, 0,
////                                    0.0, 0.0, 1.0};
//    float dataRotation[9] = {   1.0, 0.0, 0.0,
//                                0,  1.0, 0,
//                                0.0, 0.0, 1.0};
//    cv::Mat R = cv::Mat(3, 3, CV_32F, dataRotation);
//    float dataTranslation[3] = { 3.2, 0, 0 };
//    cv::Mat t = cv::Mat(3, 1, CV_32F, dataTranslation);


    float data[3] = {-0.465698,
    -0.075128041,
    0.53546041};
    cv::Mat point3DInImg1 = cv::Mat(3, 1, CV_32F, data);

    float dataK[9] = { 256, 0, 319.5, 0, 254.4, 239.5, 0, 0, 1};
    cv::Mat K1 = cv::Mat(3, 3, CV_32F, dataK);
    cv::Mat K2 = K1;

    float dataRotation[9] = {   0.99626195, 0.075956553, -0.041143145,
                                -0.068531439, 0.98492026, 0.1588567,
                                0.052588914, -0.15544328, 0.986444};
    cv::Mat R = cv::Mat(3, 3, CV_32F, dataRotation);
    float dataTranslation[3] = { -0.14661238, 0.098521031, 0.03928398};
    cv::Mat t = cv::Mat(3, 1, CV_32F, dataTranslation);

//    minAngle: 13.1493	 Dist: 46.6132 Reproj1: 0.303948 Reproj2: 0.505587
//    Kp1: 97.0445 204.042 Projected: 96.8529 203.806
//    Kp2 168.41 93.1627 Projected: 168.647 93.6091
//    World: [-1.4450988;
//    -0.31460476;
//    0.22954626]
//    H:
//    [1.0900733, 0.078431956, -143.68582;
//    -0.066213891, 0.93235451, 20.601089;
//    0.00032615208, -0.0002279077, 0.94161034]
//    Kp2in2 48.395 98.7493
//    Tab:
//    [0.99626195, 0.075956553, -0.041143145, -0.14661238;
//    -0.068531439, 0.98492026, 0.1588567, 0.098521031;
//    0.052588914, -0.15544328, 0.986444, 0.03928398;
//    0, 0, 0, 1]
//    point in A:
//    [-0.465698;
//    -0.075128041;
//    0.53546041]
//    Ka:
//    [256, 0, 319.5;
//    0, 254.39999, 239.5;
//    0, 0, 1] [256, 0, 319.5;
//    0, 254.39999, 239.5;
//    0, 0, 1]





    int patchSize = 5;
    patchRefinement.testHomography(image, point3DInImg1, K1, K2, R, t, patchSize);

}