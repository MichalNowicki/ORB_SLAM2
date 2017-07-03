//
// Created by michalnowicki on 03.07.17.
//

#include <iostream>
#include "Modifications/PatchRefinement.h"
#include<opencv2/core/core.hpp>
#include <opencv2/opencv.hpp>

int main(int argc, char **argv)
{

    cv::Mat imagecolor = cv::imread("lena.bmp"), image;
    cv::cvtColor(imagecolor,image,CV_BGR2GRAY);


//    uchar val[64] = {    1, 0, 0, 1, 1,
//                        1, 0, 0, 0, 1,
//                        1,1 ,1 , 1, 1,
//                         0 ,0 ,0 ,0, 0};
//    cv::Mat image = cv::Mat(8,8, CV_8UC1, val);

    PatchRefinement patchRefinement(image);

}