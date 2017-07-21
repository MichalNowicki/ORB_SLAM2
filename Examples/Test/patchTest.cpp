//
// Created by michalnowicki on 03.07.17.
//

#include <iostream>
#include "Modifications/PatchRefinement.h"
#include<opencv2/core/core.hpp>
#include <opencv2/opencv.hpp>


void testOptimization(cv::Mat img, int patchSize, cv::Point2f refP, cv::Point2f curP, double minimalStep) {
    patchSize = 5;

    std::cout <<"Ref pos " << refP.x << " " << refP.y << std::endl;
    std::cout <<"Cur pos " << curP.x << " " << curP.y << std::endl;
    std::cout <<"Patch size " << patchSize << std::endl;


    Eigen::Matrix3d H = Eigen::Matrix3d::Identity();
    cv::Point2f result;

    int largePatchSize = patchSize * 2;
    cv::Mat refPatch = img(cv::Rect(refP.y - patchSize,refP.x - patchSize, largePatchSize, largePatchSize));

    cv::Mat curPatch = img(cv::Rect(curP.y - patchSize,curP.x - patchSize, largePatchSize, largePatchSize));

    std::cout <<"Large patch sizes " << refPatch.size() << std::endl;

    PatchRefinement patchRefinement(patchSize);
    bool success = patchRefinement.optimizePosition(refPatch, refP, refP, 1.0, curPatch, curP, 1.0, H, result);
}

void testHomography(cv::Mat image2, cv::Mat point3DInImg1, cv::Mat K1, cv::Mat K2, cv::Mat R, cv::Mat t, int patchSize) {

    std::cout << "testHommography()" << std::endl;
    std::cout << "-> point3DInImg1 = " << std::endl << point3DInImg1 << std::endl;
    std::cout << "-> K1 = " << std::endl << K1 << std::endl;
    std::cout << "-> K2 = " << std::endl << K2 << std::endl;
    std::cout << "-> R = " << std::endl << R << std::endl;
    std::cout << "-> t = " << std::endl << t << std::endl;
    std::cout << "----------------------"<< std::endl;

    // We project point into image 1
    cv::Mat projPInImg1 = K1 * point3DInImg1;
    projPInImg1 = PatchRefinement::normalize2D(projPInImg1);
    std::cout << "-> projPInImg1 = " << std::endl << projPInImg1 << std::endl;

    // We find the 3D position of point in image 2
    cv::Mat point3DInImg2 = R.inv() * (point3DInImg1 - t);
    std::cout << "-> point3DInImg2 = " << std::endl << point3DInImg2 << std::endl;

    // We project point into image 2
    cv::Mat projPInImg2 = K2 * point3DInImg2;
    projPInImg2 = PatchRefinement::normalize2D(projPInImg2);
    std::cout << "-> projPInImg2 = " << std::endl << projPInImg2 << std::endl;

    //// We prepare everything for compute homography
    // T1w -> coordinate 1 in world is indentity
    cv::Mat T1w = cv::Mat::eye(4, 4, CV_32F);

    // T2w -> coirdinate 2 in world (really 1) is our trasformation
    cv::Mat T2w = cv::Mat::eye(4, 4, CV_32F);
    R.copyTo(T2w.rowRange(0,3).colRange(0,3));
    t.copyTo(T2w.rowRange(0,3).col(3));
    std::cout << "-> T2w = " << std::endl << T2w << std::endl;

    // We assume planar normal
    float dataNormal[3] = { 0, 0, 1 };
    cv::Mat n = cv::Mat(3, 1, CV_32F, dataNormal);
    n = n / cv::norm(n);
    std::cout << "-> n = " << std::endl << n << std::endl;

    // We have to find how much the plane is moved from the origin of the coordinate system
    double d = PatchRefinement::getDistanceToPlane(point3DInImg1, n);
    std::cout << "-> d = " << d << std::endl;

    // We compute the homography
    // Pose of c.s. A in B
    cv::Mat Tba = T2w.inv() * T1w;
    cv::Mat H = PatchRefinement::computeHomography(Tba, n, d, K1, K2);
    std::cout << "-> H = " << std::endl << H << std::endl;

    //// Let's verify the homography
    std::cout << "----------------------"<< std::endl;
    std::cout << "Verifying homography" << std::endl;

    double diff = cv::norm(projPInImg2 - PatchRefinement::normalize2D(H*projPInImg1));
    double diff2 = cv::norm(projPInImg1 - PatchRefinement::normalize2D(H.inv()*projPInImg2));
    std :: cout << "Img1: " << projPInImg1.t() << " " << PatchRefinement::normalize2D(H.inv()*projPInImg2) << std::endl;
    std :: cout << "Img2: " << projPInImg2.t() << " " << PatchRefinement::normalize2D(H*projPInImg1) << std::endl;

    if ( diff < 0.0001 && diff2 < 0.0001)
        std::cout << "\t Homography is OK!" << std::endl;
    else
        std::cout << "\t Homography is WROOOONG!" << std::endl;

    // This is how probably image 1 looks like! (Based on H) - this is only for TEST
    cv::Mat image1;
    cv::warpPerspective(image2, image1, H, image2.size(), cv::INTER_LINEAR, cv::BORDER_CONSTANT, cv::Scalar());

//    PatchRefinement patchRefinement(patchSize);
//
//    cv::Point2f refKpInA, refKpInB, result;
//    refKpInA.x = projPInImg1.at<float>(0);
//    refKpInA.y = projPInImg1.at<float>(1);
//    refKpInB.x = projPInImg2.at<float>(0);
//    refKpInB.y = projPInImg2.at<float>(1);
//    Eigen::Matrix3d Heigen;
//    cv::cv2eigen(H, Heigen);


}

int main(int argc, char **argv)
{
    // We test on real image
    cv::Mat imagecolor = cv::imread("lena.bmp"), image;
    cv::cvtColor(imagecolor,image,CV_BGR2GRAY);

    //
    int patchSize = 11;

    // Verify the optimization
    cv::Point2f refP = cv::Point2f(100,100);
    cv::Point2f curP = cv::Point2f(101,101);
    testOptimization(image, patchSize, refP, curP, 1e-6);

    // Let's check the homography with data from real-life system
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

//    testHomography(image, point3DInImg1, K1, K2, R, t, patchSize);
}