//
// Created by michalnowicki on 03.07.17.
//

#ifndef ORB_SLAM2_PATCHREFINEMENT_H
#define ORB_SLAM2_PATCHREFINEMENT_H

#include <Eigen/Eigen>
#include <opencv2/core/core.hpp>
#include <iostream>
#include <vector>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv/cv.hpp>
#include <opencv2/core/eigen.hpp>

class PatchRefinement {

public:
    PatchRefinement() {
    };

    /*
     * Computes the difference between patches and weights the error with the provided gradient vector
     */
    Eigen::Vector2d computePatchDifference(std::vector<double> patch, std::vector<double> optimizedPatch,
                                           std::vector<Eigen::Vector2d> imageGradient);

    /*
    * Computes the difference between patches
    */
    double computePatchDifference(std::vector<double> patch, std::vector<double> optimizedPatch);


    /*
     * Prints the values of the patch to the console
     */
    void printPatch(std::vector<double> patch, int patchSize);

    /*
     * Creates a window with the provided name to visualize the patch
     */
    void showPatch(std::vector<double> patch, int patchSize, std::string windowName);

    /*
     * Methods that computes the homography based on:
     * - Taw -> Pose of a coordinate system A in the World coordinate system
     * - Tbw -> Pose of a coordinate system B in the World coordinate system
     * - n   -> Normal of the plane in a coordinate system A. It should be normalized!
     * - d   -> The closest distance to the plane in a coordinate system A. If you only know point use getDistanceToPlane
     * - Ka  -> Camera matrix for a coordinate system a
     * - Kb  -> Camera matrix for a coordinate system b
     */
    cv::Mat ComputeHomography(cv::Mat Taw, cv::Mat Tbw, cv::Mat n, double d, cv::Mat Ka, cv::Mat Kb);

    /*
     * Returns the rotational part of transformation
     */
    cv::Mat extractRotation(cv::Mat T);

    /*
    * Returns the translational part of transformation
    */
    cv::Mat extractTranslation(cv::Mat T);

    /*
     * Returns the closest distance to the plane when the position of 1 point and normal is known
     */
    double getDistanceToPlane(const cv::Mat &point3D, const cv::Mat &normal);

    /*
     * Takes 3D points and normalizes it to achieve 2D point
     */
    cv::Mat normalize2D(cv::Mat p);


    /*
     * Computes the patch with homography. It takes point (px,py) in image2 (and its patch surroundings),
     * transforms them to image 1 (p1 = H * p2), performs bilinear interpolation and returns patch as single vector
     * (data is stored row by row and the vector length is patchSize*patchSize)
     *
     * The method can be used to compute patch of original image with H = Identity()
     */
    std::vector<double> computePatch(cv::Mat img, double px, double py, int patchSize, Eigen::Matrix3d H);


    /*
     * Computes the patch around point (x,y) in img while also computing gradient in those subpix positions based on precomputed gradient
     *
     * size is number of elements in a row of a gradient vector (not nice)
     */
    std::vector<double>
    computePatch(cv::Mat img, double x, double y, int size, int patchSize, std::vector<Eigen::Vector2d> gradient,
                 std::vector<Eigen::Vector2d> &gd);


    /*
     * TODO
     */
    std::vector<double> computePatchOnSubImage(cv::Mat img, int patchSize, Eigen::Matrix3d H, cv::Point2f img2kp, double scaleKp2, cv::Point2f img1kp, double scaleKp1);
    Eigen::Matrix3d cv2eigen(cv::Mat H) {
        Eigen::Matrix3d Heigen;
        cv::cv2eigen(H, Heigen);
        return Heigen;
    };


    /*
     * It computes the image gradient on image "img" and stores it row by row in variable gradient. The inverse Hessian is also computed
     *
     * Attention! It assumes square image as all patches will be square!
     */
    void computeImageGradient(cv::Mat &img, std::vector<Eigen::Vector2d> &gradient, Eigen::Matrix2d &HessianInv);


    /*
     * Tests optimization on provided image (img) with selected size of a patch (patchSize). The original patch is taken at originalPoint
     * and the optimization starts at testPoint. The optimization is stopped when 100 iterations are performed or the proposed step is smaller than minimalStep
     */
    void testOptimization(cv::Mat img, int patchSize, Eigen::Vector2d originalPoint, Eigen::Vector2d testPoint,
                          double minimalStep = 1e-6);


    /*
     * Tests the homography estimation on image2, one 3D point, two camera matrices, and pose of cs B in A (R,t) given patchSize
     */
    void testHomography(cv::Mat image2, cv::Mat point3DInImg1, cv::Mat K1, cv::Mat K2, cv::Mat R, cv::Mat t, int patchSize);


};





#endif //ORB_SLAM2_PATCHREFINEMENT_H
