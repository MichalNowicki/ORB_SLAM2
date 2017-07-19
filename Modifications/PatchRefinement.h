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
    PatchRefinement(int _patchSize, double _minIterStep) {
        patchSize = _patchSize;
        minIterStep = _minIterStep;

        iterationNumber = 15;
        stepStopThreshold = 0.10;
    };


    /*
     * Optimizes the position of the feature
     */
    bool optimizePosition(cv::Mat refLargePatch, cv::Point2f refPoint, cv::Point2f refPointCorrection, float refScale, cv::Mat curLargePatch, cv::Point2f curP,
                                     float curScale, Eigen::Matrix3d H, cv::Point2f &correction);

    /*
     * TODO: Remove it
     */
    void performOptimizationForTestS(cv::Mat refImg, cv::Point2f refP, float refScale, cv::Mat curImg, cv::Point2f curP,
                                     float curScale, Eigen::Matrix3d H);


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
     *
     */
    void printGradient(std::vector<Eigen::Vector2d> gradient);

    /*
     * Methods that computes the homography based on:
     * - Tba -> Pose of a coordinate system A in the B's coordinate system (Moves point in A to B)
     * - n   -> Normal of the plane in a coordinate system A. It should be normalized!
     * - d   -> The closest distance to the plane in a coordinate system A. If you only know point use getDistanceToPlane
     * - Ka  -> Camera matrix for a coordinate system a
     * - Kb  -> Camera matrix for a coordinate system b
     */
    cv::Mat ComputeHomography(cv::Mat Tba, cv::Mat n, double d, cv::Mat Ka, cv::Mat Kb);

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
    std::vector<double> computePatch(cv::Mat img, cv::Point2f kp, Eigen::Matrix3d H);
    std::vector<double> computePatch(cv::Mat img, Eigen::Vector2d kp, Eigen::Matrix3d H);


    /*
     * Computes the patch around point (x,y) in img while also computing gradient in those subpix positions based on precomputed gradient
     *
     * size is number of elements in a row of a gradient vector (not nice)
     */
    std::vector<double>
    computePatch(cv::Mat img, Eigen::Vector2d kp, std::vector<Eigen::Vector2d> gradient,
                 std::vector<Eigen::Vector2d> &gd);


    /*
     * TODO
     */
    std::vector<double> computePatchOnSubImage(cv::Mat img, Eigen::Matrix3d H, cv::Point2f img2kp, double scaleKp2, cv::Point2f img1kp, double scaleKp1, cv::Point2f img1kpCorrection);

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
    void computeImageGradient(cv::Mat &img, cv::Point2f kp, std::vector<Eigen::Vector2d> &gradient, Eigen::Matrix2d &HessianInv);

private:
    int patchSize;
    double minIterStep;

    int iterationNumber;
    double stepStopThreshold;

};





#endif //ORB_SLAM2_PATCHREFINEMENT_H
