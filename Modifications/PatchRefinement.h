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
    PatchRefinement(int _patchSize, int verbose = 0) : verbose(verbose) {
        patchSize = _patchSize;
        halfPatchSize = (patchSize - 1)/2;

        iterationNumber = 15;
        stepStopThreshold = 0.01;
    };


    /*
     * Optimizes the position of the feature:
     *  refLargePatch  -> subImage from image 01
     *  refPoint       -> the point in the image 01
     *  refPointCorrection -> xXXx
     *  refScale       -> the scale of the point
     *
     *  curLargePatch   -> subImage from image 02
     *  curP           -> the point in the image 02
     *  curScale       -> the scale of the point
     *
     *  H -> homography between 01 and 02
     *
     *  correction -> the proposed correction to the point position based on patches
     *
     */
    bool optimizePosition(cv::Mat refLargePatch, cv::Point2f refPatchLoc, cv::Point2f refPoint, float refScale,
                          cv::Mat curLargePatch, cv::Point2f curPatchLoc, float curScale,
                          Eigen::Matrix3d H, cv::Point2f &subpixCorrection);



    /*
     * Methods that computes the homography based on:
     * - Tba -> Pose of a coordinate system A in the B's coordinate system (Moves point in A to B)
     * - n   -> Normal of the plane in a coordinate system A. It should be normalized!
     * - d   -> The closest distance to the plane in a coordinate system A. If you only know point use getDistanceToPlane
     * - Ka  -> Camera matrix for a coordinate system a
     * - Kb  -> Camera matrix for a coordinate system b
     */
    static cv::Mat computeHomography(cv::Mat Tba, cv::Mat n, double d, cv::Mat Ka, cv::Mat Kb);


    /*
     * Mat to eigen conversion
     */
    static Eigen::Matrix3d inline cv2eigen(cv::Mat H) {
        Eigen::Matrix3d Heigen;
        cv::cv2eigen(H, Heigen);
        return Heigen;
    };

    /*
     * Returns the rotational part of transformation
     */
    static cv::Mat extractRotation(cv::Mat T);

    /*
    * Returns the translational part of transformation
    */
    static cv::Mat extractTranslation(cv::Mat T);

    /*
     * Returns the closest distance to the plane when the position of 1 point and normal is known
     */
    static double getDistanceToPlane(const cv::Mat &point3D, const cv::Mat &normal);

    /*
     * Takes 3D points and normalizes it to achieve 2D point
     */
    static cv::Mat normalize2D(cv::Mat p);

    /*
    * Prints the values of the patch to the console
    */
    void printPatch(std::vector<double> patch, int patchSize);

    /*
     * Creates a window with the provided name to visualize the patch
     */
    void showPatch(std::vector<double> patch, int patchSize, std::string windowName);

    /*
     * Print the values of the gradient to the screen
     */
    void printGradient(std::vector<Eigen::Vector2d> gradient);

private:
    /*
     * Computes the patch with homography. It takes point (px,py) in image2 (and its patch surroundings),
     * transforms them to image 1 (p1 = H * p2), performs bilinear interpolation and returns patch as single vector
     * (data is stored row by row and the vector length is patchSize*patchSize)
     *
     * The method can be used to compute patch of original image with H = Identity()
     */
    std::vector<double> computePatchOnSubImage(cv::Mat curLargePatch, cv::Point2f curPatchLoc,  double curScale,
                                               cv::Point2f refPoint, double refScale, Eigen::Matrix3d H);


    /*
     * Computes the difference between patches and weights the error with the provided gradient vector
     */
    Eigen::Vector2d computePatchDifference(std::vector<double> refPatch, std::vector<double> curPatch,
                                           std::vector<Eigen::Vector2d> imageGradient);

    /*
    * Computes the SSD difference between patches
    */
    double computePatchDifference(std::vector<double> refPatch, std::vector<double> curPatch);


    /*
     * Computes the patch around point (x,y) in img while also computing gradient in those subpix positions based on precomputed gradient
     */
    std::vector<double>
    computePatch(cv::Mat img, Eigen::Vector2d kp, std::vector<Eigen::Vector2d> gradient,
                 std::vector<Eigen::Vector2d> &gd);



    /*
     * It computes the image gradient on image "img" and stores it row by row in variable gradient. The inverse Hessian is also computed
     *
     * Attention! It assumes square image as all patches will be square!
     */
    void computeImageGradient(cv::Mat &img, cv::Point2f kp, std::vector<Eigen::Vector2d> &gradient, Eigen::Matrix2d &HessianInv);


    int verbose;
    int patchSize, halfPatchSize;

    int iterationNumber;
    double stepStopThreshold;

};





#endif //ORB_SLAM2_PATCHREFINEMENT_H
