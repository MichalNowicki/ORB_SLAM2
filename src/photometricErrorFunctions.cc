#include <vector>
#include <cmath>
#include <numeric>
#include <iostream>

#include "photometricErrorFunctions.h"


namespace photo {
    // Gets the subpixel value using bilinear interpolation
    double getSubpixImageValue(double u, double v, std::vector<std::vector<float>> &image) {

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

    double getSubpixImageValue(double u, double v, cv::Mat image) {

        const double xInt = int(u), yInt = int(v);
        const double xSub = u - xInt, ySub = v - yInt;

        const double topLeft = (1.0 - xSub) * (1.0 - ySub);
        const double topRight = xSub * (1.0 - ySub);
        const double bottomLeft = (1.0 - xSub) * ySub;
        const double bottomRight = xSub * ySub;


        if (yInt < 0 || xInt < 0 || yInt + 1 >= image.rows || xInt + 1 >= image.cols) {
            return -1;
        }

        return topLeft * image.at<uchar>(yInt,xInt) +
               topRight * image.at<uchar>(yInt,xInt + 1) +
               bottomLeft * image.at<uchar>(yInt + 1,xInt) +
               bottomRight * image.at<uchar>(yInt + 1,xInt + 1);
    }


    // Computes distance from point3D to plane defined by normal
    double getDistanceToPlane(const Eigen::Vector3d &point3D, const Eigen::Vector3d &normal) {
        return -normal.transpose() * point3D;
    }

    // Creates camera matrix from (fx, y, cx, cy)
    Eigen::Matrix3d getCameraMatrix(float fx, float fy, float cx, float cy) {
        Eigen::Matrix3d cameraMatrix = Eigen::Matrix3d::Identity();
        cameraMatrix(0, 0) = fx;
        cameraMatrix(0, 2) = cx;
        cameraMatrix(1, 1) = fy;
        cameraMatrix(1, 2) = cy;
        cameraMatrix(1, 2) = cy;

        return cameraMatrix;
    }

    // Computes homography
    Eigen::Matrix3d
    computeHomography(Eigen::Matrix4d Tba, Eigen::Vector3d n, double d, Eigen::Matrix3d Ka, Eigen::Matrix3d Kb) {
        // Getting R,t
        Eigen::Matrix3d R21 = Tba.block<3, 3>(0, 0);
        Eigen::Vector3d t21 = Tba.block<3, 1>(3, 0);

        Eigen::Matrix3d H = Ka * (R21 - t21 * n.transpose() / d) * Kb.inverse();

        // Homography
        return H;
    }

    // Computes the difference between patches with affine parameters (a,b) esitmation with (a*ref + b = cur)
    double computePatchDiffAffine(const std::vector<double> &refPatch, const std::vector<double> &curPatch) {
        // Least squares solution to linear brightness change: (alfa * curPatch + beta = refPatch)
        Eigen::MatrixXf A = Eigen::MatrixXf::Ones(refPatch.size(), 2);
        Eigen::VectorXf b = Eigen::VectorXf::Zero(curPatch.size());
        for (int i = 0; i < refPatch.size(); i++) {
            A(i, 0) = curPatch[i];
            b(i) = refPatch[i];
        }
        Eigen::VectorXf brightness = A.colPivHouseholderQr().solve(b);

        std::cout << "Affine: " << brightness(0) << " " << brightness(1) << std::endl;

        // Computing the difference
        double rmse = 0;
        for (int i = 0; i < refPatch.size(); i++) {
            float diff = brightness(0) * curPatch[i] + brightness(1) - refPatch[i];
            rmse += diff * diff;
        }
        return sqrt(rmse / refPatch.size());
    }

    // Computes the difference between patches with mean subtraction
    double computePatchDiffAvg(const std::vector<double> &refPatch, const std::vector<double> &curPatch) {
        double averageRefPatch = std::accumulate(refPatch.begin(), refPatch.end(), 0.0) / refPatch.size();
        double averageCurPatch = std::accumulate(curPatch.begin(), curPatch.end(), 0.0) / curPatch.size();

        // Computing the difference
        double rmse = 0;
        for (int i = 0; i < refPatch.size(); i++) {
            float diff = (curPatch[i] - averageCurPatch) - (refPatch[i] - averageRefPatch);
            rmse += diff * diff;
        }
        return sqrt(rmse / refPatch.size());
    }
}