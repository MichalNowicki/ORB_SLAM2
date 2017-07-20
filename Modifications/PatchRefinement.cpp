//
// Created by michalnowicki on 03.07.17.
//

#include "PatchRefinement.h"
#include <limits>

Eigen::Vector2d PatchRefinement::computePatchDifference(std::vector<double> patch, std::vector<double> optimizedPatch,
                                       std::vector<Eigen::Vector2d> imageGradient) {

    float averagePatch = std::accumulate( patch.begin(), patch.end(), 0.0)/patch.size();
    float averageOptPatch = std::accumulate( optimizedPatch.begin(), optimizedPatch.end(), 0.0)/optimizedPatch.size();

// TODO: Turned off mean for tests
//    float averagePatch = 0.0;
//    float averageOptPatch = 0.0;

    Eigen::Vector2d res = Eigen::Vector2d::Zero();
    for (int i=0;i<patch.size();i++) {
        float d = (optimizedPatch[i]-averageOptPatch) - (patch[i] - averagePatch);
        res += imageGradient[i] * (-d);
    }
    return res;
}

bool PatchRefinement::optimizePosition(cv::Mat refLargePatch, cv::Point2f refPatchLoc, cv::Point2f refPoint, float refScale,
                                       cv::Mat curLargePatch, cv::Point2f curPatchLoc, cv::Point2f curPoint, float curScale,
                                       Eigen::Matrix3d H, cv::Point2f &subpixCorrection) {


    // The index of the centre of the larger patch
    double centerPosLargePatch = (refLargePatch.rows - 1) / 2;

    // Computation of the ref patch
    std::vector<double> refPatch = computePatchOnSubImage(refLargePatch, refPatchLoc, refScale, curPoint, curScale, H);

    // The warp was outside the saved neighbourhood
    if (refPatch.size()  == 0) {
        if (verbose)
            std::cout << "\tRefPatch is empty -> warp was outside saved subImg" << std::endl;
        return false;
    }

    // Compute gradient & hessian on the current image
    std::vector<Eigen::Vector2d> gradient;
    Eigen::Matrix2d HessianInv;
    cv::Point2f center(centerPosLargePatch, centerPosLargePatch);
    computeImageGradient(curLargePatch, center, gradient, HessianInv);

    // Let's find the subpix gradient
    std::vector<Eigen::Vector2d> subPosGradient(patchSize * patchSize, Eigen::Vector2d::Zero());

    // Lets start the optimization
    //  It is started in (centerPosLargePatch, centerPosLargePatch) of the saved images neighbourhood
    Eigen::Vector2d curOpt;
    curOpt[0] = center.x;
    curOpt[1] = center.y;

    double prevErr = std::numeric_limits<double>::max();
    Eigen::Vector2d prevStep = Eigen::Vector2d::Zero();

    // We perform iterations or until stop condition is met
    for (int i = 0; i < iterationNumber; i++) {
        if(verbose)
            std::cout <<"Iteration: " << i << " Position: " << curOpt[0] << ", " << curOpt[1] << std::endl;

        // The coordinates might be nan when hessian is not invertible
        if(std::isnan(curOpt[0]) || std::isnan(curOpt[1])) {
            if (verbose)
                std::cout << "\tNaN positions. Probably Hessian is not invertible" << std::endl;
            return false;
        }

        // Compute new patch
        std::vector<double> currentPatch = computePatch(curLargePatch, curOpt, gradient, subPosGradient);


        // Residuals
        double err = computePatchDifference(refPatch, currentPatch);
        Eigen::Vector2d res = computePatchDifference(refPatch, currentPatch, subPosGradient);

        // Obtained error
//        if(verbose)
//        std::cout << "Error: " << err << "   \tavg pixel err: " << err / (patchSize * patchSize) << std::endl;

//        if ( err / (patchSize * patchSize)  > 50 )
//        {
//            std::cout << "refPatch vs currentPatch" << std::endl;
//            for (int i =0;i<refPatch.size(); i++) {
//                std::cout << refPatch[i] << " vs " << currentPatch[i] << std::endl;
//            }
//
//        }

        // We stop if the error is larger than previously obtained
        if (err > prevErr) {
            curOpt = curOpt - prevStep;
            break;
        }

        // Step based on Hessian
        Eigen::Vector2d step = HessianInv * res;

        // Obtained step
        if(verbose)
            std::cout << "Proposed step: " << std::endl << step.transpose() << std::endl;

        // Stop condition
//        if (step.squaredNorm() < stepStopThreshold)
//            break;

        // Do the step
        curOpt = curOpt + step;

        // There is a limit how much we can and should move less than 1/4 of the saved patch
        if(curOpt[0] > (centerPosLargePatch + halfPatchSize / 2) || curOpt[0] < (centerPosLargePatch - halfPatchSize / 2) || curOpt[1] > (centerPosLargePatch + halfPatchSize / 2) || curOpt[1] < (centerPosLargePatch - halfPatchSize / 2)) {
            if (verbose)
                std::cout << "\tWe moved too much! optX = " << curOpt[0] << " optY = " << curOpt[1] << std::endl;
            return false;
        }

        prevErr = err;
        prevStep = step;
    }
    double dx = (curOpt[0] - centerPosLargePatch) * curScale;
    double dy = (curOpt[1] - centerPosLargePatch) * curScale;

    if(verbose) {
        std::cout << "Base image changed by: " << (curOpt[0] - centerPosLargePatch) * curScale << ", "
                  << (curOpt[1] - centerPosLargePatch) * curScale << " scale: " << curScale << std::endl;
    }

    std::vector<double> currentPatch = computePatch(curLargePatch, curOpt, gradient, subPosGradient);
    double errorBefore = computePatchDifference(refPatch, currentPatch);

    curOpt[0] = center.x;
    curOpt[1] = center.y;
    currentPatch = computePatch(curLargePatch, curOpt, gradient, subPosGradient);
    double errorAfter = computePatchDifference(refPatch, currentPatch);
//    std::cout << "CORRECTION: " << dx << ", " << dy << " Curr err: " << errorBefore << " Old err: " << errorAfter << std::endl;


    subpixCorrection = cv::Point2f(dx, dy);
    return true;
}


double PatchRefinement::computePatchDifference(std::vector<double> patch, std::vector<double> optimizedPatch) {
    float averagePatch = std::accumulate( patch.begin(), patch.end(), 0.0)/patch.size();
    float averageOptPatch = std::accumulate( optimizedPatch.begin(), optimizedPatch.end(), 0.0)/optimizedPatch.size();

    double res = 0;
    for (int i=0;i<patch.size();i++) {
        res += pow((optimizedPatch[i]-averageOptPatch) - (patch[i] - averagePatch),2);
    }
    return std::sqrt(res);
}

void PatchRefinement::printPatch(std::vector<double> patch, int patchSize) {
    for (int i=0;i<patch.size();i++)
    {
        std::cout << patch[i] <<" " ;
        if (i%patchSize == patchSize - 1)
            std::cout<< std::endl;
    }
}


void PatchRefinement::showPatch(std::vector<double> patch, int patchSize, std::string windowName) {
    cv::Mat patchToShow = cv::Mat(patchSize, patchSize, CV_8U);
    for (int i=0;i<patch.size();i++)
        patchToShow.at<uchar>(i/patchSize,i%patchSize) = uchar(patch[i]);
    imshow(windowName, patchToShow);
}

cv::Mat PatchRefinement::computeHomography(cv::Mat Tba, cv::Mat n, double d, cv::Mat Ka, cv::Mat Kb)
{
    // Getting R,t
    cv::Mat R21 = extractRotation(Tba);
    cv::Mat t21 = extractTranslation(Tba);

    // Homography
    return Ka * (R21 -  t21 * n.t() / d )*Kb.inv();
}

cv::Mat PatchRefinement::extractRotation(cv::Mat T)
{
    return T.rowRange(0,3).colRange(0,3).clone();
}

cv::Mat PatchRefinement::extractTranslation(cv::Mat T)
{
    return T.rowRange(0,3).col(3).clone();
}

double PatchRefinement::getDistanceToPlane(const cv::Mat &point3D, const cv::Mat &normal)  {
    cv::Mat nTimesPoint = normal.t() * point3D;
    return -nTimesPoint.at<float>(0,0);
}

cv::Mat PatchRefinement::normalize2D(cv::Mat p) {
    double z = p.at<float>(2,0);
    p.col(0) = p.col(0) / z;
    return p;
}

//std::vector<double> PatchRefinement::computePatch(cv::Mat img, cv::Point2f kp, Eigen::Matrix3d H) {
//    Eigen::Vector2d p;
//    p[0] = kp.x;
//    p[1] = kp.y;
//    return computePatch(img, p, H);
//}
//
//std::vector<double> PatchRefinement::computePatch(cv::Mat img, Eigen::Vector2d kp, Eigen::Matrix3d H) {
//
//    std::vector<double> patch (patchSize*patchSize, 0.0);
//
//    int halfPatchSize = (patchSize - 1)/2;
//    int index = 0;
//    for (int j = kp[1] - halfPatchSize; j < kp[1] + halfPatchSize + 1; j++) {
//        for (int i = kp[0] - halfPatchSize; i < kp[0] + halfPatchSize + 1; i++) {
//
//            // We find the position on the image1 based on (px,py) in image2 and H
//            Eigen::Vector3d pIn2(i,j,1);
//            Eigen::Vector3d pIn1 = H * pIn2;
//            pIn1 = pIn1 / pIn1(2);
//
//            const double xInt = int(pIn1(0)), yInt = int(pIn1(1));
//            const double xSub = pIn1(0) - xInt, ySub = pIn1(1) - yInt;
//
//            // From wiki: http://upload.wikimedia.org/math/9/b/4/9b4e1064436ecccd069ea238b656c063.png
//            const double topLeft = (1.0 - xSub) * (1.0 - ySub);
//            const double topRight = xSub * (1.0 - ySub);
//            const double bottomLeft = (1.0 - xSub) * ySub;
//            const double bottomRight = xSub * ySub;
//
//
//            float value = topLeft * img.at<uchar>(yInt, xInt) + topRight * img.at<uchar>(yInt, xInt + 1) +
//                          bottomLeft * img.at<uchar>(yInt + 1, xInt) + bottomRight * img.at<uchar>(yInt + 1, xInt + 1);
//            patch[index] = value;
//
//            index++;
//        }
//    }
//
//    return patch;
//}


std::vector<double> PatchRefinement::computePatchOnSubImage(cv::Mat refLargePatch, cv::Point2f refPatchLoc,  double refScale,
                                                            cv::Point2f curPoint, double curScale, Eigen::Matrix3d H) {

    // The patch to compute
    std::vector<double> patch (patchSize*patchSize, 0.0);

    // Let's iterate over all points in patch for point in image B and use homography to find them in image A
    int index = 0;
    for (int j = - halfPatchSize; j <  halfPatchSize + 1; j++) {
        for (int i = - halfPatchSize; i < halfPatchSize + 1; i++) {

            // Computation of the points of patch in image2 taking into account the scale
            double x = curPoint.x + i * curScale;
            double y = curPoint.y + j * curScale;

            // Projecting (x,y) from image 2 into image 1 with  H
            Eigen::Vector3d pIn2(x,y,1);
            Eigen::Vector3d pIn1 = H * pIn2;
            pIn1 = pIn1 / pIn1(2);

            // We only have a image part in image 1!
            // Therefore, we find the point w.r.t the patch location (center), scale that distance by octave's scale and just moves it to the center of stored patch
            int centerPos = (refLargePatch.rows - 1) / 2;
            pIn1(0) = (pIn1(0) - refPatchLoc.x)/refScale + centerPos;
            pIn1(1) = (pIn1(1) - refPatchLoc.y)/refScale + centerPos;

            // If warped outside of the stored large patch
            if ( pIn1(0) < 0 || pIn1(0) > refLargePatch.cols - 1 || pIn1(1) < 0 || pIn1(1) > refLargePatch.rows - 1)
                return std::vector<double>();

            // Finding the subpix value using bilinear interpolation
            const double xInt = int(pIn1(0)), yInt = int(pIn1(1));
            const double xSub = pIn1(0) - xInt, ySub = pIn1(1) - yInt;

            // From wiki: http://upload.wikimedia.org/math/9/b/4/9b4e1064436ecccd069ea238b656c063.png
            const double topLeft = (1.0 - xSub) * (1.0 - ySub);
            const double topRight = xSub * (1.0 - ySub);
            const double bottomLeft = (1.0 - xSub) * ySub;
            const double bottomRight = xSub * ySub;

            // Value of that patch
            patch[index] = topLeft * refLargePatch.at<uchar>(yInt, xInt) + topRight * refLargePatch.at<uchar>(yInt, xInt + 1) +
                          bottomLeft * refLargePatch.at<uchar>(yInt + 1, xInt) + bottomRight * refLargePatch.at<uchar>(yInt + 1, xInt + 1);

            index++;
        }
    }

    return patch;
}

std::vector<double> PatchRefinement::computePatch(cv::Mat img, Eigen::Vector2d kp, std::vector<Eigen::Vector2d> gradient, std::vector<Eigen::Vector2d> &gd) {

    const double xInt = int(kp[0]), yInt = int(kp[1]);
    const double xSub = kp[0] - xInt, ySub = kp[1] - yInt;

    // From wiki: http://upload.wikimedia.org/math/9/b/4/9b4e1064436ecccd069ea238b656c063.png
    const double topLeft = (1.0 - xSub) * (1.0 - ySub);
    const double topRight = xSub * (1.0 - ySub);
    const double bottomLeft = (1.0 - xSub) * ySub;
    const double bottomRight = xSub * ySub;

    std::vector<double> patch (patchSize*patchSize, 0.0);

    int halfPatchSize = (patchSize - 1)/2;
    int index = 0;
    int gradientStride = patchSize + 1;

    for (int y = 0; y < patchSize; y++) {
        for (int x = 0; x < patchSize ; x++) {
            int i = xInt - halfPatchSize + x;
            int j = yInt - halfPatchSize + y;

            patch[index] = topLeft * img.at<uchar>(j, i) + topRight * img.at<uchar>(j, i + 1) +
                           bottomLeft * img.at<uchar>(j + 1, i) + bottomRight * img.at<uchar>(j + 1, i + 1);

            gd[index] = topLeft * gradient[index+y] + topRight * gradient[index+1+y] +
                        bottomLeft * gradient[index+gradientStride+y] + bottomRight * gradient[index+gradientStride+1+y];

            index++;
        }
    }

    return patch;
}

void PatchRefinement::computeImageGradient(cv::Mat &img, cv::Point2f kp, std::vector<Eigen::Vector2d> &gradient, Eigen::Matrix2d &HessianInv) {

    // If patchSize is p, then we need gradient in (p+1)*(p+1) rectangle. The gradient is computed with central diff so (p+3)*(p+3)
    int optXInt = floor(kp.x), optYInt = floor(kp.y);
    int minX = optXInt - (patchSize - 1)/2, maxX = optXInt + (patchSize - 1)/2;
    int minY = optYInt - (patchSize - 1)/2, maxY = optYInt + (patchSize - 1)/2;
    cv::Mat patchForGradient = img.rowRange(minY-1, maxY+3).colRange(minX-1, maxX+3);

    int size = patchForGradient.rows;
    gradient = std::vector<Eigen::Vector2d>((size-2)*(size-2), Eigen::Vector2d::Zero());
    Eigen::Matrix2d Hessian = Eigen::Matrix2d::Zero();

    for (int i = 1; i<size-1;i++)
    {
        for (int j = 1;j<size-1;j++) {
            Eigen::Vector2d Jacobian;
            Jacobian[0] = 0.5 * (patchForGradient.at<uchar>(i,j+1) - patchForGradient.at<uchar>(i,j-1));
            Jacobian[1] = 0.5 * (patchForGradient.at<uchar>(i+1,j) - patchForGradient.at<uchar>(i-1,j));
            gradient[(j-1)+(i-1)*(size-2)] = Jacobian;

//            std :: cout << "Jacobian ! " << std::endl << Jacobian << std::endl;
//
//            std::cout << "Hessian : " << std::endl << Hessian << std::endl;
//            std::cout << "J*J^t : " << std::endl << Jacobian * Jacobian.transpose() << std::endl;

            Hessian += Jacobian * Jacobian.transpose();
        }
    }
    HessianInv = Hessian.inverse();
};

void PatchRefinement::printGradient(std::vector<Eigen::Vector2d> gradient) {
    int noInLine = std::sqrt(gradient.size());

    for (int j=0;j<2;j++)
    {
        if (j == 0)
            std::cout<<"X" << std::endl;
        else
            std::cout<<"Y" << std::endl;

        for (int i=0;i<gradient.size();i++) {
            std::cout << gradient[i][j] << " ";
            if (i % noInLine == noInLine - 1)
                std::cout << std::endl;
        }
    }

}


//void PatchRefinement::performOptimizationForTestS(cv::Mat refImg, cv::Point2f refP, float refScale, cv::Mat curImg, cv::Point2f curP,
//                                 float curScale, Eigen::Matrix3d H) {
//
//    // Ref patch
//    std::vector<double> refPatch = computePatch(refImg, refP, Eigen::Matrix3d::Identity());
//
////        std::vector<double> refPatch = computePatchOnSubImage(refImg, H.inverse(), curP, curScale,
////                                                              refP, refScale);
//
//    // Compute gradient in the current image
//    std::vector<Eigen::Vector2d> gradient;
//    Eigen::Matrix2d HessianInv;
//    computeImageGradient(curImg, curP, gradient, HessianInv);
//
//    std::cout << "PATCH SIZE: " << patchSize << " Gradient size: " << gradient.size() << std::endl;
//
//    // Iteration counter
//    Eigen::Vector2d currentP;
//    currentP[0] = curP.x;
//    currentP[1] = curP.y;
//
//    std::vector<Eigen::Vector2d> subPosGradient(patchSize*patchSize, Eigen::Vector2d::Zero());
//    for (int i = 0; i < iterationNumber; i++) {
//        std::cout <<"Iteration: " << i << " Position: " << currentP[0] << ", " << currentP[1] << std::endl;
//
//        // Compute new patch
//        std::vector<double> movedPatch = computePatch(curImg, currentP, gradient, subPosGradient);
//
//        // Estimate residuals
//        Eigen::Vector2d Jres = computePatchDifference(refPatch, movedPatch, subPosGradient);
//
//        // Obtained error
//        std::cout << "Error: " << std::endl << Jres.transpose() << std::endl;
//
//        // Hessian
//        std::cout << "HessianInv: " << std::endl << HessianInv << std::endl;
//
//        // Step based on Hessian
//        Eigen::Vector2d step = HessianInv * Jres;
//
//        // Obtained error
//        std::cout << "Step: " << std::endl << step.transpose() << std::endl;
//
//        // Stop condition
////        if (step.squaredNorm() < minIterStep)
////            break;
//
//        // We do the step
//        currentP = currentP + step;
//    }
//
//    std::cout << "Final position : " << currentP[0] << ", " << currentP[1] << std::endl;
//
//
//}


