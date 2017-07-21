//
// Created by michalnowicki on 03.07.17.
//

#include "PatchRefinement.h"
#include <limits>

Eigen::Vector2d PatchRefinement::computePatchDifference(std::vector<double> refPatch, std::vector<double> curPatch,
                                       std::vector<Eigen::Vector2d> imageGradient) {

//    float averageRefPatch = std::accumulate( refPatch.begin(), refPatch.end(), 0.0)/refPatch.size();
//    float averageCurPatch = std::accumulate( curPatch.begin(), curPatch.end(), 0.0)/curPatch.size();
//    float averageRefPatch = 0.0;
//    float averageCurPatch = 0.0;

    // Least squares solution to linear brightness change: (alfa * curPatch + beta = refPatch)
    Eigen::MatrixXf A = Eigen::MatrixXf::Ones(refPatch.size(), 2);
    Eigen::VectorXf b = Eigen::VectorXf::Zero(curPatch.size());
    for (int i=0;i<refPatch.size();i++) {
        A(i,0) = curPatch[i];
        b(i) = refPatch[i];
    }
    Eigen::VectorXf brightness = A.colPivHouseholderQr().solve(b);

    Eigen::Vector2d res = Eigen::Vector2d::Zero();
    for (int i=0;i<refPatch.size();i++) {
//        float d = (curPatch[i]-averageCurPatch) - (refPatch[i] - averageRefPatch);
        float d = brightness(0) * curPatch[i] + brightness(1) - refPatch[i];
        res += imageGradient[i] * d;
    }
    return res;
}

double PatchRefinement::computePatchDifference(std::vector<double> refPatch, std::vector<double> curPatch) {
//    float averagePatch = std::accumulate( patch.begin(), patch.end(), 0.0)/patch.size();
//    float averageOptPatch = std::accumulate( optimizedPatch.begin(), optimizedPatch.end(), 0.0)/optimizedPatch.size();

    // TODO: Turned off mean for tests - can do more when affine brightness is also estimated
    float averageRefPatch = 0.0;
    float averageCurPatch = 0.0;

    double res = 0;
    for (int i=0;i<refPatch.size();i++) {
        res += pow((curPatch[i]-averageCurPatch) - (refPatch[i] - averageRefPatch),2);
    }
    return std::sqrt(res);
}


bool PatchRefinement::optimizePosition(cv::Mat refLargePatch, cv::Point2f refPatchLoc, cv::Point2f refPoint, float refScale,
                                       cv::Mat curLargePatch, cv::Point2f curPatchLoc, float curScale,
                                       Eigen::Matrix3d H, cv::Point2f &subpixCorrection) {

    std::cout << " -------- " << std::endl;
    std::cout << " refPatchLoc = " << refPatchLoc.x << " " << refPatchLoc.y << std::endl;
    std::cout << " refPoint = " << refPoint.x << " " << refPoint.y << std::endl;
    std::cout << " curPatchLoc = " << curPatchLoc.x << " " << curPatchLoc.y << std::endl;
    std::cout << " -------- " << std::endl;

    // The helpful variable
    double centerPosLargePatch = (refLargePatch.rows - 1) / 2;

    // The center of the reference patch in large patch c.s.
    cv::Point2f refPosInLargePatch;
    refPosInLargePatch.x = centerPosLargePatch + (refPoint.x - refPatchLoc.x) / refScale;
    refPosInLargePatch.y = centerPosLargePatch + (refPoint.y - refPatchLoc.y) / refScale;

    // Compute gradient & hessian on the current image
    std::vector<Eigen::Vector2d> gradient;
    Eigen::Matrix2d HessianInv;
    computeImageGradient(refLargePatch, refPosInLargePatch, gradient, HessianInv);

    // Let's find the subpix gradient
    std::vector<Eigen::Vector2d> subPosGradient(patchSize * patchSize, Eigen::Vector2d::Zero());

    // Optimization starts at the refPosInLargePatch
    Eigen::Vector2d curOpt;
    curOpt[0] = refPosInLargePatch.x ;
    curOpt[1] = refPosInLargePatch.y ;

    // Compute ref patch
    std::vector<double> refPatch = computePatch(refLargePatch, curOpt, gradient, subPosGradient);

    // Lets start the optimization
    double prevErr = std::numeric_limits<double>::max();
    Eigen::Vector2d prevStep = Eigen::Vector2d::Zero();

    // Iterate until stop condition is met
    for (int i = 0; i < iterationNumber; i++) {
//        if(verbose)
            std::cout <<"***\tIteration: " << i << " Position: " << curOpt[0] << ", " << curOpt[1] << std::endl;

        // The coordinates might be nan when hessian is not invertible
        if(std::isnan(curOpt[0]) || std::isnan(curOpt[1])) {
            if (verbose)
                std::cout << "\tNaN positions. Probably Hessian is not invertible" << std::endl;
            return false;
        }

        // Computation of the cur patch - the original points is moved by the predicted correction
        cv::Point2f centerToWarp;
        centerToWarp.x = refPoint.x + (curOpt[0]-centerPosLargePatch)*curScale;
        centerToWarp.y = refPoint.y + (curOpt[1]-centerPosLargePatch)*curScale;
        std::cout<<"Original position: " << refPoint.x << " " << refPoint.y << " Proposed: " << centerToWarp.x << " " << centerToWarp.y << std::endl;


        std::vector<double> curPatch = computePatchOnSubImage(curLargePatch, curPatchLoc, curScale, centerToWarp, refScale, H);


        // The warp was outside the saved neighbourhood
        if (curPatch.size()  == 0) {
//            if (verbose)
                std::cout << "\tRefPatch is empty -> warp was outside saved subImg" << std::endl;
            return false;
        }


        // Residuals
        double err = computePatchDifference(refPatch, curPatch);
        Eigen::Vector2d res = computePatchDifference(refPatch, curPatch, subPosGradient);

        // Obtained error
//        if(verbose)
        std::cout << "Error: " << err << "   \tavg pixel err: " << err / (patchSize * patchSize) << std::endl;

        // Optimization is stopped if the error is larger and the previous step is retracked
        if (err > prevErr) {
            curOpt = curOpt + prevStep;
            break;
        }

        // Step based on Hessian
        Eigen::Vector2d step = HessianInv * res;

        // Obtained step
//        if(verbose)
            std::cout << "Proposed step: " << std::endl << step.transpose()  << " Step norm: " << step.squaredNorm() << std::endl;

        // Stop condition
        if (step.squaredNorm() < stepStopThreshold*stepStopThreshold)
            break;

        // Do the step
        curOpt = curOpt - step;

        // There is a limit how much we can and should move less than 1/2 of the saved patch
        if(curOpt[0] > (refPosInLargePatch.x + halfPatchSize) || curOpt[0] < (refPosInLargePatch.x - halfPatchSize) ||
                curOpt[1] > (refPosInLargePatch.y + halfPatchSize) || curOpt[1] < (refPosInLargePatch.y - halfPatchSize)) {
//            if (verbose)
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

    cv::Point2f centerToWarp;
    centerToWarp.x = refPoint.x + dx;
    centerToWarp.y = refPoint.y + dy;
    std::vector<double> curPatchAfter = computePatchOnSubImage(curLargePatch, curPatchLoc, curScale, centerToWarp, refScale, H);
    if (curPatchAfter.size()  == 0) {
        std::cout << "\tRefPatch is empty -> warp was outside saved subImg" << std::endl;
        return false;
    }

    double errorAfter = computePatchDifference(refPatch, curPatchAfter);


    std::vector<double> curPatchBefore = computePatchOnSubImage(curLargePatch, curPatchLoc, curScale, refPoint, refScale, H);

    double errorBefore = computePatchDifference(refPatch, curPatchBefore);


    std::cout << "CORRECTION: " << dx << ", " << dy << " Err at start: " << errorBefore << " Err at end: " << errorAfter << std::endl;


//    subpixCorrection = cv::Point2f(dx, dy);

    Eigen::Vector3d corrInB = Eigen::Vector3d::Zero();
    corrInB[0] = -dx;
    corrInB[1] = -dy;

    corrInB = H*corrInB;

    std::cout << "Correction in B: " << corrInB[0] << " " << corrInB[1] << " " << corrInB[2] << std::endl;


    Eigen::Vector3d pInAO(refPatchLoc.x,refPatchLoc.y,1);
    Eigen::Vector3d pInBO = H * pInAO;
    pInBO = pInBO / pInBO(2);

    Eigen::Vector3d pInA(centerToWarp.x,centerToWarp.y,1);
    Eigen::Vector3d pInB = H * pInA;
    pInB = pInB / pInB(2);

    std::cout << "Original A: " << refPatchLoc.x << " " << refPatchLoc.y << std::endl;
    std::cout << "Corrected A: " << centerToWarp.x << " " << centerToWarp.y << std::endl;

    std::cout << "From Match B: " << curPatchLoc.x << " " << curPatchLoc.y << std::endl;
    std::cout << "Homography B: " << pInBO[0] << " " << pInBO[1] << " -> err = " << errorBefore << std::endl;
    std::cout << "Corrected B: " << pInB[0] << " " << pInB[1] << " -> err = " << errorAfter << std::endl;
    std::cout << "Corrected B ver 2: " << pInBO[0]+corrInB[0] << " " << pInBO[1]+corrInB[1] << " -> err = " << errorAfter << std::endl;


    subpixCorrection = cv::Point2f(corrInB[0], corrInB[1]);


    return true;
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



std::vector<double> PatchRefinement::computePatchOnSubImage(cv::Mat curLargePatch, cv::Point2f curPatchLoc,  double curScale,
                                                            cv::Point2f refPoint, double refScale, Eigen::Matrix3d H) {

    // The patch to compute
    std::vector<double> patch (patchSize*patchSize, 0.0);

    // Let's iterate over all points in patch for point in image B and use homography to find them in image A
    int index = 0;
    for (int j = - halfPatchSize; j <  halfPatchSize + 1; j++) {
        for (int i = - halfPatchSize; i < halfPatchSize + 1; i++) {

            // Computation of the points of patch in image A taking into account the scale
            double x = refPoint.x + i * refScale;
            double y = refPoint.y + j * refScale;

            // Projecting (x,y) from image A into image B with  H
            Eigen::Vector3d pInA(x,y,1);
            Eigen::Vector3d pInB = H * pInA;
            pInB = pInB / pInB(2);


            // We only have a image part in image 1!
            // Therefore, we find the point w.r.t the patch location (center), scale that distance by octave's scale and just moves it to the center of stored patch
            int centerPos = curLargePatch.rows/ 2;
            pInB(0) = (pInB(0) - curPatchLoc.x)/curScale + centerPos;
            pInB(1) = (pInB(1) - curPatchLoc.y)/curScale + centerPos;


            // If warped outside of the stored large patch
            if ( pInB(0) < 0 || pInB(0) > curLargePatch.cols - 1 || pInB(1) < 0 || pInB(1) > curLargePatch.rows - 1)
                return std::vector<double>();

            // Finding the subpix value using bilinear interpolation
            const double xInt = int(pInB(0)), yInt = int(pInB(1));
            const double xSub = pInB(0) - xInt, ySub = pInB(1) - yInt;

            // From wiki: http://upload.wikimedia.org/math/9/b/4/9b4e1064436ecccd069ea238b656c063.png
            const double topLeft = (1.0 - xSub) * (1.0 - ySub);
            const double topRight = xSub * (1.0 - ySub);
            const double bottomLeft = (1.0 - xSub) * ySub;
            const double bottomRight = xSub * ySub;

            // Value of that patch
            patch[index] = topLeft * curLargePatch.at<uchar>(yInt, xInt) + topRight * curLargePatch.at<uchar>(yInt, xInt + 1) +
                          bottomLeft * curLargePatch.at<uchar>(yInt + 1, xInt) + bottomRight * curLargePatch.at<uchar>(yInt + 1, xInt + 1);

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