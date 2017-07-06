//
// Created by michalnowicki on 03.07.17.
//

#include "PatchRefinement.h"


Eigen::Vector2d PatchRefinement::computePatchDifference(std::vector<double> patch, std::vector<double> optimizedPatch,
                                       std::vector<Eigen::Vector2d> imageGradient) {
    float averagePatch = std::accumulate( patch.begin(), patch.end(), 0.0)/patch.size();
    float averageOptPatch = std::accumulate( optimizedPatch.begin(), optimizedPatch.end(), 0.0)/optimizedPatch.size();

    Eigen::Vector2d res = Eigen::Vector2d::Zero();
    for (int i=0;i<patch.size();i++) {
        float d = (optimizedPatch[i]-averageOptPatch) - (patch[i] - averagePatch);
        res += imageGradient[i] * (-d);
    }
    return res;
}

double PatchRefinement::computePatchDifference(std::vector<double> patch, std::vector<double> optimizedPatch) {
    float averagePatch = std::accumulate( patch.begin(), patch.end(), 0.0)/patch.size();
    float averageOptPatch = std::accumulate( optimizedPatch.begin(), optimizedPatch.end(), 0.0)/optimizedPatch.size();

    double res = 0;
    for (int i=0;i<patch.size();i++) {
        res += (optimizedPatch[i]-averageOptPatch) - (patch[i] - averagePatch);
    }
    return res;
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

cv::Mat PatchRefinement::ComputeHomography(cv::Mat Tba, cv::Mat n, double d, cv::Mat Ka, cv::Mat Kb)
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

std::vector<double> PatchRefinement::computePatch(cv::Mat img, double px, double py, int patchSize, Eigen::Matrix3d H) {

    std::vector<double> patch (patchSize*patchSize, 0.0);

    int halfPatchSize = (patchSize - 1)/2;
    int index = 0;
    for (int j = py - halfPatchSize; j < py + halfPatchSize + 1; j++) {
        for (int i = px - halfPatchSize; i < px + halfPatchSize + 1; i++) {

            // We find the position on the image1 based on (px,py) in image2 and H
            Eigen::Vector3d pIn2(i,j,1);
            Eigen::Vector3d pIn1 = H * pIn2;
            pIn1 = pIn1 / pIn1(2);

            const double xInt = int(pIn1(0)), yInt = int(pIn1(1));
            const double xSub = pIn1(0) - xInt, ySub = pIn1(1) - yInt;

            // From wiki: http://upload.wikimedia.org/math/9/b/4/9b4e1064436ecccd069ea238b656c063.png
            const double topLeft = (1.0 - xSub) * (1.0 - ySub);
            const double topRight = xSub * (1.0 - ySub);
            const double bottomLeft = (1.0 - xSub) * ySub;
            const double bottomRight = xSub * ySub;


            float value = topLeft * img.at<uchar>(yInt, xInt) + topRight * img.at<uchar>(yInt, xInt + 1) +
                          bottomLeft * img.at<uchar>(yInt + 1, xInt) + bottomRight * img.at<uchar>(yInt + 1, xInt + 1);
            patch[index] = value;

            index++;
        }
    }

    return patch;
}

std::vector<double> PatchRefinement::computePatchOnSubImage(cv::Mat img, int patchSize, Eigen::Matrix3d H, cv::Point2f img2kp, double scaleKp2, cv::Point2f img1kp, double scaleKp1) {

    std::vector<double> patch (patchSize*patchSize, 0.0);

    int halfPatchSize = (patchSize - 1)/2;
    int index = 0;



    for (int j = - halfPatchSize; j <  halfPatchSize + 1; j++) {
        for (int i = - halfPatchSize; i < halfPatchSize + 1; i++) {

            // The patch size is defined in detection octave -> it will be different size on full img
            double x = img2kp.x + i * scaleKp2;
            double y = img2kp.y + j * scaleKp2;

            // We find the position on the image1 based on (x,y) in image2 and H
            Eigen::Vector3d pIn2(x,y,1);
            Eigen::Vector3d pIn1 = H * pIn2;
            pIn1 = pIn1 / pIn1(2);

            // We only have a image part in image 1 !
            // Therefore, we find the point with respect to keypoint, scale that distance by octave's scale and just moves it to the center of stored patch
            int centerPos = (img.rows - 1) / 2;
            pIn1(0) = (pIn1(0) - img1kp.x)/scaleKp1 + centerPos;
            pIn1(1) = (pIn1(1) - img1kp.y)/scaleKp1 + centerPos;

            // Warped outside of the image
            if ( pIn1(0) < 0 || pIn1(0) > img.cols - 1 || pIn1(1) < 0 || pIn1(1) > img.rows - 1)
            {
               //    printf("\t Our subimage is too small to compute warped patch!\n");
                return std::vector<double>();
            }

            // We can do warping
            const double xInt = int(pIn1(0)), yInt = int(pIn1(1));
            const double xSub = pIn1(0) - xInt, ySub = pIn1(1) - yInt;

            // From wiki: http://upload.wikimedia.org/math/9/b/4/9b4e1064436ecccd069ea238b656c063.png
            const double topLeft = (1.0 - xSub) * (1.0 - ySub);
            const double topRight = xSub * (1.0 - ySub);
            const double bottomLeft = (1.0 - xSub) * ySub;
            const double bottomRight = xSub * ySub;


            float value = topLeft * img.at<uchar>(yInt, xInt) + topRight * img.at<uchar>(yInt, xInt + 1) +
                          bottomLeft * img.at<uchar>(yInt + 1, xInt) + bottomRight * img.at<uchar>(yInt + 1, xInt + 1);
            patch[index] = value;

            index++;
        }
    }

    return patch;
}

std::vector<double> PatchRefinement::computePatch(cv::Mat img, double x, double y, int size, int patchSize, std::vector<Eigen::Vector2d> gradient, std::vector<Eigen::Vector2d> &gd) {

    const double xInt = int(x), yInt = int(y);
    const double xSub = x - xInt, ySub = y - yInt;

    // From wiki: http://upload.wikimedia.org/math/9/b/4/9b4e1064436ecccd069ea238b656c063.png
    const double topLeft = (1.0 - xSub) * (1.0 - ySub);
    const double topRight = xSub * (1.0 - ySub);
    const double bottomLeft = (1.0 - xSub) * ySub;
    const double bottomRight = xSub * ySub;

    std::vector<double> patch (patchSize*patchSize, 0.0);

    int halfPatchSize = (patchSize - 1)/2;
    int index = 0;
    for (int j = yInt - halfPatchSize; j < yInt + halfPatchSize + 1; j++) {
        for (int i = xInt - halfPatchSize; i < xInt + halfPatchSize + 1; i++) {

            patch[index] = topLeft * img.at<uchar>(j, i) + topRight * img.at<uchar>(j, i + 1) +
                           bottomLeft * img.at<uchar>(j + 1, i) + bottomRight * img.at<uchar>(j + 1, i + 1);

            gd[index] = topLeft * gradient[j*size+i] + topRight * gradient[j*size+i+1] +
                        bottomLeft * gradient[(j+1)*size+i] + bottomRight * gradient[(j+1)*size+i+1];

            index++;
        }
    }

    return patch;
}

void PatchRefinement::computeImageGradient(cv::Mat &img, std::vector<Eigen::Vector2d> &gradient, Eigen::Matrix2d &HessianInv) {
    int size = img.rows;
    gradient = std::vector<Eigen::Vector2d>(size*size, Eigen::Vector2d::Zero());
    Eigen::Matrix2d Hessian = Eigen::Matrix2d::Zero();
    for (int i = 1; i<size-1;i++)
    {
        for (int j = 1;j<size-1;j++) {
            Eigen::Vector2d Jacobian;
            Jacobian[0] = 0.5 * (img.at<uchar>(i,j+1) - img.at<uchar>(i,j-1));
            Jacobian[1] = 0.5 * (img.at<uchar>(i+1,j) - img.at<uchar>(i-1,j));
            gradient[j+i*size] = Jacobian;

            Hessian += Jacobian * Jacobian.transpose();
        }
    }
    HessianInv = Hessian.inverse();
};


void PatchRefinement::testOptimization(cv::Mat img, int patchSize, Eigen::Vector2d originalPoint, Eigen::Vector2d testPoint, double minimalStep) {
    std::cout <<"Original patch for " << originalPoint.transpose() << std::endl;
    std::cout <<"Test patch for " << testPoint.transpose() << std::endl;

    // Compute gradient
    std::vector<Eigen::Vector2d> gradient;
    Eigen::Matrix2d HessianInv;
    computeImageGradient(img, gradient, HessianInv);

    // Subpix gradient - will be used later
    std::vector<Eigen::Vector2d> subPosGradient(patchSize*patchSize, Eigen::Vector2d::Zero());

    // Original patch
    std::vector<double> originalPatch = computePatch(img, originalPoint[0], originalPoint[1], img.rows, patchSize, gradient, subPosGradient);


    // Lets start the optimization from provided position
    float optX = testPoint[0], optY = testPoint[1];

    // We perform 100 iterations or until stop condition is met
    for (int i = 0; i < 100; i++) {
        std::cout <<"Iteration: " << i << " Position: " << optX << ", " << optY << std::endl;

        // Compute new patch
        std::vector<double> movedPatch = computePatch(img, optX, optY, img.rows, patchSize, gradient, subPosGradient);

        // Estimate residuals
        Eigen::Vector2d res = computePatchDifference(originalPatch, movedPatch, subPosGradient);

        // Obtained error
        std::cout << "Error: " << std::endl << res.transpose() << std::endl;

        // Step based on Hessian
        Eigen::Vector2d step = HessianInv * res;

        // Stop condition
        if (step.squaredNorm() < minimalStep)
            break;

        // We do the step
        optX = optX + step[0];
        optY = optY + step[1];
    }

    std::cout << "Final position : " << optX << ", " << optY << std::endl;
}


void PatchRefinement::testHomography(cv::Mat image2, cv::Mat point3DInImg1, cv::Mat K1, cv::Mat K2, cv::Mat R, cv::Mat t, int patchSize) {

    std::cout << "testHommography()" << std::endl;
    std::cout << "-> point3DInImg1 = " << std::endl << point3DInImg1 << std::endl;
    std::cout << "-> K1 = " << std::endl << K1 << std::endl;
    std::cout << "-> K2 = " << std::endl << K2 << std::endl;
    std::cout << "-> R = " << std::endl << R << std::endl;
    std::cout << "-> t = " << std::endl << t << std::endl;
    std::cout << "----------------------"<< std::endl;

    // We project point into image 1
    cv::Mat projPInImg1 = K1 * point3DInImg1;
    projPInImg1 = normalize2D(projPInImg1);
    std::cout << "-> projPInImg1 = " << std::endl << projPInImg1 << std::endl;

    // We find the 3D position of point in image 2
    cv::Mat point3DInImg2 = R.inv() * (point3DInImg1 - t);
    std::cout << "-> point3DInImg2 = " << std::endl << point3DInImg2 << std::endl;

    // We project point into image 2
    cv::Mat projPInImg2 = K2 * point3DInImg2;
    projPInImg2 = normalize2D(projPInImg2);
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
    double d = getDistanceToPlane(point3DInImg1, n);
    std::cout << "-> d = " << d << std::endl;

    // We compute the homography
    // Pose of c.s. A in B
    cv::Mat Tba = T2w.inv() * T1w;
    cv::Mat H = ComputeHomography(Tba, n, d, K1, K2);
    std::cout << "-> H = " << std::endl << H << std::endl;

    //// Let's verify the homography
    std::cout << "----------------------"<< std::endl;
    std::cout << "Verifying homography" << std::endl;

    double diff = cv::norm(projPInImg2 - normalize2D(H*projPInImg1));
    double diff2 = cv::norm(projPInImg1 - normalize2D(H.inv()*projPInImg2));
    std :: cout << "Img1: " << projPInImg1.t() << " " << normalize2D(H.inv()*projPInImg2) << std::endl;
    std :: cout << "Img2: " << projPInImg2.t() << " " << normalize2D(H*projPInImg1) << std::endl;

    if ( diff < 0.0001 && diff2 < 0.0001)
        std::cout << "\t Homography is OK!" << std::endl;
    else
        std::cout << "\t Homography is WROOOONG!" << std::endl;

    // This is how probably image 1 looks like! (Based on H) - this is only for TEST
    cv::Mat image1;
    cv::warpPerspective(image2, image1, H, image2.size(), cv::INTER_LINEAR, cv::BORDER_CONSTANT, cv::Scalar());


    //// Let's verify the homography
    std::cout << "----------------------"<< std::endl;
    std::cout << "Verifying warping patches" << std::endl;

    // Patch parameters
    int halfPatchSize = (patchSize - 1) / 2;

    // Patch we observe in image 2
    int x = int(projPInImg2.at<float>(0));
    int y = int(projPInImg2.at<float>(1));
    std::vector<double> patchImg2 = computePatch(image2, x, y, patchSize, Eigen::Matrix3d::Identity());

    // Patch we observe in image 1
    Eigen::Matrix3d Heigen;
    cv::cv2eigen(H, Heigen);
    std::vector<double> patchImg1 = computePatch(image1, x, y, patchSize, Heigen.inverse());

    std::cout << "Patch 1!" << std::endl;
    printPatch(patchImg1, patchSize);
//    showPatch(patchImg1, patchSize, "Patch 1");

    std::cout << "Patch 2!" << std::endl;
    printPatch(patchImg2, patchSize);
//    showPatch(patchImg2, patchSize, "Patch 2");

    std::cout << "Patch diff: " << computePatchDifference(patchImg1, patchImg2) << std::endl;

    // Verify original images
//    cv::Mat toShow(image2.rows, image2.cols+image2.cols, CV_8U);
//    cv::Mat left(toShow, cv::Rect(0, 0, image2.cols, image2.rows));
//    image1.copyTo(left);
//    cv::Mat right(toShow, cv::Rect(image2.cols, 0, image2.cols, image2.rows));
//    image2.copyTo(right);
//    cv::cvtColor(toShow, toShow, cv::COLOR_GRAY2BGR);

//    halfPatchSize = 100;
//
//    int cx = projPInImg2.at<float>(0), cy = projPInImg2.at<float>(1);
//    cv::circle(toShow, cv::Point2f(cx+image2.cols, cy), 3, cv::Scalar(0,255,0), 2);
//    for (int j = cy - halfPatchSize; j < cy + halfPatchSize + 1; j=j+2*halfPatchSize-1)
//        for (int i = cx - halfPatchSize; i < cx + halfPatchSize + 1; i=i+2*halfPatchSize-1)
//            cv::circle(toShow, cv::Point2f(image2.cols+i, j), 2, cv::Scalar(0,0,255), 2);
//
//
//    cv::circle(toShow, cv::Point2f(projPInImg1.at<float>(0), projPInImg1.at<float>(1)), 3, cv::Scalar(0,255,0), 2);
//    for (int j = cy - halfPatchSize; j < cy + halfPatchSize + 1; j=j+2*halfPatchSize-1) {
//        for (int i = cx - halfPatchSize; i < cx + halfPatchSize + 1; i=i+2*halfPatchSize-1) {
//
//            Eigen::Vector3f pIn2(i, j, 1);
//            Eigen::Vector3f pIn1 = Heigen * pIn2;
//            pIn1 = pIn1 / pIn1(2);
//
//            cv::circle(toShow, cv::Point2f(pIn1(0), pIn1(1)), 2, cv::Scalar(0,0,255), 2);
//        }
//    }

//    cv::imshow("Image 1 and Image 2", toShow);
//    cv::waitKey(0);
}