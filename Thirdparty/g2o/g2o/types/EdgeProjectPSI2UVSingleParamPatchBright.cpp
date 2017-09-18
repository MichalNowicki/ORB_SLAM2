//
// Created by michalnowicki on 06.09.17.
//


#include "EdgeProjectPSI2UVSingleParamPatchBright.h"


#include "../core/factory.h"
#include "../stuff/macros.h"

namespace g2o {
    using namespace std;

    bool EdgeProjectPSI2UVSingleParamPatchBright::write(std::ostream &os) const {
        os << _cam->id() << " ";
        for (int i = 0; i < errorSize; i++) {
            os << measurement()[i] << " ";
        }

        for (int i = 0; i < errorSize; i++)
            for (int j = i; j < errorSize; j++) {
                os << " " << information()(i, j);
            }
        return os.good();
    }

    bool EdgeProjectPSI2UVSingleParamPatchBright::read(std::istream &is) {
        int paramId;
        is >> paramId;
        setParameterId(0, paramId);

        for (int i = 0; i < errorSize; i++) {
            is >> _measurement[i];
        }
        for (int i = 0; i < errorSize; i++)
            for (int j = i; j < errorSize; j++) {
                is >> information()(i, j);
                if (i != j)
                    information()(j, i) = information()(i, j);
            }
        return true;
    }


    void EdgeProjectPSI2UVSingleParamPatchBright::computeError() {
        const VertexSBAPointInvD *pointInvD = static_cast<const VertexSBAPointInvD *>(_vertices[0]);
        const VertexSE3ExpmapBright *T_p_from_world = static_cast<const VertexSE3ExpmapBright *>(_vertices[1]);
        const VertexSE3ExpmapBright *T_anchor_from_world = static_cast<const VertexSE3ExpmapBright *>(_vertices[2]);
        const CameraParameters *cam = static_cast<const CameraParameters *>(parameter(0));


        SE3QuatBright T_p_est = T_p_from_world->estimate();
        SE3QuatBright T_anchor_est = T_anchor_from_world->estimate();


        double cx = cam->principle_point[0], cy = cam->principle_point[1];
        double fx = cam->focal_length_x, fy = cam->focal_length_y;


//        Matrix<double, 9, 1, Eigen::ColMajor> computedError;
        Matrix<double, 25, 1, Eigen::ColMajor> computedError;


        for (int i=0;i<neighbours.size();i++)
        {
            // Getting the patch value in anchor
            int refU = largePatchCenter + neighbours[i].first, refV = largePatchCenter + neighbours[i].second;
            double refValue = largePatchAnchor[refV * largePatchStride + refU];

            // Getting the projected point in obs
            Eigen::Vector3d pointInFirst;
            pointInFirst[2] = 1. / pointInvD->estimate();
            pointInFirst[0] = (pointInvD->u0 - cx + neighbours[i].first * pointAnchorScale) * pointInFirst[2] / fx;
            pointInFirst[1] = (pointInvD->v0 - cy + neighbours[i].second * pointAnchorScale) * pointInFirst[2] / fy;
            Eigen::Vector3d pointInGlobal = T_anchor_est.se3quat.inverse().map(pointInFirst);
            Vector2d projectedPoint = cam->cam_map(T_p_est.se3quat.map(pointInGlobal));

            // Find where the projected point is on stored patch
            double patchOffsetU = (projectedPoint[0] - _measurement[0]) / pointObsScale;
            double patchOffsetV = (projectedPoint[1] - _measurement[1]) / pointObsScale;

            double obsU = largePatchCenter + patchOffsetU;
            double obsV = largePatchCenter + patchOffsetV;

            // The obs value will not be integer
            const double xInt = int(obsU), yInt = int(obsV);
            const double xSub = obsU - xInt, ySub = obsV - yInt;

            const double topLeft = (1.0 - xSub) * (1.0 - ySub);
            const double topRight = xSub * (1.0 - ySub);
            const double bottomLeft = (1.0 - xSub) * ySub;
            const double bottomRight = xSub * ySub;


            if (yInt < 0 || xInt < 0 || yInt + 1 > largePatchStride || xInt + 1 > largePatchStride)
            {
//                std::cout << "Projected: " << projectedPoint[0] << " " << projectedPoint[1] << " InvDepth: " << pointInvD->estimate() << std::endl;
//                std::cout << "Outside patch: " << yInt << " " << xInt << "  PatchOffset: " << patchOffsetU << " " << patchOffsetV << std::endl;
//                std::cout << "\t" << projectedPoint[0] << " " << _measurement[0] << " " << projectedPoint[1] << " " << _measurement[1] << std::endl;

//                for (int j=0;j<9;j++)
                for (int j=0;j<25;j++)
                    computedError(j,0) = 255;
                break;
            }

            double obsValue = topLeft * largePatchObs[yInt * largePatchStride + xInt] +
                              topRight * largePatchObs[yInt * largePatchStride + xInt + 1] +
                              bottomLeft * largePatchObs[(yInt + 1) * largePatchStride + xInt] +
                              bottomRight * largePatchObs[(yInt+1) * largePatchStride + xInt + 1];

            // a_1 / a_2 * (I_1 - b_1) - (I_2 - b_2)
//            computedError(i,0) = T_anchor_est.a / T_p_est.a * (refValue - T_anchor_est.b) - (obsValue - T_p_est.b);


            computedError(i,0) = exp(T_p_est.a) / exp(T_anchor_est.a)  * (refValue - T_anchor_est.b) - (obsValue - T_p_est.b);

//            computedError(i,0) = refValue - obsValue;

//            if (refValue - obsValue > 10) {
//                std::cout << "Projected and measured: (" << _measurement[0] << ", " << _measurement[1] << ") vs ("
//                          << projectedPoint[0] << ". " << projectedPoint[1] << ")" << std::endl;
//                std::cout << "computedError(i,0) = " << refValue << " " << obsValue << std::endl;
//            }
        }

        _error = computedError;
    }

    inline Eigen::Matrix<double, 1, 2> EdgeProjectPSI2UVSingleParamPatchBright::d_inten_d_proj(const double obsU, const double obsV) {

        // The obs value will not be integer
        const double xInt = int(obsU), yInt = int(obsV);
        const double xSub = obsU - xInt, ySub = obsV - yInt;

        const double topLeft = (1.0 - xSub) * (1.0 - ySub);
        const double topRight = xSub * (1.0 - ySub);
        const double bottomLeft = (1.0 - xSub) * ySub;
        const double bottomRight = xSub * ySub;


        Eigen::Matrix<double, 1, 2> G;
        for (int i=0;i<2;i++) {
            G(0,i) = topLeft * largePatchObsGradient[yInt * largePatchStride + xInt][i] +
                     topRight * largePatchObsGradient[yInt * largePatchStride + xInt + 1][i] +
                     bottomLeft * largePatchObsGradient[(yInt+1) * largePatchStride + xInt][i] +
                     bottomRight * largePatchObsGradient[(yInt+1) * largePatchStride + xInt + 1][i];
        }
        return G;
    };

    inline Matrix<double, 2, 3, Eigen::ColMajor>
    EdgeProjectPSI2UVSingleParamPatchBright::d_proj_d_y(const double &fx, const double &fy, const Vector3D &xyz) {
        double z_sq = xyz[2] * xyz[2];
        Matrix<double, 2, 3, Eigen::ColMajor> J;
        J << fx / xyz[2], 0, -(fx * xyz[0]) / z_sq,
                0, fy / xyz[2], -(fy * xyz[1]) / z_sq;
        return J;
    }

    inline Matrix<double, 3, 6, Eigen::ColMajor> EdgeProjectPSI2UVSingleParamPatchBright::d_expy_d_y(const Vector3D &y) {
        Matrix<double, 3, 6, Eigen::ColMajor> J;
        J.topLeftCorner<3, 3>() = -skew(y);
        J.bottomRightCorner<3, 3>().setIdentity();

        return J;
    }

    inline Matrix<double, 3, 1, Eigen::ColMajor>
    EdgeProjectPSI2UVSingleParamPatchBright::d_Tinvpsi_d_psi(const SE3Quat &T, const Vector3D &psi) {
        Matrix3D R = T.rotation().toRotationMatrix();
        Vector3D x = invert_depth(psi);
        Matrix<double, 3, 1, Eigen::ColMajor> J;
        J = -R * x;
        J *= 1. / psi.z();
        return J;
    }

//    void EdgeProjectPSI2UVSingleParamPatchBright::linearizeOplus() {
//
//        // Estiamted values
//        VertexSBAPointInvD *pointInvD = static_cast<VertexSBAPointInvD *>(_vertices[0]);
//        VertexSE3ExpmapBright *vpose = static_cast<VertexSE3ExpmapBright *>(_vertices[1]);
//        VertexSE3ExpmapBright *vanchor = static_cast<VertexSE3ExpmapBright *>(_vertices[2]);
//
//        SE3QuatBright vpose_est = vpose->estimate();
//        SE3QuatBright vanchor_est = vanchor->estimate();
//
//        // Camera parameters
//        const CameraParameters *cam
//                = static_cast<const CameraParameters *>(parameter(0));
//        double cx = cam->principle_point[0], cy = cam->principle_point[1];
//        double fx = cam->focal_length_x, fy = cam->focal_length_y;
//
//
//        // Empty jacobians
//        _jacobianOplus[0] = Matrix<double, 9, 1, Eigen::ColMajor>::Zero();
//        _jacobianOplus[1] = Matrix<double, 9, 8, Eigen::ColMajor>::Zero();
//        _jacobianOplus[2] = Matrix<double, 9, 8, Eigen::ColMajor>::Zero();
//
//
//        // Getting current pose estimate
//        SE3Quat T_cw = vpose_est.se3quat;
//
//        // Getting the anchor pose estimate
//        SE3Quat T_aw = vanchor_est.se3quat;
//
//        // From anchor to current
//        SE3Quat T_ca = T_cw * T_aw.inverse();
//
//        // For all points in neighbourhood
//        for (int i=0;i<9;i++) {
//            // Getting the patch value in anchor
//            int refU = largePatchCenter + neighbours[i].first, refV = largePatchCenter + neighbours[i].second;
//            double refValue = largePatchAnchor[refV * largePatchStride + refU];
//            double Ianchor = refValue - vanchor_est.b;
//
//            // Getting the projected point in obs
//            Eigen::Vector3d pointInFirst;
//            pointInFirst[2] = 1. / pointInvD->estimate();
//            pointInFirst[0] = (pointInvD->u0 - cx + neighbours[i].first * pointAnchorScale) * pointInFirst[2] / fx;
//            pointInFirst[1] = (pointInvD->v0 - cy + neighbours[i].second * pointAnchorScale) * pointInFirst[2] / fy;
//
//            // Global point position
//            Eigen::Vector3d pointInGlobal = T_aw.inverse().map(pointInFirst);
//
//            // 3D point in obs
//            Vector3D pointInObs = T_cw.map(pointInGlobal);
//
//            // 2D projected point in anchor
//            Vector2d projectedPoint = cam->cam_map(pointInObs);
//
//            // Point in anchor in inverse depth parametrization
//            Vector3D psi_a = invert_depth(pointInFirst);
//
//            // Jacobian of camera
//            Matrix<double, 2, 3, Eigen::ColMajor> Jcam
//                    = d_proj_d_y(cam->focal_length_x, cam->focal_length_y, pointInObs);
//
//            // Find where the projected point is on stored patch
//            double patchOffsetU = (projectedPoint[0] - _measurement[0]) / pointObsScale;
//            double patchOffsetV = (projectedPoint[1] - _measurement[1]) / pointObsScale;
//
//            // Observation on current frame in largePatch CS
//            double obsU = largePatchCenter + patchOffsetU;
//            double obsV = largePatchCenter + patchOffsetV;
//
//            // Image gradient
//            Matrix<double, 1, 2> Ji = d_inten_d_proj(obsU, obsV);
//
//            // Jacobians of point, observation pose and anchor pose
//            _jacobianOplus[0].row(i) = - Ji * Jcam * d_Tinvpsi_d_psi(T_ca, psi_a);
//
//            _jacobianOplus[1].block<1,6>(i,0) = - Ji * Jcam * d_expy_d_y(pointInObs);
//            _jacobianOplus[1](i,6) = exp (vpose_est.a) / exp( vanchor_est.a) * (Ianchor);
//            _jacobianOplus[1](i,7) = 1;
//
//            _jacobianOplus[2].block<1,6>(i,0) = Ji * Jcam * T_ca.rotation().toRotationMatrix() * d_expy_d_y(pointInFirst);
//            _jacobianOplus[2](i,6) = - exp (vpose_est.a) / exp( vanchor_est.a) * (Ianchor);
//            _jacobianOplus[2](i,7) = - exp (vpose_est.a) / exp( vanchor_est.a);
//
////            std::cout << "Ji: " << std::endl<< Ji << std::endl;
////            std::cout << "Jcam: " << std::endl<< Jcam << std::endl;
////            std::cout << "d_Tinvpsi_d_psi: "<< std::endl << d_Tinvpsi_d_psi(T_ca, psi_a) << std::endl;
////            std::cout << "d_expy_d_y(pointInObs): " << std::endl<< d_expy_d_y(pointInObs) << std::endl;
////            std::cout << "d_expy_d_y(pointInFirst): " << std::endl << d_expy_d_y(pointInFirst) << std::endl;
//        }
//
////        std::cout << "_jacobianOplus[0]: " << std::endl << _jacobianOplus[0] << std::endl;
////
////        std::cout << "_jacobianOplus[1]: " << std::endl << _jacobianOplus[1] << std::endl;
////
////        std::cout << "_jacobianOplus[1]: " << std::endl << _jacobianOplus[2] << std::endl;
//
//    }

    bool EdgeProjectPSI2UVSingleParamPatchBright::isDepthPositive() {

        const VertexSBAPointInvD *pointInvD = static_cast<const VertexSBAPointInvD *>(_vertices[0]);
        const VertexSE3ExpmapBright *T_p_from_world = static_cast<const VertexSE3ExpmapBright *>(_vertices[1]);
        const VertexSE3ExpmapBright *T_anchor_from_world = static_cast<const VertexSE3ExpmapBright *>(_vertices[2]);
        const CameraParameters *cam = static_cast<const CameraParameters *>(parameter(0));

        SE3QuatBright T_p_est = T_p_from_world->estimate();
        SE3QuatBright T_anchor_est = T_anchor_from_world->estimate();

        double cx = cam->principle_point[0], cy = cam->principle_point[1];
        double fx = cam->focal_length_x, fy = cam->focal_length_y;

        Eigen::Vector3d pointInAnchor;
        pointInAnchor[2] = 1. / pointInvD->estimate();
        pointInAnchor[0] = (pointInvD->u0 - cx) * pointInAnchor[2] / fx;
        pointInAnchor[1] = (pointInvD->v0 - cy) * pointInAnchor[2] / fy;

        Eigen::Vector3d pointInGlobal = T_anchor_est.se3quat.inverse().map(pointInAnchor);
        Eigen::Vector3d pointInObs = T_p_est.se3quat.map(pointInGlobal);

        return pointInAnchor(2) > 0.0 && pointInObs(2) > 0.0;
    }


}