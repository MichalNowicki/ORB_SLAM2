//
// Created by michalnowicki on 30.08.17.
//

#include "EdgeProjectPSI2UVSingleParam.h"

#include "../core/factory.h"
#include "../stuff/macros.h"

namespace g2o {
    using namespace std;

    bool EdgeProjectPSI2UVSingleParam::write(std::ostream &os) const {
        os << _cam->id() << " ";
        for (int i = 0; i < 2; i++) {
            os << measurement()[i] << " ";
        }

        for (int i = 0; i < 2; i++)
            for (int j = i; j < 2; j++) {
                os << " " << information()(i, j);
            }
        return os.good();
    }

    bool EdgeProjectPSI2UVSingleParam::read(std::istream &is) {
        int paramId;
        is >> paramId;
        setParameterId(0, paramId);

        for (int i = 0; i < 2; i++) {
            is >> _measurement[i];
        }
        for (int i = 0; i < 2; i++)
            for (int j = i; j < 2; j++) {
                is >> information()(i, j);
                if (i != j)
                    information()(j, i) = information()(i, j);
            }
        return true;
    }


    void EdgeProjectPSI2UVSingleParam::computeError() {

        const VertexSBAPointInvD *pointInvD = static_cast<const VertexSBAPointInvD *>(_vertices[0]);
        const VertexSE3Expmap *T_p_from_world = static_cast<const VertexSE3Expmap *>(_vertices[1]);
        const VertexSE3Expmap *T_anchor_from_world = static_cast<const VertexSE3Expmap *>(_vertices[2]);
        const CameraParameters *cam = static_cast<const CameraParameters *>(parameter(0));

        double cx = cam->principle_point[0], cy = cam->principle_point[1];
        double fx = cam->focal_length_x, fy = cam->focal_length_y;

        Eigen::Vector3d pointInFirst;
        pointInFirst[2] = 1. / pointInvD->estimate();
        pointInFirst[0] = (pointInvD->u0 - cx) * pointInFirst[2] / fx;
        pointInFirst[1] = (pointInvD->v0 - cy) * pointInFirst[2] / fy;

        Eigen::Vector3d pointInGlobal = T_anchor_from_world->estimate().inverse().map(pointInFirst);

        Vector2d obs(_measurement);
        Vector2d projectedPoint = cam->cam_map(T_p_from_world->estimate().map(pointInGlobal));
        _error = obs - projectedPoint;
    }

    inline Matrix<double, 2, 3, Eigen::ColMajor>
    EdgeProjectPSI2UVSingleParam::d_proj_d_y(const double &fx, const double &fy, const Vector3D &xyz) {
        double z_sq = xyz[2] * xyz[2];
        Matrix<double, 2, 3, Eigen::ColMajor> J;
        J << fx / xyz[2], 0, -(fx * xyz[0]) / z_sq,
                0, fy / xyz[2], -(fy * xyz[1]) / z_sq;
        return J;
    }

    inline Matrix<double, 3, 6, Eigen::ColMajor> EdgeProjectPSI2UVSingleParam::d_expy_d_y(const Vector3D &y) {
        Matrix<double, 3, 6, Eigen::ColMajor> J;
        J.topLeftCorner<3, 3>() = -skew(y);
        J.bottomRightCorner<3, 3>().setIdentity();

        return J;
    }

    inline Matrix<double, 3, 1, Eigen::ColMajor>
    EdgeProjectPSI2UVSingleParam::d_Tinvpsi_d_psi(const SE3Quat &T, const Vector3D &psi) {
        Matrix3D R = T.rotation().toRotationMatrix();
        Vector3D x = invert_depth(psi);
        Matrix<double, 3, 1, Eigen::ColMajor> J;
        J = -R * x;
        J *= 1. / psi.z();
        return J;
    }

    void EdgeProjectPSI2UVSingleParam::linearizeOplus() {

        VertexSBAPointInvD *pointInvD = static_cast<VertexSBAPointInvD *>(_vertices[0]);

        VertexSE3Expmap *vpose = static_cast<VertexSE3Expmap *>(_vertices[1]);
        SE3Quat T_cw = vpose->estimate();
        VertexSE3Expmap *vanchor = static_cast<VertexSE3Expmap *>(_vertices[2]);
        const CameraParameters *cam
                = static_cast<const CameraParameters *>(parameter(0));

        double cx = cam->principle_point[0], cy = cam->principle_point[1];
        double fx = cam->focal_length_x, fy = cam->focal_length_y;

        Eigen::Vector3d x_a;
        x_a[2] = 1. / pointInvD->estimate();
        x_a[0] = (pointInvD->u0 - cx) * x_a[2] / fx;
        x_a[1] = (pointInvD->v0 - cy) * x_a[2] / fy;

        Vector3D psi_a = invert_depth(x_a);

        SE3Quat A_aw = vanchor->estimate();
        SE3Quat T_ca = T_cw * A_aw.inverse();

        Vector3D y = T_ca * x_a;
        Matrix<double, 2, 3, Eigen::ColMajor> Jcam
                = d_proj_d_y(cam->focal_length_x, cam->focal_length_y, y);
        _jacobianOplus[0] = -Jcam * d_Tinvpsi_d_psi(T_ca, psi_a);
        _jacobianOplus[1] = -Jcam * d_expy_d_y(y);
        _jacobianOplus[2] = Jcam * T_ca.rotation().toRotationMatrix() * d_expy_d_y(x_a);
    }

    bool EdgeProjectPSI2UVSingleParam::isDepthPositive() {

        const VertexSBAPointInvD *pointInvD = static_cast<const VertexSBAPointInvD *>(_vertices[0]);
        const VertexSE3Expmap *T_p_from_world = static_cast<const VertexSE3Expmap *>(_vertices[1]);
        const VertexSE3Expmap *T_anchor_from_world = static_cast<const VertexSE3Expmap *>(_vertices[2]);
        const CameraParameters *cam = static_cast<const CameraParameters *>(parameter(0));

        double cx = cam->principle_point[0], cy = cam->principle_point[1];
        double fx = cam->focal_length_x, fy = cam->focal_length_y;

        Eigen::Vector3d pointInAnchor;
        pointInAnchor[2] = 1. / pointInvD->estimate();
        pointInAnchor[0] = (pointInvD->u0 - cx) * pointInAnchor[2] / fx;
        pointInAnchor[1] = (pointInvD->v0 - cy) * pointInAnchor[2] / fy;

        Eigen::Vector3d pointInGlobal = T_anchor_from_world->estimate().inverse().map(pointInAnchor);
        Eigen::Vector3d pointInObs = T_p_from_world->estimate().map(pointInGlobal);

        return pointInAnchor(2) > 0.0 && pointInObs(2) > 0.0;
    }


}