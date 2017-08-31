//
// Created by michalnowicki on 30.08.17.
//

#include "EdgeProjectPSI2UV.h"

#include "../core/factory.h"
#include "../stuff/macros.h"


namespace g2o {
    using namespace std;


    bool EdgeProjectPSI2UV::write(std::ostream &os) const {
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

    bool EdgeProjectPSI2UV::read(std::istream &is) {
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

    void EdgeProjectPSI2UV::computeError() {
        const VertexSBAPointXYZ *psi = static_cast<const VertexSBAPointXYZ *>(_vertices[0]);
        const VertexSE3Expmap *T_p_from_world = static_cast<const VertexSE3Expmap *>(_vertices[1]);
        const VertexSE3Expmap *T_anchor_from_world = static_cast<const VertexSE3Expmap *>(_vertices[2]);
        const CameraParameters *cam = static_cast<const CameraParameters *>(parameter(0));

        Vector2D obs(_measurement);
        _error = obs - cam->cam_map(T_p_from_world->estimate()
                                    * T_anchor_from_world->estimate().inverse()
                                    * invert_depth(psi->estimate()));
    }

    inline Matrix<double, 2, 3, Eigen::ColMajor>
    EdgeProjectPSI2UV::d_proj_d_y(const double &fx, const double &fy, const Vector3D &xyz) {
        double z_sq = xyz[2] * xyz[2];
        Matrix<double, 2, 3, Eigen::ColMajor> J;
        J << fx / xyz[2], 0, -(fx * xyz[0]) / z_sq,
                0, fy / xyz[2], -(fy * xyz[1]) / z_sq;
        return J;
    }

    inline Matrix<double, 3, 6, Eigen::ColMajor> EdgeProjectPSI2UV::d_expy_d_y(const Vector3D &y) {
        Matrix<double, 3, 6, Eigen::ColMajor> J;
        J.topLeftCorner<3, 3>() = -skew(y);
        J.bottomRightCorner<3, 3>().setIdentity();

        return J;
    }

    inline Matrix3D EdgeProjectPSI2UV::d_Tinvpsi_d_psi(const SE3Quat &T, const Vector3D &psi) {
        Matrix3D R = T.rotation().toRotationMatrix();
        Vector3D x = invert_depth(psi);
        Vector3D r1 = R.col(0);
        Vector3D r2 = R.col(1);
        Matrix3D J;
        J.col(0) = r1;
        J.col(1) = r2;
        J.col(2) = -R * x;
        J *= 1. / psi.z();
        return J;
    }

    void EdgeProjectPSI2UV::linearizeOplus() {

        VertexSBAPointXYZ *vpoint = static_cast<VertexSBAPointXYZ *>(_vertices[0]);
        Vector3D psi_a = vpoint->estimate();
        VertexSE3Expmap *vpose = static_cast<VertexSE3Expmap *>(_vertices[1]);
        SE3Quat T_cw = vpose->estimate();
        VertexSE3Expmap *vanchor = static_cast<VertexSE3Expmap *>(_vertices[2]);
        const CameraParameters *cam
                = static_cast<const CameraParameters *>(parameter(0));

        SE3Quat A_aw = vanchor->estimate();
        SE3Quat T_ca = T_cw * A_aw.inverse();
        Vector3D x_a = invert_depth(psi_a);
        Vector3D y = T_ca * x_a;
        Matrix<double, 2, 3, Eigen::ColMajor> Jcam
                = d_proj_d_y(cam->focal_length_x, cam->focal_length_y, y);
        _jacobianOplus[0] = -Jcam * d_Tinvpsi_d_psi(T_ca, psi_a);
        _jacobianOplus[1] = -Jcam * d_expy_d_y(y);
        _jacobianOplus[2] = Jcam * T_ca.rotation().toRotationMatrix() * d_expy_d_y(x_a);
    }

    bool EdgeProjectPSI2UV::isDepthPositive() {
        const VertexSBAPointXYZ *pointInvD = static_cast<const VertexSBAPointXYZ *>(_vertices[0]);
        const VertexSE3Expmap *firstPose = static_cast<const VertexSE3Expmap *>(_vertices[1]);
        const VertexSE3Expmap *secondPose = static_cast<const VertexSE3Expmap *>(_vertices[2]);

        Vector3D pointInAnchor = invert_depth(pointInvD->estimate());
        Vector3D pointInGlobal = secondPose->estimate().inverse().map(pointInAnchor);
        Vector3D pointInFirst = firstPose->estimate().map(pointInGlobal);

        return pointInAnchor(2) > 0.0 && pointInFirst(2) > 0.0;
    }

}
