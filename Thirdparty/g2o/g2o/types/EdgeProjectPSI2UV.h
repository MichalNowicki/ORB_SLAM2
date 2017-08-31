//
// Created by michalnowicki on 30.08.17.
//

#ifndef ORB_SLAM2_EDGEPROJECTPSI2UV_H
#define ORB_SLAM2_EDGEPROJECTPSI2UV_H

#include "../core/base_vertex.h"
#include "../core/base_binary_edge.h"
#include "../core/base_unary_edge.h"
#include "../core/base_multi_edge.h"
#include "../types/types_six_dof_expmap.h"
#include "se3_ops.h"
#include "se3quat.h"
#include "types_sba.h"
#include <Eigen/Geometry>

namespace g2o {
    using namespace std;

    class EdgeProjectPSI2UV : public g2o::BaseMultiEdge<2, Vector2D> {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        EdgeProjectPSI2UV() {
            resizeParameters(1);
            installParameter(_cam, 0);
        }

        virtual bool read(std::istream &is);

        virtual bool write(std::ostream &os) const;

        void computeError();

        virtual void linearizeOplus();

        inline Matrix<double, 2, 3, Eigen::ColMajor>
        d_proj_d_y(const double &fx, const double &fy, const Vector3D &xyz);

        inline Matrix<double, 3, 6, Eigen::ColMajor> d_expy_d_y(const Vector3D &y);

        inline Matrix3D d_Tinvpsi_d_psi(const SE3Quat &T, const Vector3D &psi);

        bool isDepthPositive();

        CameraParameters *_cam;

    private:
        inline Vector3D invert_depth(const Vector3D &x) {
            return unproject2d(x.head<2>()) / x[2];
        }

        Vector3d unproject2d(const Vector2d &v) {
            Vector3d res;
            res(0) = v(0);
            res(1) = v(1);
            res(2) = 1;
            return res;
        }
    };
}

#endif //ORB_SLAM2_EDGEPROJECTPSI2UV_H
