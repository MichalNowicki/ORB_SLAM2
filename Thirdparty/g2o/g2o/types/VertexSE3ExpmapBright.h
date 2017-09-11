//
// Created by michalnowicki on 06.09.17.
//

#ifndef ORB_SLAM2_VERTEXSE3EXPMAPBRIGHT_H
#define ORB_SLAM2_VERTEXSE3EXPMAPBRIGHT_H

#include "../core/base_vertex.h"
#include "../core/base_binary_edge.h"
#include "../core/base_unary_edge.h"
#include "../core/base_multi_edge.h"
#include "se3_ops.h"
#include "se3quat.h"
#include "types_sba.h"
#include <Eigen/Geometry>

namespace g2o {

    using namespace Eigen;
    typedef Matrix<double, 8, 1> Vector8d;

    struct SE3QuatBright {
        SE3Quat se3quat;
        double a, b;
    };

    class VertexSE3ExpmapBright : public BaseVertex<8, SE3QuatBright> {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        VertexSE3ExpmapBright();

        bool read(std::istream &is);

        bool write(std::ostream &os) const;

        virtual void setToOriginImpl() {
            SE3QuatBright est;
            est.se3quat = SE3Quat();
            est.a = 0;
            est.b = 0;
            _estimate = est;
        }

        virtual void oplusImpl(const double *update_) {
            Eigen::Map<const Vector8d> update(update_);
            Eigen::Map<const Vector6d> updatePose(update.head(6).data());

            SE3QuatBright est = estimate();
            est.se3quat = SE3Quat::exp(updatePose) * est.se3quat;
            est.a = est.a + update[7];
            est.b = est.b + update[8];
            setEstimate(est);
        }
    };
}
#endif //ORB_SLAM2_VERTEXSE3EXPMAPBRIGHT_H
