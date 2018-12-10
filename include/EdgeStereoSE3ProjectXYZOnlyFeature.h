//
// Created by mnowicki on 10.12.18.
//

#ifndef ORB_SLAM2_EDGESTEREOSE3PROJECTXYZONLYFEATURE_H
#define ORB_SLAM2_EDGESTEREOSE3PROJECTXYZONLYFEATURE_H


#include "Thirdparty/g2o/g2o/core/base_unary_edge.h"
#include "Thirdparty/g2o/g2o/types/se3quat.h"
#include "vertex_pointxyz.h"

namespace g2o {

    class EdgeStereoSE3ProjectXYZOnlyFeature : public BaseUnaryEdge<3, Vector3d, VertexPointXYZ> {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        EdgeStereoSE3ProjectXYZOnlyFeature() {}

        bool read(std::istream &is) {
            for (int i=0; i<3; i++){
                is >> _measurement[i];
            }
            for (int i=0; i<3; i++)
                for (int j=i; j<3; j++) {
                    is >> information()(i,j);
                    if (i!=j)
                        information()(j,i)=information()(i,j);
                }
            return true;
        }

        bool write(std::ostream &os) const {
            for (int i=0; i<3; i++){
                os << measurement()[i] << " ";
            }

            for (int i=0; i<3; i++)
                for (int j=i; j<3; j++){
                    os << " " <<  information()(i,j);
                }
            return os.good();
        }

        void computeError() {
            const VertexPointXYZ *v1 = static_cast<const VertexPointXYZ *>(_vertices[0]);
            Vector3d obs(_measurement);
            // TODO!
            _error = obs - cam_project(kfPose.map(v1->estimate()), bf);
        }

        bool isDepthPositive() {
            const VertexPointXYZ *v1 = static_cast<const VertexPointXYZ *>(_vertices[0]);
            // TODO!
            return (kfPose.map(v1->estimate()))(2) > 0.0;
        }


        virtual void linearizeOplus() {
            VertexPointXYZ *mpPosition = static_cast<g2o::VertexPointXYZ *>(_vertices[0]);
            Vector3d xyz_trans = kfPose.map(mpPosition->estimate());

            const Matrix3d R =  kfPose.rotation().toRotationMatrix();

            double x = xyz_trans[0];
            double y = xyz_trans[1];
            double z = xyz_trans[2];
            double z_2 = z*z;

            _jacobianOplusXi(0,0) = -fx*R(0,0)/z+fx*x*R(2,0)/z_2;
            _jacobianOplusXi(0,1) = -fx*R(0,1)/z+fx*x*R(2,1)/z_2;
            _jacobianOplusXi(0,2) = -fx*R(0,2)/z+fx*x*R(2,2)/z_2;

            _jacobianOplusXi(1,0) = -fy*R(1,0)/z+fy*y*R(2,0)/z_2;
            _jacobianOplusXi(1,1) = -fy*R(1,1)/z+fy*y*R(2,1)/z_2;
            _jacobianOplusXi(1,2) = -fy*R(1,2)/z+fy*y*R(2,2)/z_2;

            _jacobianOplusXi(2,0) = _jacobianOplusXi(0,0)-bf*R(2,0)/z_2;
            _jacobianOplusXi(2,1) = _jacobianOplusXi(0,1)-bf*R(2,1)/z_2;
            _jacobianOplusXi(2,2) = _jacobianOplusXi(0,2)-bf*R(2,2)/z_2;
        }



        Vector3d cam_project(const Vector3d &trans_xyz, const float &bf) {
            const float invz = 1.0f/trans_xyz[2];
            Vector3d res;
            res[0] = trans_xyz[0]*invz*fx + cx;
            res[1] = trans_xyz[1]*invz*fy + cy;
            res[2] = res[0] - bf*invz;
            return res;
        }

        g2o::SE3Quat kfPose;
        double fx, fy, cx, cy, bf;
    };
}


#endif //ORB_SLAM2_EDGESTEREOSE3PROJECTXYZONLYFEATURE_H
