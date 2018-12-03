//
// Created by mnowicki on 02.12.18.
//

#ifndef ORB_SLAM2_EDGESE3PROJECTXYZONLYFEATURE_H
#define ORB_SLAM2_EDGESE3PROJECTXYZONLYFEATURE_H

#include "Thirdparty/g2o/g2o/core/base_unary_edge.h"
#include "Thirdparty/g2o/g2o/types/se3quat.h"
#include "vertex_pointxyz.h"

namespace g2o {

    class EdgeSE3ProjectXYZOnlyFeature : public BaseUnaryEdge<2, Vector2d, VertexPointXYZ> {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        EdgeSE3ProjectXYZOnlyFeature() {}

        bool read(std::istream &is) {
            for (int i=0; i<2; i++){
                is >> _measurement[i];
            }
            for (int i=0; i<2; i++)
                for (int j=i; j<2; j++) {
                    is >> information()(i,j);
                    if (i!=j)
                        information()(j,i)=information()(i,j);
                }
            return true;
        }

        bool write(std::ostream &os) const {
            for (int i=0; i<2; i++){
                os << measurement()[i] << " ";
            }

            for (int i=0; i<2; i++)
                for (int j=i; j<2; j++){
                    os << " " <<  information()(i,j);
                }
            return os.good();
        }

        void computeError() {
            const VertexPointXYZ *v1 = static_cast<const VertexPointXYZ *>(_vertices[0]);
            Vector2d obs(_measurement);
            // TODO!
            _error = obs - cam_project(kfPose.map(v1->estimate()));
        }

        bool isDepthPositive() {
            const VertexPointXYZ *v1 = static_cast<const VertexPointXYZ *>(_vertices[0]);
            // TODO!
            return (kfPose.map(v1->estimate()))(2) > 0.0;
        }


        virtual void linearizeOplus() {
            VertexPointXYZ *mpPosition = static_cast<g2o::VertexPointXYZ *>(_vertices[0]);
            Vector3d xyz_trans = kfPose.map(mpPosition->estimate());

            double x = xyz_trans[0];
            double y = xyz_trans[1];
            double z = xyz_trans[2];
            double z_2 = z*z;

            Matrix<double,2,3> tmp;
            tmp(0,0) = fx;
            tmp(0,1) = 0;
            tmp(0,2) = -x/z*fx;

            tmp(1,0) = 0;
            tmp(1,1) = fy;
            tmp(1,2) = -y/z*fy;

            _jacobianOplusXi =  -1./z * tmp * kfPose.rotation().toRotationMatrix();
        }


        Vector2d project2d(const Vector3d& v)  {
            Vector2d res;
            res(0) = v(0)/v(2);
            res(1) = v(1)/v(2);
            return res;
        }

        Vector2d cam_project(const Vector3d &trans_xyz) {
            Vector2d proj = project2d(trans_xyz);
            Vector2d res;
            res[0] = proj[0]*fx + cx;
            res[1] = proj[1]*fy + cy;
            return res;
        }

        g2o::SE3Quat kfPose;
        double fx, fy, cx, cy;
    };
}

#endif //ORB_SLAM2_EDGESE3PROJECTXYZONLYFEATURE_H
