//
// Author: Michal Nowicki
//

#ifndef ORB_SLAM2_TYPES_SIX_DOF_PHOTO_H
#define ORB_SLAM2_TYPES_SIX_DOF_PHOTO_H

#include "Thirdparty/g2o/g2o/core/base_vertex.h"
#include "Thirdparty/g2o/g2o/core/base_binary_edge.h"
#include "Thirdparty/g2o/g2o/core/base_unary_edge.h"
#include "Thirdparty/g2o/g2o/core/base_multi_edge.h"
#include "Thirdparty/g2o/g2o/types/vertexSE3ExpmapBright.h"
#include "Thirdparty/g2o/g2o/types/types_six_dof_expmap.h"
#include "Thirdparty/g2o/g2o/types/se3_ops.h"
#include "Thirdparty/g2o/g2o/types/se3quat.h"
#include "Thirdparty/g2o/g2o/types/types_sba.h"
#include <Eigen/Geometry>
#include <Eigen/Core>

#include "photometricErrorFunctions.h"

namespace g2o {
    using namespace std;

    typedef Eigen::Matrix<double,13,1,Eigen::ColMajor> Vector13D;
    typedef Eigen::Matrix<double, 8, 1,Eigen::ColMajor> Vector8D;

    struct SE3QuatAB {
        SE3Quat se3quat;
        double aL, bL;
    };

    class VertexSE3ExpmapAB : public BaseVertex<8, SE3QuatAB> {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        VertexSE3ExpmapAB() : BaseVertex<8, SE3QuatAB>() {

        };

        bool read(std::istream &is) {
            Vector7d est;
            for (int i = 0; i < 7; i++)
                is >> est[i];
            SE3Quat cam2world;
            cam2world.fromVector(est);


            SE3QuatAB fullEstimate;
            fullEstimate.se3quat = cam2world.inverse();
            is >> fullEstimate.aL;
            is >> fullEstimate.bL;

            setEstimate(fullEstimate);
            return true;
        }

        bool write(std::ostream &os) const {
            SE3QuatAB est = estimate();

            SE3Quat cam2world(est.se3quat.inverse());
            for (int i = 0; i < 7; i++)
                os << cam2world[i] << " ";
            os << est.aL << " " << est.bL << " ";
            return os.good();
        }

        virtual void setToOriginImpl() {
            SE3QuatAB est;
            est.se3quat = SE3Quat();
            est.aL = 0;
            est.bL = 0;
            _estimate = est;
        }

        virtual void oplusImpl(const double *update_) {
            Eigen::Map<const Vector8D> update(update_);
            Eigen::Map<const Vector6d> updatePose(update.head(6).data());

            SE3QuatAB est = estimate();
            est.se3quat = SE3Quat::exp(updatePose) * est.se3quat;
            est.aL = est.aL + update[6];
            est.bL = est.bL + update[7];
            setEstimate(est);
        }
    };

    class  EdgeSE3ProjectXYZOnlyPoseAB: public  BaseUnaryEdge<2, Vector2d, VertexSE3ExpmapAB>{
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        EdgeSE3ProjectXYZOnlyPoseAB(){}

        bool read(std::istream& is) {
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

        bool write(std::ostream& os) const {
            for (int i=0; i<2; i++){
                os << measurement()[i] << " ";
            }

            for (int i=0; i<2; i++)
                for (int j=i; j<2; j++){
                    os << " " <<  information()(i,j);
                }
            return os.good();
        }

        void computeError()  {
            const VertexSE3ExpmapAB* v1 = static_cast<const VertexSE3ExpmapAB*>(_vertices[0]);
            Vector2d obs(_measurement);
            SE3QuatAB est = v1->estimate();
            _error = obs-cam_project(est.se3quat.map(Xw));
        }

        bool isDepthPositive() {
            const VertexSE3ExpmapAB* v1 = static_cast<const VertexSE3ExpmapAB*>(_vertices[0]);
            return (v1->estimate().se3quat.map(Xw))(2)>0.0;
        }


//        virtual void linearizeOplus();

        Vector2d project2d(const Vector3d& v)  {
            Vector2d res;
            res(0) = v(0)/v(2);
            res(1) = v(1)/v(2);
            return res;
        }

        Vector2d cam_project(const Vector3d & trans_xyz) {
            Vector2d proj = project2d(trans_xyz);
            Vector2d res;
            res[0] = proj[0]*fx + cx;
            res[1] = proj[1]*fy + cy;
            return res;
        }

        Vector3d Xw;
        double fx, fy, cx, cy;
    };

    class EdgeSE3PhotoOnlyPoseAB: public BaseUnaryEdge<13, Vector13D, VertexSE3ExpmapAB> {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        EdgeSE3PhotoOnlyPoseAB()  {

            //   x
            //  xxx
            // xxxxx - neighbours used in photo check
            //  xxx
            //   x
            neighbours.push_back(make_pair(0, -2));

            neighbours.push_back(make_pair(-1, -1));
            neighbours.push_back(make_pair(0, -1));
            neighbours.push_back(make_pair(1, -1));

            neighbours.push_back(make_pair(-2, 0));
            neighbours.push_back(make_pair(-1, 0));
            neighbours.push_back(make_pair(0, 0));
            neighbours.push_back(make_pair(1, 0));
            neighbours.push_back(make_pair(2, 0));

            neighbours.push_back(make_pair(-1, 1));
            neighbours.push_back(make_pair(0, 1));
            neighbours.push_back(make_pair(1, 1));

            neighbours.push_back(make_pair(0, 2));
        }

        virtual bool read  (std::istream& is);
        virtual bool write (std::ostream& os) const;
        void computeError();
//        virtual void linearizeOplus ();

//        inline Eigen::Matrix<double, 1, 2> d_inten_d_proj(const double u, const double v) ;
//        inline Matrix<double, 2, 3, Eigen::ColMajor> d_proj_d_y(const double &fx, const double &fy, const Vector3D &xyz, const double &baseline);
//        inline Matrix<double, 3, 6, Eigen::ColMajor> d_expy_d_y(const Vector3D &y);
//        inline Matrix<double, 3, 1, Eigen::ColMajor> d_Tinvpsi_d_psi(const SE3Quat &T, const Vector3D &psi);

//        bool isDepthPositive();


        // Variables to set!
        int pyramidIndex;
        std::vector< photo::imgStr *> imgAnchor;
        std::vector< photo::imgStr *> imgObs;
        Eigen::Vector3d featureInWorld;
        Eigen::Matrix4d poseA;
        double lastU, lastV, currentU, currentV;
        double fx, fy, cx, cy;

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


        std::vector< std::pair<double, double> > neighbours;
    };

}
#endif //ORB_SLAM2_TYPES_SIX_DOF_PHOTO_H
