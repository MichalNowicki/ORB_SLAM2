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



    class EdgeSE3PhotoOnlyPose: public BaseUnaryEdge<13, Vector13D, VertexSE3Expmap> {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        EdgeSE3PhotoOnlyPose()  {

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


        int errorSize;

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
