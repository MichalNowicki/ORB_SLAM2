////
//// Author: Michal Nowicki
////
//
//#ifndef ORB_SLAM2_TYPES_SIX_DOF_PHOTO_H
//#define ORB_SLAM2_TYPES_SIX_DOF_PHOTO_H
//
//#include "../core/base_vertex.h"
//#include "../core/base_binary_edge.h"
//#include "../core/base_unary_edge.h"
//#include "../core/base_multi_edge.h"
//#include "../types/vertexSE3ExpmapBright.h"
//#include "../types/types_six_dof_expmap.h"
//#include "se3_ops.h"
//#include "se3quat.h"
//#include "types_sba.h"
//#include <Eigen/Geometry>
//#include <Eigen/Core>
//
//namespace g2o {
//    using namespace std;
//
//    typedef Eigen::Matrix<double,9,1,Eigen::ColMajor> Vector9D;
//
//    struct imgStr {
//        float imageScale;
//        std::vector< std::vector< float> > image;
//        std::vector< std::vector< Eigen::Vector2f > > gradient;
//    };
//
//    class EdgeSE3PhotoOnlyPose: public BaseUnaryEdge<9, Vector9D, VertexSE3Expmap> {
//    public:
//        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
//
//        EdgeSE3PhotoOnlyPose()  {
//            //   x
//            //  x x
//            // x x x - neighbours used in optimization
//            //  x x
//            //   x
//            neighbours.push_back(make_pair(0,0));
//            neighbours.push_back(make_pair(0,2));
//            neighbours.push_back(make_pair(1,1));
//            neighbours.push_back(make_pair(2,0));
//            neighbours.push_back(make_pair(1,-1));
//            neighbours.push_back(make_pair(0,-2));
//            neighbours.push_back(make_pair(-1,-1));
//            neighbours.push_back(make_pair(-2,0));
//            neighbours.push_back(make_pair(-1,1));
//        }
//
//        void setAdditionalData(std::vector< imgStr *> &imageAnchor,
//                               std::vector< imgStr *> &imageObs,
//                               double _baseline) {
//            imgAnchor = imageAnchor;
//            imgObs = imageObs;
//        }
//
//        void selectPyramidIndex(int _pyramidIndex) {
//            pyramidIndex = _pyramidIndex;
//        }
//
//
//        virtual bool read  (std::istream& is);
//        virtual bool write (std::ostream& os) const;
//        void computeError();
////        virtual void linearizeOplus ();
//
////        inline Eigen::Matrix<double, 1, 2> d_inten_d_proj(const double u, const double v) ;
////        inline Matrix<double, 2, 3, Eigen::ColMajor> d_proj_d_y(const double &fx, const double &fy, const Vector3D &xyz, const double &baseline);
////        inline Matrix<double, 3, 6, Eigen::ColMajor> d_expy_d_y(const Vector3D &y);
////        inline Matrix<double, 3, 1, Eigen::ColMajor> d_Tinvpsi_d_psi(const SE3Quat &T, const Vector3D &psi);
//
////        bool isDepthPositive();
//
//        // Variables to set!
//        int pyramidIndex;
//        std::vector< imgStr *> imgAnchor;
//        std::vector< imgStr *> imgObs;
//        Eigen::Vector3d featureInWorld;
//        Eigen::Matrix4d poseA;
//        double lastU, lastV, currentU, currentV;
//        double fx, fy, cx, cy;
//
//    private:
//        inline Vector3D invert_depth(const Vector3D &x) {
//            return unproject2d(x.head<2>()) / x[2];
//        }
//
//        Vector3d unproject2d(const Vector2d &v) {
//            Vector3d res;
//            res(0) = v(0);
//            res(1) = v(1);
//            res(2) = 1;
//            return res;
//        }
//
//        inline double getSubpixImageValue(double u, double v, std::vector< std::vector< float> > &image);
//
//        double getDistanceToPlane(const Eigen::Vector3d &point3D, const Eigen::Vector3d &normal) {
//            return -normal.transpose() * point3D;
//        }
//
//        Eigen::Matrix3d computeHomography(Eigen::Matrix4d Tba, Eigen::Vector3d n, double d, Eigen::Matrix3d Ka, Eigen::Matrix3d Kb) {
//            // Getting R,t
//            Eigen::Matrix3d R21 = Tba.block<3,3>(0,0);
//            Eigen::Vector3d t21 = Tba.block<3,1>(3,0);
//
//            Eigen::Matrix3d H = Ka * (R21 - t21 * n.transpose() / d) * Kb.inverse();
//
//            // Homography
//            return H;
//        }
//
//        Eigen::Matrix3d getCameraMatrix(float fx, float fy, float cx, float cy) {
//            Eigen::Matrix3d cameraMatrix = Eigen::Matrix3d::Identity();
//            cameraMatrix(0,0) = fx;
//            cameraMatrix(0,2) = cx;
//            cameraMatrix(1,1) = fy;
//            cameraMatrix(1,2) = cy;
//            cameraMatrix(1,2) = cy;
//
//            return cameraMatrix;
//        }
//
//        std::vector< std::pair<double, double> > neighbours;
//    };
//
//}
//#endif //ORB_SLAM2_TYPES_SIX_DOF_PHOTO_H
