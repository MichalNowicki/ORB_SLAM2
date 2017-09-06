//
// Created by michalnowicki on 06.09.17.
//

#ifndef ORB_SLAM2_EDGEPROJECTPSI2UVSINGLEPARAMPATCHBRIGHT_H
#define ORB_SLAM2_EDGEPROJECTPSI2UVSINGLEPARAMPATCHBRIGHT_H

#include "../core/base_vertex.h"
#include "../core/base_binary_edge.h"
#include "../core/base_unary_edge.h"
#include "../core/base_multi_edge.h"
#include "../types/types_six_dof_expmap.h"
#include "se3_ops.h"
#include "se3quat.h"
#include "types_sba.h"
#include <Eigen/Geometry>
#include <Eigen/Core>
#include "../types/VertexSE3ExpmapBright.h"

namespace g2o {
    using namespace std;

    typedef Eigen::Matrix<double,9,1,Eigen::ColMajor> Vector9D;

    class EdgeProjectPSI2UVSingleParamPatchBright : public g2o::BaseMultiEdge<9, Vector9D> {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        EdgeProjectPSI2UVSingleParamPatchBright()  {
            resizeParameters(1);
            installParameter(_cam, 0);

            //   x
            //  x x
            // x x x - neighbours used in optimization
            //  x x
            //   x
            neighbours.push_back(make_pair(0,0));
            neighbours.push_back(make_pair(0,2));
            neighbours.push_back(make_pair(1,1));
            neighbours.push_back(make_pair(2,0));
            neighbours.push_back(make_pair(1,-1));
            neighbours.push_back(make_pair(0,-2));
            neighbours.push_back(make_pair(-1,-1));
            neighbours.push_back(make_pair(-2,0));
            neighbours.push_back(make_pair(-1,1));
        }

        void setAdditionalData(std::vector<double> _largePatchAnchor, std::vector<double> _largePatchObs,
                               double _pointAnchorScale, double _pointObsScale) {
            largePatchAnchor = _largePatchAnchor;
            largePatchObs = _largePatchObs;
            pointAnchorScale = _pointAnchorScale;
            pointObsScale = _pointObsScale;

            int largePatchElements = largePatchAnchor.size();

            largePatchStride = sqrt(largePatchElements);
            largePatchCenter = (largePatchStride - 1) / 2;


            for (int i=0;i<largePatchElements; i++)
                largePatchObsGradient.push_back(Vector2D(0,0));

            // Selecting every non-border element
            for (int i = 1, index = largePatchStride + 1; i<largePatchStride-1;i++ )
            {
                for (int j = 1;j<largePatchStride-1;j++) {

                    // Gradients
                    Vector2D imageGradient;
                    imageGradient[0] = 0.5 * (largePatchObs[i * largePatchStride  + j+1] - largePatchObs[i*largePatchStride +j-1]);
                    imageGradient[1] = 0.5 * (largePatchObs[(i+1)*largePatchStride +j] - largePatchObs[(i-1)*largePatchStride+j]);
                    largePatchObsGradient[index] = imageGradient;
                    index++;
                }
                index = index + 2; // The first and last value in each row is equal to 0
            }

//            std::cout << "Edge::setAdditionalData: origScale: " << pointAnchorScale << " pointObsScale: " << pointObsScale << std::endl;

//            std::cout << "PATCH" << std::endl;
//            for (int i=0;i<largePatchElements;i++) {
//                std :: cout << largePatchObs[i] << " ";
//                if (i%largePatchStride == largePatchStride - 1)
//                    std::cout << std::endl;
//            }
//
//            std::cout << "GRADIENT x" << std::endl;
//            for (int i=0;i<largePatchElements;i++) {
//                std :: cout << largePatchObsGradient[i][0] << " ";
//                if (i%largePatchStride == largePatchStride - 1)
//                    std::cout << std::endl;
//            }
//
//            std::cout << "GRADIENT Y" << std::endl;
//            for (int i=0;i<largePatchElements;i++) {
//                std :: cout << largePatchObsGradient[i][1] << " ";
//                if (i%largePatchStride == largePatchStride - 1)
//                    std::cout << std::endl;
//            }

        }

        virtual bool read  (std::istream& is);
        virtual bool write (std::ostream& os) const;
        void computeError();
//        virtual void linearizeOplus ();


        inline Eigen::Matrix<double, 1, 2> d_inten_d_proj(const double u, const double v) ;
        inline Matrix<double, 2, 3, Eigen::ColMajor> d_proj_d_y(const double &fx, const double &fy, const Vector3D &xyz);
        inline Matrix<double, 3, 6, Eigen::ColMajor> d_expy_d_y(const Vector3D &y);
        inline Matrix<double, 3, 1, Eigen::ColMajor> d_Tinvpsi_d_psi(const SE3Quat &T, const Vector3D &psi);

        bool isDepthPositive();

        CameraParameters * _cam;

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

        std::vector<double> largePatchAnchor, largePatchObs;
        std::vector< Vector2D > largePatchObsGradient;
        double pointAnchorScale, pointObsScale;
        int largePatchCenter, largePatchStride;
    };
}

#endif //ORB_SLAM2_EDGEPROJECTPSI2UVSINGLEPARAMPATCHBRIGHT_H
