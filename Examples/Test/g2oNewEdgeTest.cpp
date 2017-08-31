//
// Created by michalnowicki on 25.08.17.
//
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <fstream>
#include <map>
#include <utility>


#include "Thirdparty/g2o/g2o/core/block_solver.h"
#include "Thirdparty/g2o/g2o/core/optimization_algorithm_levenberg.h"
#include "Thirdparty/g2o/g2o/solvers/linear_solver_eigen.h"
#include "Thirdparty/g2o/g2o/types/types_six_dof_expmap.h"
#include "Thirdparty/g2o/g2o/core/robust_kernel_impl.h"
#include "Thirdparty/g2o/g2o/solvers/linear_solver_dense.h"
#include "Thirdparty/g2o/g2o/types/types_seven_dof_expmap.h"
#include "Thirdparty/g2o/g2o/types/EdgeProjectPSI2UV.h"
#include "Thirdparty/g2o/g2o/types/EdgeProjectPSI2UVSingleParam.h"


#include "Converter.h"

#include <opencv2/calib3d.hpp>


using namespace std;


class ORBSLAMBA {
public:
    struct camera {
        int id;
        double r1, r2, r3, x, y, z;
    };

    struct point {
        int id;
        Eigen::Matrix<double, 3, 1> pos;
    };

    struct measurement {
        int idC, idP;
        double u, v, invSigma2;
    };

    ORBSLAMBA(const std::string &filename) {
        std::ifstream readStream(filename.c_str());
        readStream >> camNum >> pNum >> obsNum;
        std::cout << "camNum = " << camNum << " pNum = " << pNum << " obsNum = " << obsNum << std::endl;

        for (int i = 0; i < camNum; i++) {
            camera c;
            readStream >> c.id >> c.r1 >> c.r2 >> c.r3 >> c.x >> c.y >> c.z >> fx >> fy >> cx >> cy;
            cameras.push_back(c);
        }


        for (int i = 0; i < pNum; i++) {
            point p;
            double x, y, z;
            readStream >> p.id >> x >> y >> z;
            p.pos << x, y, z;

            points.push_back(p);
        }

        for (int i = 0, fid = 0; i < obsNum; i++) {
            measurement m;
            readStream >> m.idC >> m.idP;
            readStream >> m.u >> m.v >> m.invSigma2;

            if (firstMeasurement.find(m.idP) == firstMeasurement.end()) {
                firstMeasurement[m.idP] = fid++;
                firstObs.push_back(m);
            } else
                measurements.push_back(m);
        }

    }

    std::vector<camera> cameras;
    std::vector<point> points;
    std::vector<measurement> measurements;
    std::vector<measurement> firstObs;
    std::map<int, int> firstMeasurement;

    int camNum, pNum, obsNum;
    double fx, fy, cx, cy;

};


Eigen::Vector3d unproject2d(const Eigen::Vector2d& v){
    Eigen::Vector3d res;
    res(0) = v(0);
    res(1) = v(1);
    res(2) = 1;
    return res;
}

inline Eigen::Vector3d invert_depth(const Eigen::Vector3d & x){
    return unproject2d(x.head<2>())/x[2];
}


int main(int argc, char * argv[]) {
    ORBSLAMBA orbslamba(argv[1]);

    // Setup optimizer
    g2o::SparseOptimizer optimizer;
    g2o::BlockSolver_6_3::LinearSolverType *linearSolver;

    linearSolver = new g2o::LinearSolverEigen<g2o::BlockSolver_6_3::PoseMatrixType>();

    g2o::BlockSolver_6_3 *solver_ptr = new g2o::BlockSolver_6_3(linearSolver);

    g2o::OptimizationAlgorithmLevenberg *solver = new g2o::OptimizationAlgorithmLevenberg(solver_ptr);
    optimizer.setAlgorithm(solver);


    Eigen::Vector2d principal_point(orbslamba.cx, orbslamba.cy);
    g2o::CameraParameters * cam_params
            = new g2o::CameraParameters (orbslamba.fx, orbslamba.fy, principal_point, 0.);
    cam_params->setId(0);

    if (!optimizer.addParameter(cam_params)){
        assert(false);
    }

    // Creating poses
    bool first = true;
    for (auto c : orbslamba.cameras) {
        g2o::VertexSE3Expmap *vSE3 = new g2o::VertexSE3Expmap();


        // Expmap -> rotation matrix and translation vector -> Quaternion
        cv::Mat Rexp = cv::Mat(3, 1, CV_32F), R = cv::Mat(3, 3, CV_32F), T = cv::Mat(4, 4, CV_32F);
        Rexp.at<float>(0) = c.r1;
        Rexp.at<float>(1) = c.r2;
        Rexp.at<float>(2) = c.r3;

        cv::Rodrigues(Rexp, R);
        T.rowRange(0, 3).colRange(0, 3) = R;
        T.row(0).col(3) = c.x;
        T.row(1).col(3) = c.y;
        T.row(2).col(3) = c.z;

        vSE3->setEstimate(ORB_SLAM2::Converter::toSE3Quat(T));
        vSE3->setId(c.id);
        if (first) {
            vSE3->setFixed(true);
            first = false;
        }

        optimizer.addVertex(vSE3);


        std::cout <<"Added pose " << c.id << " : (x,y,z) = (" << c.x << ", " << c.y << ", " << c.z << ")" << std::endl;
    }


    // Creating points and setting its initial inverse depth
    std::map<int, g2o::VertexSBAPointXYZ *> vertexMap;
    for (auto p : orbslamba.points) {

        g2o::VertexSBAPointXYZ *vPoint = new g2o::VertexSBAPointXYZ();

        vPoint->setEstimate( invert_depth(p.pos) );

        vPoint->setId(p.id);
        vPoint->setMarginalized(true);
        optimizer.addVertex(vPoint);

        vertexMap[p.id]=vPoint;

        std::cout <<"Added point " << p.id << " : (invD) = (" << 1/ p.pos[2] << ")" << std::endl;
    }


    const float thHuberMono = sqrt(5.991);

    vector<g2o::EdgeProjectPSI2UV *> vpEdgesMono;
    vpEdgesMono.reserve(orbslamba.obsNum);

    // Adding edges that join 3 nodes - point (inv. depth), pose with first obs (pose) and current pose (pos)
    // The measurement of the edge is the (u,v) in the second pose
    for (auto m : orbslamba.measurements) {

        Eigen::Matrix<double, 2, 1> obs;
        obs << m.u, m.v;

        g2o::EdgeProjectPSI2UV *e = new g2o::EdgeProjectPSI2UV();
        e->resize(3);

        int firstId = orbslamba.firstMeasurement[m.idP];
        auto firstObs = orbslamba.firstObs[firstId];

        e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex *>(optimizer.vertex(m.idP)));
        e->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex *>(optimizer.vertex(firstObs.idC)));
        e->setVertex(2, dynamic_cast<g2o::OptimizableGraph::Vertex *>(optimizer.vertex(m.idC)));



        e->setMeasurement(obs);

        e->setInformation(Eigen::Matrix2d::Identity() * m.invSigma2);

        g2o::RobustKernelHuber* rk = new g2o::RobustKernelHuber;
        e->setRobustKernel(rk);
        rk->setDelta(thHuberMono);



        e->setParameterId(0, 0);
        optimizer.addEdge(e);
        vpEdgesMono.push_back(e);

        std::cout <<"Measurement - p: " << m.idP  << " firstObsId: " << firstObs.idC << " pose: " <<  m.idC << std::endl;
        std::cout <<"\t Original (u,v) = (" << firstObs.u << ", " << firstObs.v << ")  Measured: " << m.u << ", "<< m.v << std::endl;
        std::cout <<"\t fx = " << orbslamba.fx << " fy = " << orbslamba.fy << " cx = " << orbslamba.cx << " cy = " << orbslamba.cy << std::endl;

        // TODO: For debugging
//        e->computeError();
    }


    optimizer.initializeOptimization(0);

    optimizer.computeActiveErrors();
    double loadChi=optimizer.activeChi2();
    std::cout << "Pre average chi2: " << loadChi / orbslamba.obsNum << std::endl;


    optimizer.optimize(10);


    double avgchi2 = 0;
    for (auto e : vpEdgesMono) {
        avgchi2 += e->chi2();
    }
    std::cout << "Post average chi2: " << avgchi2 / orbslamba.obsNum << std::endl;


//    double avgReprojError = 0;
//    for (auto e : vpEdgesMono) {
//
//        const g2o::VertexSBAPointXYZ *point = static_cast<const g2o::VertexSBAPointXYZ *>(e->vertex(0));
//        const g2o::VertexSE3Expmap *firstPose = static_cast<const g2o::VertexSE3Expmap *>(e->vertex(1));
//        const g2o::VertexSE3Expmap *secondPose = static_cast<const g2o::VertexSE3Expmap *>(e->vertex(2));
//
//        Eigen::Vector3d pointInFirst;
//
//
//
//        Eigen::Vector3d pointInGlobal = firstPose->estimate().inverse().map(invert_depth(point->estimate()));
//        Eigen::Vector2d measurement = e->measurement();
//        Eigen::Vector2d projectedPoint = e->cam_project(secondPose->estimate().map(pointInGlobal));
//
//
//        std::cout << measurement[0] << " " << measurement[1] << " " << projectedPoint[0] << " " << projectedPoint[1] << std::endl;
//
//        avgReprojError += sqrt(pow(measurement[0] - projectedPoint[0],2) + pow(measurement[1] - projectedPoint[1],2));
//    }
//
//    std::cout <<"Avg reprojection error: " << avgReprojError / vpEdgesMono.size() << std::endl;
}