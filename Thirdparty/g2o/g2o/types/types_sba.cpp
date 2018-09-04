// g2o - General Graph Optimization
// Copyright (C) 2011 Kurt Konolige
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright
//   notice, this list of conditions and the following disclaimer in the
//   documentation and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
// IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
// TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
// TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "types_sba.h"
#include <iostream>

namespace g2o {

  using namespace std;


  VertexSBAPointXYZ::VertexSBAPointXYZ() : BaseVertex<3, Vector3d>()
  {
  }

  bool VertexSBAPointXYZ::read(std::istream& is)
  {
    Vector3d lv;
    for (int i=0; i<3; i++)
      is >> _estimate[i];
    return true;
  }

  bool VertexSBAPointXYZ::write(std::ostream& os) const
  {
    Vector3d lv=estimate();
    for (int i=0; i<3; i++){
      os << lv[i] << " ";
    }
    return os.good();
  }

    CameraParameters
    ::CameraParameters()
            : focal_length_x(1.),
              focal_length_y(1.),
              principle_point(Vector2D(0., 0.)),
              baseline(0.5) {
    }

    Vector2D project2dt(const Vector3D& v)  {
      Vector2D res;
      res(0) = v(0)/v(2);
      res(1) = v(1)/v(2);
      return res;
    }

    Vector2D CameraParameters::cam_map(const Vector3D &trans_xyz) const {
      Vector2D proj = project2dt(trans_xyz);
      Vector2D res;
      res[0] = proj[0] * focal_length_x + principle_point[0];
      res[1] = proj[1] * focal_length_y + principle_point[1];
      return res;
    }

    Vector2D CameraParameters::mostcam_map(const Vector3D &trans_xyz, const double base) const {

      Vector2D res;
      res[0] = (trans_xyz[0] - base) / trans_xyz[2] * focal_length_x + principle_point[0];
      res[1] = (trans_xyz[1]) / trans_xyz[2] * focal_length_y + principle_point[1];
      return res;
    }

    Vector3D CameraParameters::stereocam_uvu_map(const Vector3D &trans_xyz) const {
      Vector2D uv_left = cam_map(trans_xyz);
      double proj_x_right = (trans_xyz[0] - baseline) / trans_xyz[2];
      double u_right = proj_x_right * focal_length_x + principle_point[0];
      return Vector3D(uv_left[0], uv_left[1], u_right);
    }

    VertexSBAPointInvD::VertexSBAPointInvD() : BaseVertex<1, double>() {

    };

    bool VertexSBAPointInvD::read(std::istream &is) {
      return false;
    };

    bool VertexSBAPointInvD::write(std::ostream &os) const {
      return false;
    };

} // end namespace
