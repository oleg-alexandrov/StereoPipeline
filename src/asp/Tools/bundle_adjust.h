// __BEGIN_LICENSE__
//  Copyright (c) 2009-2013, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NGT platform is licensed under the Apache License, Version 2.0 (the
//  "License"); you may not use this file except in compliance with the
//  License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
// __END_LICENSE__


#ifndef __ASP_TOOLS_BUNDLEADJUST_H__
#define __ASP_TOOLS_BUNDLEADJUST_H__

#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>

#include <vw/Camera/CAHVORModel.h>
#include <vw/Camera/PinholeModel.h>
#include <vw/BundleAdjustment.h>
#include <vw/Math.h>

#include <stdlib.h>
#include <iostream>

#include <asp/Sessions.h>
#include <asp/Core/BundleAdjustUtils.h>
#include <asp/Core/StereoSettings.h>

namespace asp{
  template<class camera_vector_t>
  void get_cam_params(camera_vector_t const& a_j,
                      vw::Vector3 &position_correction,
                      vw::Quat &pose_correction) {
    position_correction = subvector(a_j, 0, 3);
    pose_correction = axis_angle_to_quaternion( subvector(a_j, 3, 3) );
    //std::cout << "1a quat is " << pose_correction << std::endl;
    //std::cout << "1a axis " << subvector(a_j, 3, 3) << std::endl;
    //std::cout << "1b axis " << pose_correction.axis_angle() << std::endl;
    //std::cout << "1b quat is " << axis_angle_to_quaternion(pose_correction.axis_angle()) << std::endl;
    //std::cout << "---1 quat is " << pose_correction << std::endl;
  }

  template<class camera_vector_t>
  void set_cam_params(camera_vector_t & a_j,
                      vw::Vector3 const& position_correction,
                      vw::Quat const& pose_correction) {
    subvector(a_j, 0, 3) = position_correction;
    subvector(a_j, 3, 3) = pose_correction.axis_angle();
    //std::cout << "---2 quat is " << pose_correction << std::endl;

    //std::cout << "2a quat is " << pose_correction << std::endl;
    //std::cout << "2a axis " << subvector(a_j, 3, 3) << std::endl;
    //std::cout << "2b axis " << pose_correction.axis_angle() << std::endl;
    //std::cout << "2b quat is " << axis_angle_to_quaternion(pose_correction.axis_angle()) << std::endl;

  }

  template<class ModelT>
  void concat_extrinsics_intrinsics(const double* const extrinsics,
                                    const double* const intrinsics,
                                    typename ModelT::camera_intr_vector_t & concat){
    int intr_len = ModelT::intrinsic_params_n, cam_len = ModelT::camera_params_n;
    for (int c = 0; c < cam_len; c++)
      concat[c] = extrinsics[c];
    for (int i = 0; i < intr_len; i++)
      concat[cam_len + i] = intrinsics[i];
  }

}

// Bundle adjustment functor
class BundleAdjustmentModel:
  public vw::ba::ModelBase<BundleAdjustmentModel, 6, 3> {

public:
  typedef vw::Vector<double,camera_params_n> camera_vector_t;
  typedef vw::Vector<double,point_params_n> point_vector_t;
  typedef boost::shared_ptr<vw::camera::CameraModel> cam_ptr_t;
  // No intrinsic params
  const static size_t intrinsic_params_n = 0;
  typedef vw::Vector<double,intrinsic_params_n> intrinsic_vector_t;
  typedef vw::Vector<double,camera_params_n
                     +intrinsic_params_n> camera_intr_vector_t;
private:
  std::vector<cam_ptr_t> m_cameras;
  boost::shared_ptr<vw::ba::ControlNetwork> m_network;

  std::vector<camera_intr_vector_t> a;
  std::vector<point_vector_t> b;
  std::vector<camera_vector_t> a_target;
  std::vector<point_vector_t> b_target;
  int m_num_pixel_observations;

public:
  BundleAdjustmentModel(std::vector<cam_ptr_t> const& cameras,
                        boost::shared_ptr<vw::ba::ControlNetwork> network) :
    m_cameras(cameras), m_network(network), a(cameras.size()),
    b(network->size()), a_target(cameras.size()), b_target(network->size()) {

    // Compute the number of observations from the bundle.
    m_num_pixel_observations = 0;
    for (unsigned i = 0; i < network->size(); ++i)
      m_num_pixel_observations += (*network)[i].size();

    // Set up the 'b' vectors, storing the initial values.
    // The 'a' vectors however just start out zero.
    for (unsigned i = 0; i < network->size(); ++i) {
      b[i] = (*m_network)[i].position();
      b_target[i] = b[i];
    }
  }

  // Return a reference to the camera and point parameters.
  camera_vector_t A_parameters(int j) const { return a[j]; }
  point_vector_t B_parameters(int i) const { return b[i]; }
  void set_A_parameters(int j, camera_intr_vector_t const& a_j) {
    a[j] = a_j;
  }
  void set_B_parameters(int i, point_vector_t const& b_i) {
    b[i] = b_i;
  }

  // Return the initial parameters
  camera_vector_t A_target(int j) const { return a_target[j]; }
  point_vector_t B_target(int i) const { return b_target[i]; }

  unsigned num_cameras() const { return a.size(); }
  unsigned num_points() const { return b.size(); }
  unsigned num_pixel_observations() const { return m_num_pixel_observations; }

  // Return the covariance of the camera parameters for camera j.
  inline vw::Matrix<double,camera_params_n,camera_params_n>
  A_inverse_covariance ( unsigned /*j*/ ) const {
    vw::Matrix<double,camera_params_n,camera_params_n> result;
    result(0,0) = 1/100.0;
    result(1,1) = 1/100.0;
    result(2,2) = 1/100.0;
    result(3,3) = 1/1e-1;
    result(4,4) = 1/1e-1;
    result(5,5) = 1/1e-1;
    return result;
  }

  // Return the covariance of the point parameters for point i.
  inline vw::Matrix<double,point_params_n,point_params_n>
  B_inverse_covariance ( unsigned /*i*/ ) const {
    vw::Matrix<double,point_params_n,point_params_n> result;
    result(0,0) = 1/20;
    result(1,1) = 1/20;
    result(2,2) = 1/20;
    return result;
  }

  // Given the 'a' vector (camera model parameters) for the j'th
  // image, and the 'b' vector (3D point location) for the i'th
  // point, return the location of b_i on imager j in pixel
  // coordinates.
  vw::Vector2 operator() ( unsigned /*i*/, unsigned j,
                           camera_vector_t const& a_j,
                           point_vector_t const& b_i) const {
    vw::Vector3 position_correction;
    vw::Quat pose_correction;
    asp::get_cam_params(a_j, position_correction, pose_correction);
    vw::camera::AdjustedCameraModel cam(m_cameras[j],
                                        position_correction,
                                        pose_correction);
    return cam.point_to_pixel(b_i);
  }

  void write_adjustment(int j, std::string const& filename) const {
    vw::Vector3 position_correction;
    vw::Quat pose_correction;
    asp::get_cam_params(a[j], position_correction, pose_correction);
    write_adjustments(filename, position_correction, pose_correction);
  }

  std::vector<cam_ptr_t> adjusted_cameras() const {
    std::vector<cam_ptr_t> result(m_cameras.size());
    for (unsigned j = 0; j < result.size(); ++j) {
      vw::Vector3 position_correction;
      vw::Quat pose_correction;
      asp::get_cam_params(a[j], position_correction, pose_correction);
      result[j] = cam_ptr_t( new vw::camera::AdjustedCameraModel( m_cameras[j], position_correction, pose_correction ) );
    }
    return result;
  }

  inline double image_compare( vw::Vector2 const& meas,
                               vw::Vector2 const& obj ) {
    return norm_2( meas - obj );
  }

  inline double position_compare( camera_vector_t const& meas,
                                  camera_vector_t const& obj ) {
    return norm_2( subvector(meas,0,3) - subvector(obj,0,3) );
  }

  inline double pose_compare( camera_vector_t const& meas,
                              camera_vector_t const& obj ) {
    return norm_2( subvector(meas,3,3) - subvector(obj,3,3) );
  }

  inline double gcp_compare( point_vector_t const& meas,
                             point_vector_t const& obj ) {
    return norm_2(meas - obj);
  }

  // Give access to the control network
  boost::shared_ptr<vw::ba::ControlNetwork> control_network() const {
    return m_network;
  }

  void bundlevis_cameras_append(std::string const& filename) const {
    std::ofstream ostr(filename.c_str(),std::ios::app);
    for ( unsigned j = 0; j < a.size(); j++ ) {
      vw::Vector3 position_correction;
      vw::Quat pose_correction;
      asp::get_cam_params(a[j], position_correction, pose_correction);
      vw::camera::AdjustedCameraModel cam(m_cameras[j],
                                          position_correction,
                                          pose_correction);
      vw::Vector3 position = cam.camera_center( vw::Vector2() );
      vw::Quat pose = cam.camera_pose( vw::Vector2() );
      ostr << std::setprecision(18) << j << "\t" << position[0] << "\t"
           << position[1] << "\t" << position[2] << "\t";
      ostr << pose[0] << "\t" << pose[1] << "\t"
           << pose[2] << "\t" << pose[3] << "\n";
    }
  }

  void bundlevis_points_append(std::string const& filename) const {
    std::ofstream ostr(filename.c_str(),std::ios::app);
    unsigned i = 0;
    BOOST_FOREACH( point_vector_t const& p, b ) {
      ostr << i++ << std::setprecision(18) << "\t" << p[0] << "\t"
           << p[1] << "\t" << p[2] << "\n";
    }
  }
};

// Model to be used to float all parameters of a pinhole model.  There
// are 6 camera parameters, corresponding to: camera center (3), and
// camera orientation (3). Also there are three intrinsic parameters:
// focal length (1), and pixel offsets (2), which are shared
// among the cameras.
class BAPinholeModel:
  public vw::ba::ModelBase<BAPinholeModel, 6, 3> {

public:
  typedef vw::Vector<double,camera_params_n> camera_vector_t;
  typedef vw::Vector<double,point_params_n> point_vector_t;
  typedef boost::shared_ptr<vw::camera::CameraModel> cam_ptr_t;
  // Three intrinsic params
  const static size_t intrinsic_params_n = 3;
  typedef vw::Vector<double,intrinsic_params_n> intrinsic_vector_t;
  typedef vw::Vector<double,camera_params_n
                     +intrinsic_params_n> camera_intr_vector_t;
  // Need this scale to force the rotations to not change that wildly
  // when determining pinhole cameras from scratch.
  double pose_scale; //  = 1.0e+6;
private:
  std::vector<cam_ptr_t> m_cameras;
  boost::shared_ptr<vw::ba::ControlNetwork> m_network;
  std::vector<camera_intr_vector_t> a;

public:
  BAPinholeModel(std::vector<cam_ptr_t> const& cameras,
                 boost::shared_ptr<vw::ba::ControlNetwork> network) :
    m_cameras(cameras), m_network(network), a(cameras.size()){
    pose_scale = atof(getenv("P"));
    std::cout << "---pose scale is " << pose_scale << std::endl;
    if (pose_scale <= 0) {
      std::cout << "---error in pose_scale!" << std::endl;
      exit(0);
    }
  }

  unsigned num_cameras() const { return m_cameras.size(); }
  unsigned num_points() const { return m_network->size(); }

  void set_A_parameters(int j, camera_intr_vector_t const& a_j) {
    // Set the camera parameters.
    // Get rid of this function!!!!
    a[j] = a_j;
  }

  // Get pinhole model from its internal representation of arrays
  // used for optimization.
  vw::camera::PinholeModel from_internal(camera_intr_vector_t const& a_j){

    camera_intr_vector_t b = a_j;

    // Undo the scale of the rotation variables
    subvector(b,3,3) /= pose_scale;

    vw::Vector3 position;
    vw::Quat pose;
    asp::get_cam_params(b, position, pose);

    return vw::camera::PinholeModel(position,
                                    pose.rotation_matrix(),
                                    b[6], b[6],  // focal lengths
                                    b[7], b[8]); // pixel offsets
  }

  // Copy pinhole models into arrays used for optimization.
  // This reverses the process in from_internal().

  void to_internal(vw::camera::PinholeModel & model,
                              camera_intr_vector_t & a_j){

    vw::Vector3 position = model.camera_center();
    vw::Quat pose = model.camera_pose();

    // Set the extrinsics
    asp::set_cam_params(a_j, position, pose);
    // Scale the rotation variables
    subvector(a_j,3,3) *= pose_scale;

    // Set the intrinsic_params_n
    vw::Vector2 fl = model.focal_length();
    vw::Vector2 po = model.point_offset();
    a_j[6] = fl[0];
    a_j[7] = po[0];
    a_j[8] = po[1];
  }

  void from_internal(std::vector<boost::shared_ptr<vw::camera::CameraModel> > & cameras){
    VW_ASSERT(cameras.size() == a.size(),
              vw::ArgumentErr() << "Book-keeping error in number of cameras.\n");
    int num_cams = a.size();
    for (int icam = 0; icam < num_cams; icam++){
      vw::camera::PinholeModel * pincam
        = dynamic_cast<vw::camera::PinholeModel*>(cameras[icam].get());
      VW_ASSERT(pincam != NULL,
                vw::ArgumentErr() << "A pinhole camera expected.\n");

      *pincam = from_internal(a[icam]);
      //std::cout << "---input val "  << ' ' << icam  << ' ' << a[icam] << std::endl;
    }
  }

  void to_internal(std::vector<boost::shared_ptr<vw::camera::CameraModel> > const& cameras){
    VW_ASSERT(cameras.size() == a.size(),
              vw::ArgumentErr() << "Book-keeping error in number of cameras.\n");
    int num_cams = a.size();
    for (int icam = 0; icam < num_cams; icam++){
      vw::camera::PinholeModel * pincam
        = dynamic_cast<vw::camera::PinholeModel*>(cameras[icam].get());
      VW_ASSERT(pincam != NULL,
                vw::ArgumentErr() << "A pinhole camera expected.\n");

      to_internal(*pincam, a[icam]);
    }
  }

  // Given the 'a' vector (camera model parameters) for the j'th
  // image, and the 'b' vector (3D point location) for the i'th
  // point, return the location of b_i on imager j in pixel
  // coordinates.
  vw::Vector2 operator() ( unsigned /*i*/, unsigned j,
                           camera_intr_vector_t const& a_j,
                           point_vector_t const& b_i) {

    return from_internal(a_j).point_to_pixel(b_i);
  }

  // Give access to the control network
  boost::shared_ptr<vw::ba::ControlNetwork> control_network() const {
    return m_network;
  }

  void write_camera_models(std::vector<std::string> const& cam_files){

    VW_ASSERT(cam_files.size() == a.size(),
              vw::ArgumentErr() << "Must have as many camera files as cameras.\n");

    for (int icam = 0; icam < (int)cam_files.size(); icam++){
      vw::vw_out() << "Writing: " << cam_files[icam] << std::endl;
      from_internal(a[icam]).write(cam_files[icam]);
    }

  }


};

#endif
