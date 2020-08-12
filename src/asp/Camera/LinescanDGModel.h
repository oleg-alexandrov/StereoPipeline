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


/// \file LinescanDGModel.h
///
/// A generic linescan camera model object
///
///
#ifndef __STEREO_CAMERA_LINESCAN_DG_MODEL_H__
#define __STEREO_CAMERA_LINESCAN_DG_MODEL_H__

#include <vw/Camera/LinescanModel.h>
#include <vw/Camera/PinholeModel.h>
#include <vw/Camera/Extrinsics.h>
#include <vw/Camera/CameraSolve.h>

#include <asp/Core/StereoSettings.h>

#include <vw/Math/EulerAngles.h>
#include <vw/Camera/CameraSolve.h>
#include <asp/Core/StereoSettings.h>
#include <asp/Camera/RPCModel.h>
#include <asp/Camera/RPC_XML.h>
#include <boost/date_time/posix_time/posix_time.hpp>

namespace vw {
namespace camera {

  /// This is a generic line scan camera model that can be derived
  /// from to help implement specific cameras.  Some parts (velocity
  /// and atmospheric correction) currently only work for Earth.
  ///
  /// This expects the pose to be a rotation from the camera frame to
  /// the world frame. The position is a the camera's location in the
  /// world frame.
  ///
  /// The intrinisic model expects +Z to point out the camera. +X is
  /// the column direction of the image and is perpendicular to
  /// direction of flight. +Y is the row direction of the image (down
  /// the image); it is also the flight direction.  If this is not 
  /// accurate for your camera you can apply a rotation in PoseFuncT.

  inline vw::Quat get_adj(double s, std::vector<double> const& coeffs, int mode) {

    //std::cout << "--size is " << coeffs.size() << std::endl;
    double angle = coeffs[0];
    //std::cout << "--fist coeff is " << coeffs[0] << std::endl;
    for (int it = 0; it < (coeffs.size() - 1)/2; it++) {
      //std::cout << "--it is " << it << std::endl;
      
      double cos_coeff = coeffs[ 2 * it + 1];
      double sin_coeff = coeffs[ 2 * it + 2];
      
      //std::cout << "---coeff " << cos_coeff << ' ' << sin_coeff << std::endl;
      
      angle += cos_coeff * cos(s * (it + 1.0)) + sin_coeff * sin(s * (it + 1.0));
    }
    
    vw::Matrix3x3 M = vw::math::identity_matrix<3>();
    if (mode == 1) {
      // Along track jitter (rotate around x)
      M(1, 1) = cos(angle); M(1, 2) = -sin(angle);
      M(2, 1) = sin(angle); M(2, 2) =  cos(angle);
    }else if (mode == 2){
      // Across track jitter (rotate around y)
      M(0, 0) = cos(angle); M(0, 2) = -sin(angle);
      M(2, 0) = sin(angle); M(2, 2) =  cos(angle);
    }else if (mode == 3){
      // rotate around z
      M(0, 0) = cos(angle); M(0, 1) = -sin(angle);
      M(1, 0) = sin(angle); M(1, 1) =  cos(angle);
    }
    
    //std::cout << "M = " << M << std::endl;
    
    return vw::Quat(M);
  }
  
  inline vw::Quat get_adj(double t, double tb, double te,
                          std::vector<double> const& coeffsx,
                          std::vector<double> const& coeffsy,
                          std::vector<double> const& coeffsz) {
    
    // Take care of when coeffs are not initialized
    if (coeffsx.empty() || coeffsy.empty() || coeffsz.empty()) {
      return vw::Quat(vw::math::identity_matrix<3>());
    }

    if (coeffsx.size() != coeffsy.size() || coeffsx.size() != coeffsz.size() ) {
      vw::vw_throw( vw::ArgumentErr() << "Must have same number of coeffs in x, y, and z.\n");
    }

    //std::cout << "t is " << t << std::endl;
    //std::cout << "--tb and te are " << tb << ' ' << te << std::endl;
    
    // Normalize the value to [0, pi]. It may still go a little beyond it, depending on t.
    double s = M_PI * (t - tb)/(te - tb);
    
    return get_adj(s, coeffsx, 1) * get_adj(s, coeffsy, 2) * get_adj(s, coeffsz, 3);
  }

  class LinescanModelASP : public vw::camera::CameraModel {

  public:
    //------------------------------------------------------------------
    // Constructors / Destructors
    //------------------------------------------------------------------
    LinescanModelASP(Vector2i const& image_size,
                     bool            correct_velocity_aberration,
                     bool            correct_atmospheric_refraction,
                     double ms_offset):
      m_image_size(image_size), 
      m_correct_velocity_aberration(correct_velocity_aberration),
      m_correct_atmospheric_refraction(correct_atmospheric_refraction),
      m_ms_offset(ms_offset) {
//       std::cout << "--correct velocity aberration is " << m_correct_velocity_aberration
//                 << std::endl;
//       std::cout << "--correct atmospheric refraction " << m_correct_atmospheric_refraction
//                 << std::endl;
//       std::cout << "m_ms_offset " << m_ms_offset << std::endl;
      
      // Set default values for these constants which can be overridden later on.
      const double DEFAULT_EARTH_RADIUS      = 6371000.0;  // In meters.
      const double DEFAULT_SURFACE_ELEVATION = 0.0;
      m_mean_earth_radius      = DEFAULT_EARTH_RADIUS;
      m_mean_surface_elevation = DEFAULT_SURFACE_ELEVATION;
    }

    virtual ~LinescanModelASP() {}
    virtual std::string type() const { return "Linescan"; }

    //------------------------------------------------------------------
    // Interface
    //------------------------------------------------------------------
    
    // -- This set of functions implements virtual functions from CameraModel.h --
    
    /// Gives the camera position in world coordinates.
    virtual Vector3 camera_center(Vector2 const& pix) const {
      return get_camera_center_at_time(get_time_at_line(pix.y()));
    }

    /// Gives a pose vector which represents the rotation from camera to world units
    virtual vw::Quat camera_pose(Vector2 const& pix) const {
      return get_camera_pose_at_time(get_time_at_line(pix.y()));
    }

    // -- These are new functions --

    /// Returns the image size in pixels
    Vector2i get_image_size() const { return m_image_size; }
    
    int samples_per_line() const { return m_image_size[0]; }
    int number_of_lines () const { return m_image_size[1]; }
    
    /// Gives the camera velocity in world coordinates.
    Vector3 camera_velocity(vw::Vector2 const& pix) const {
      return get_camera_velocity_at_time(get_time_at_line(pix.y()));
    }
    
    // New functions for Linescan derived classes.
    // - Most of these deal with the fact that the camera is moving while
    //   acquiring the image.
    // - Derived classes need to implement these so that the functions inherited from
    //   CameraModel will work.  Consider using the functors in Extrinsics.h.

    /// Gives the camera center at a time.
    virtual Vector3 get_camera_center_at_time  (double time) const = 0;
    /// Gives the camera velocity at a time.
    virtual Vector3 get_camera_velocity_at_time(double time) const = 0;
    /// Get the pose at a time.
    virtual vw::Quat    get_camera_pose_at_time    (double time) const = 0;
    /// Return the computed time for a given line.
    virtual double  get_time_at_line           (double line) const = 0;

    // Here we use an initial guess for the line number
    // - This class provides a generic implementation but specific
    //   linescan cameras may be able to use more specific implementation.


    /// Get the pixel that observes a point in world coordinates.
    virtual Vector2 point_to_pixel(Vector3 const& point) const {
      return point_to_pixel(point, -1); // Redirect to the function with no guess
    }

    virtual Vector2 point_to_pixel(Vector3 const& point, double starty) const = 0;

/*
std::ostream& operator<<( std::ostream& os, LinescanModelASP const& camera_model) {
  os << "\n-------------------- Linescan Camera Model -------------------\n\n";
  os << " Number of Lines        :   " << camera_model.number_of_lines()        << "\n";
  os << " Samples per Line       :   " << camera_model.samples_per_line()       << "\n";
  os << "\n------------------------------------------------------------------------\n\n";
  return os;
}
*/

  public:

    double m_ms_offset; // the multispectral sensor line offset relative to pan
    
  public:

    /// Image size in pixels: [num lines, num samples]
    Vector2i m_image_size;      

    double m_mean_earth_radius;
    double m_mean_surface_elevation;

    /// Set this flag to enable velocity aberration correction.
    /// - For satellites this makes a big difference, make sure it is set!
    bool m_correct_velocity_aberration;
    
    /// Set this flag to enable atmospheric refraction correction.
    bool m_correct_atmospheric_refraction;

    /// Returns the radius of the Earth under the current camera position.
    double get_earth_radius() const;

  }; // End class LinescanModelASP
  
/*
  /// Output stream method for printing a summary of the linear
  /// pushbroom camera model parameters.
  std::ostream& operator<<( std::ostream& os, LinescanModelASP const& camera_model);
*/

}}      // namespace vw::camera

namespace asp {


  // The intrinisic model expects +Z to be point out the camera. +X is
  // the column direction of the image and is perpendicular to
  // direction of flight. +Y is the row direction of the image (down
  // the image); it is also the flight direction. This is different
  // from Digital Globe model, but you can rotate pose beforehand.

  // The standard class variant is:
  //      typedef LinescanDGModel<vw::camera::PiecewiseAPositionInterpolation,
  //                      			  vw::camera::SLERPPoseInterpolation> DGCameraModel;

  // The useful load_dg_camera_model() function is at the end of the file.

  /// Specialization of the generic LinescanModelASP for Digital Globe satellites.
  /// - Two template types are left floating so that AdjustedLinescanDGModel can modify them.
  template <class PositionFuncT, class PoseFuncT>
  class LinescanDGModel : public vw::camera::LinescanModelASP {
  public:
    //------------------------------------------------------------------
    // Constructors / Destructors
    //------------------------------------------------------------------
    LinescanDGModel(PositionFuncT const& position,
                    vw::camera::LinearPiecewisePositionInterpolation const& velocity,
                    PoseFuncT     const& pose,
                    vw::camera::TLCTimeInterpolation                 const& time,
                    vw::Vector2i  const& image_size,
                    vw::Vector2   const& detector_origin,
                    double        const  focal_length,
                    double        const  mean_ground_elevation=0,
                    bool                 correct_velocity=true,
                    bool                 correct_atmosphere=true,
                    double ms_offset = 0.0,
                    std::vector<double> coeffsx = std::vector<double>(),
                    std::vector<double> coeffsy = std::vector<double>(),
                    std::vector<double> coeffsz = std::vector<double>()):
      vw::camera::LinescanModelASP(image_size, correct_velocity, correct_atmosphere, ms_offset),
      m_position_func(position), m_velocity_func(velocity),
      m_pose_func(pose),         m_time_func(time),
      m_detector_origin(detector_origin),
      m_focal_length(focal_length),
      m_coeffsx(coeffsx), m_coeffsy(coeffsy), m_coeffsz(coeffsz) {
      m_mean_surface_elevation = mean_ground_elevation; // Set base class value

//       std::cout << "---correct_velocity = " << correct_velocity << std::endl;
//       std::cout << "--correct atmosphere = " << correct_atmosphere << std::endl;
      
    } 
    virtual ~LinescanDGModel() {}
    virtual std::string type() const { return "LinescanDG"; }
    
    // -- This set of functions implements virtual functions from LinescanModelASP.h --
    
    // Implement the functions from the LinescanModelASP class using functors
    virtual vw::Vector3 get_camera_center_at_time  (double time) const { return m_position_func(time); }
    virtual vw::Vector3 get_camera_velocity_at_time(double time) const { return m_velocity_func(time); }

    // A time-dependent rotation around the x axis
    // There are two copies of this!
    vw::Quat get_adj(double t) const {
      return vw::camera::get_adj(t, m_time_func(0), m_time_func(number_of_lines() - 1),
                                 m_coeffsx, m_coeffsy, m_coeffsz);
    }
    
    virtual vw::Quat     get_camera_pose_at_time    (double time) const { return m_pose_func(time) * get_adj(time); }

    virtual double      get_time_at_line           (double line) const { return m_time_func    (line - m_ms_offset); }
    
    /// As pixel_to_vector, but in the local camera frame.
    virtual vw::Vector3 get_local_pixel_vector(vw::Vector2 const& pix) const;
    
    // -- These are new functions --
    
    double       get_focal_length   () const {return m_focal_length;   } ///< Returns the focal length in pixels
    vw::Vector2  get_detector_origin() const {return m_detector_origin;} ///< Returns the detector origin in pixels    
    
    /// Create a fake pinhole model. It will return the same results
    /// as the linescan camera at current line y, but we will use it
    /// by extension at neighboring lines as well.
    /// TODO(oalexan1): Make this work for the MS band.
    vw::camera::PinholeModel linescan_to_pinhole(double y) const;
    
    
    PositionFuncT const& get_position_func() const {return m_position_func;} ///< Access the position function
    vw::camera::LinearPiecewisePositionInterpolation
    const& get_velocity_func() const {return m_velocity_func;} ///< Access the velocity function
    PoseFuncT     const& get_pose_func    () const {return m_pose_func;    } ///< Access the pose function
    vw::camera::TLCTimeInterpolation
    const& get_time_func    () const {return m_time_func;    } ///< Access the time     function
    
    
    virtual vw::Vector2 point_to_pixel(vw::Vector3 const& point) const;

    virtual vw::Vector2 point_to_pixel(vw::Vector3 const& point, double starty) const;
    
    virtual vw::Vector3 pixel_to_vector(vw::Vector2 const& pixel) const;
    
  public: // Functions
    
    /// Low accuracy function used by point_to_pixel to get a good solver starting seed.
    vw::Vector2 point_to_pixel_uncorrected(vw::Vector3 const& point, double starty) const;

  public: // Variables
  
    // Extrinsics
    PositionFuncT                                    m_position_func; ///< Yields position at time T
    vw::camera::LinearPiecewisePositionInterpolation m_velocity_func; ///< Yields velocity at time T
    PoseFuncT                                        m_pose_func;     ///< Yields pose     at time T
    vw::camera::TLCTimeInterpolation                 m_time_func;     ///< Yields time at a given line.

    std::vector<double> m_coeffsx, m_coeffsy, m_coeffsz; // Fourier coefficients
    
    // Intrinsics
    
    /// Location of (0,0) coordinate of the detector relative to the center of
    ///  the origin of the camera coordinate system.
    /// - Stored internally in pixels.
    vw::Vector2  m_detector_origin; 
    double       m_focal_length;    ///< The focal length, also stored in pixels.

    // Levenberg Marquardt solver for linescan number
    //
    // We solve for the line number of the image that position the
    // camera so that the projection into the camera model actually
    // hits the detector. The detector is normally offset in the y
    // direction on the optical plane. Once we have the line we don't
    // need to use a solver to compute the sample.
    // - This solver is used by the point_to_pixel_uncorrected function.
    class LinescanLMA : public vw::math::LeastSquaresModelBase<LinescanLMA> {
      const LinescanDGModel* m_model;
      vw::Vector3 m_point;
      double m_ms_offset;
    public:
      typedef vw::Vector<double> result_type;   // 1D error on the optical plane.
      typedef result_type        domain_type;   // 1D linescan number
      typedef vw::Matrix<double> jacobian_type;

      LinescanLMA(const LinescanDGModel* model, const vw::Vector3& pt, double ms_offset):
        m_model(model), m_point(pt), m_ms_offset(ms_offset) {}

      inline result_type operator()(domain_type const& y) const;
    };

  }; // End class LinescanDGModel


  /// This is the standard DG implementation
  typedef LinescanDGModel<vw::camera::PiecewiseAPositionInterpolation,
                  			  vw::camera::SLERPPoseInterpolation> DGCameraModel;

  /// Load a DG camera model from an XML file.
  /// - This function does not take care of Xerces XML init/de-init, the caller must
  ///   make sure this is done before/after this function is called!
  inline boost::shared_ptr<DGCameraModel> load_dg_camera_model_from_xml(std::string const& path);


// -----------------------------------------------------------------
// LinescanDGModel class functions

template <class PositionFuncT, class PoseFuncT>
vw::camera::PinholeModel LinescanDGModel<PositionFuncT, PoseFuncT>
::linescan_to_pinhole(double y) const {

  double t = this->get_time_at_line(y);
  return vw::camera::PinholeModel(this->get_camera_center_at_time(t),
                                  this->get_camera_pose_at_time(t).rotation_matrix(),
				  this->m_focal_length, -this->m_focal_length,
				  -this->m_detector_origin[0], y - this->m_detector_origin[1]
				  );
}


template <class PositionFuncT, class PoseFuncT>
vw::Vector3 LinescanDGModel<PositionFuncT, PoseFuncT>
::get_local_pixel_vector(vw::Vector2 const& pix) const {
  vw::Vector3 local_vec(pix[0] + m_detector_origin[0],
                        m_detector_origin[1] + m_ms_offset,
                        m_focal_length);
  return normalize(local_vec);
}

// Here we use an initial guess for the line number
template <class PositionFuncT, class PoseFuncT>
vw::Vector2 LinescanDGModel<PositionFuncT, PoseFuncT>
::point_to_pixel(vw::Vector3 const& point, double starty) const {

  // Use the uncorrected function to get a fast but good starting seed.
  vw::camera::CameraGenericLMA model( this, point );
  int status;
  vw::Vector2 start = point_to_pixel_uncorrected(point, starty);

  // Run the solver
  vw::Vector3 objective(0, 0, 0);
  const double ABS_TOL = 1e-16;
  const double REL_TOL = 1e-16;
  const int    MAX_ITERATIONS = 1e+5;
  vw::Vector2 solution = vw::math::levenberg_marquardtFixed<vw::camera::CameraGenericLMA, 2,3>
    (model, start, objective, status, ABS_TOL, REL_TOL, MAX_ITERATIONS);
  VW_ASSERT( status > 0,
          vw::camera::PointToPixelErr() << "Unable to project point into LinescanDG model" );

  return solution;
}

// Redirect to the function with no guess
template <class PositionFuncT, class PoseFuncT>
vw::Vector2 LinescanDGModel<PositionFuncT, PoseFuncT>
::point_to_pixel(vw::Vector3 const& point) const {
  return point_to_pixel(point, -1); 
}
  
// -----------------------------------------------------------------
// LinescanDGModel solver functions

// Function to minimize with the no-correction LMA optimizer.
template <class PositionFuncT, class PoseFuncT>
typename LinescanDGModel<PositionFuncT, PoseFuncT>::LinescanLMA::result_type
LinescanDGModel<PositionFuncT, PoseFuncT>::LinescanLMA::operator()( domain_type const& y ) const {

  // Get point in camera's frame and rescale to pixel units
  double       t        = m_model->get_time_at_line(y[0]);
  vw::Quat     pose     = m_model->get_camera_pose_at_time(t); // camera to world
  vw::Vector3  position = m_model->get_camera_center_at_time(t);
  vw::Vector3 pt        = inverse(pose).rotate(m_point - position);
  pt                   *= m_model->m_focal_length / pt.z();

  result_type result(1);
  // Error against the location of the detector.
  result[0] = pt.y() - m_model->m_detector_origin[1] - m_ms_offset; 
  return result;
}


// Computing the uncorrected pixel location is much faster.
template <class PositionFuncT, class PoseFuncT>
vw::Vector2 LinescanDGModel<PositionFuncT, PoseFuncT>::
point_to_pixel_uncorrected(vw::Vector3 const& point, double starty) const {

  // Solve for the correct line number to use
  LinescanLMA model(this, point, m_ms_offset);
  int status;
  vw::Vector<double> objective(1), start(1);
  start[0] = m_image_size.y()/2; 
  // Use a refined guess, if available, otherwise the center line.
  if (starty >= 0)
    start[0] = starty;

  // Run the solver
  const double ABS_TOL = 1e-16;
  const double REL_TOL = 1e-16;
  const int    MAX_ITERATIONS = 1e+5;
  vw::Vector<double> solution = vw::math::levenberg_marquardt(model, start, objective, status,
                                                              ABS_TOL, REL_TOL, MAX_ITERATIONS);

  VW_ASSERT( status > 0, vw::camera::PointToPixelErr()
             << "Unable to project point into LinescanDG model" );

  // Solve for sample location now that we know the correct line
  double       t        = this->get_time_at_line(solution[0]);
  vw::Quat     pose     = this->get_camera_pose_at_time(t); // camera to world
  vw::Vector3  position = this->get_camera_center_at_time(t);
  vw::Vector3 pt        = inverse(pose).rotate(point - position);
  pt                   *= this->m_focal_length / pt.z();

  return vw::Vector2(pt.x() - m_detector_origin[0], solution[0]);
}

template <class PositionFuncT, class PoseFuncT>
vw::Vector3 LinescanDGModel<PositionFuncT, PoseFuncT>::
pixel_to_vector(vw::Vector2 const& pixel) const {
  try {
        // Compute local vector from the pixel out of the sensor
        // - m_detector_origin and m_focal_length have been converted into units of pixels
        vw::Vector3 local_vec = get_local_pixel_vector(pixel);
        // Put the local vector in world coordinates using the pose information.
        double       t        = this->get_time_at_line(pixel.y());
        vw::Quat     pose     = this->get_camera_pose_at_time(t); // camera to world
        vw::Vector3 output_vector = pose.rotate(local_vec);


        vw::Vector3 cam_ctr = camera_center(pixel);
        if (!m_correct_atmospheric_refraction) 
          output_vector
            = vw::camera::apply_atmospheric_refraction_correction(cam_ctr,
                                                                  m_mean_earth_radius,
                                                                  m_mean_surface_elevation,
                                                                  output_vector);
        
        if (!m_correct_velocity_aberration) 
          return output_vector;
        else
          return vw::camera::apply_velocity_aberration_correction(cam_ctr, camera_velocity(pixel),
                                                                  m_mean_earth_radius,
                                                                  output_vector);

      } catch(const vw::Exception &e) {
        // Repackage any of our exceptions thrown below this point as a 
        //  pixel to ray exception that other code will be able to handle.
        vw_throw(vw::camera::PixelToRayErr() << e.what());
      }
    }
  
// -----------------------------------------------------------------
// LinescanDGModel supporting functions

class SecondsFrom
{
  boost::posix_time::ptime m_reference;
public:
  inline SecondsFrom( boost::posix_time::ptime const& time ) : m_reference(time) {}

  inline double operator()( boost::posix_time::ptime const& time ) const {
    return double( (time - m_reference).total_microseconds() ) / 1e6;
  }
};

inline boost::posix_time::ptime parse_time(std::string str){
  try{
    return boost::posix_time::time_from_string(str);
  }catch(...){
    vw::vw_throw(vw::ArgumentErr() << "Failed to parse time from string: " << str
		 << ". If you are not using Digital Globe images, you may "
                 << "need to specify the session type, such as -t rpc, "
                 << "-t rpcmaprpc, -t aster, etc.\n");
  }
  return boost::posix_time::time_from_string(str); // Never reached!
}

inline void parse_coeffs(std::vector<double> & coeffs, std::string const& var,
                           std::string const& band_id) {

  coeffs.clear();
  char * var_ptr = getenv(var.c_str());
  if (var_ptr != NULL && band_id != "NoAdj") {
    std::istringstream ifs(var_ptr);
    double val;
    std::cout << var << " = ";
    while (ifs >> val) {
      coeffs.push_back(val);
      std::cout << val << " ";
    }
    std::cout << std::endl;
  }
  std::cout << "size of coeffs for " << var << " = " << coeffs.size() << std::endl;
}
  
boost::shared_ptr<DGCameraModel> load_dg_camera_model_from_xml(std::string const& path){
  //vw_out() << "DEBUG - Loading DG camera file: " << camera_file << std::endl;

  // Parse the Digital Globe XML file
  GeometricXML geo;
  AttitudeXML  att;
  EphemerisXML eph;
  ImageXML     img;
  RPCXML       rpc;

  double ms_offset = 0;

  try {
    read_xml(path, geo, att, eph, img, rpc);
  } catch ( const std::exception& e ){
    vw::vw_throw(vw::ArgumentErr() << "Invalid Digital Globe XML file: " << path
		 << ". If you are not using Digital Globe images, you may need "
                 << "to specify the session type, such as -t rpc, -t rpcmaprpc, -t aster, etc.\n"
                 << e.what() << "\n");
  }

  std::string band_id = img.band_id;
  if (band_id == "Multi") {
    char * ptr = getenv("MS_OFFSET");
    if (ptr == NULL) {
      vw::vw_throw(vw::ArgumentErr() << "Must set MS_OFFSET.");
    }
    ms_offset = atof(ptr);
  }
  
  if (band_id == "Multib7") {
    char * ptr = getenv("MS_OFFSET");
    if (ptr == NULL) {
      vw::vw_throw(vw::ArgumentErr() << "Must set MS_OFFSET.");
    }
    ms_offset = -atof(ptr);
  }
  
  if (band_id == "Multib8") {
    char * ptr = getenv("MS_OFFSET");
    if (ptr == NULL) {
      vw::vw_throw(vw::ArgumentErr() << "Must set MS_OFFSET.");
    }
    ms_offset = atof(ptr);
  }
  
  std::cout << "--band id is " << band_id << std::endl;
  std::cout << "--ms offset is " << ms_offset << std::endl;

  std::vector<double> coeffsx, coeffsy, coeffsz;
  parse_coeffs(coeffsx, "ADJX", band_id);
  parse_coeffs(coeffsy, "ADJY", band_id);
  parse_coeffs(coeffsz, "ADJZ", band_id);
  
  // Get an estimate of the surface elevation from the corners specified in the file.
  // - Not every file has this information, in which case we will just use zero.
  double mean_ground_elevation = 0;
  vw::BBox3 bbox = rpc.get_lon_lat_height_box();
  if (!bbox.empty())
    mean_ground_elevation = (bbox.min()[2] + bbox.max()[2]) / 2.0;
  
  // Convert measurements in millimeters to pixels.
  geo.principal_distance /= geo.detector_pixel_pitch;
  geo.detector_origin    /= geo.detector_pixel_pitch;

//   std::cout << "--det origin in pixels " << geo.detector_origin << std::endl;
  
  // Convert all time measurements to something that boost::date_time can read.
  boost::replace_all( eph.start_time,            "T", " " );
  boost::replace_all( img.tlc_start_time,        "T", " " );
  boost::replace_all( img.first_line_start_time, "T", " " );
  boost::replace_all( att.start_time,            "T", " " );

  // Convert UTC time measurements to line measurements. Ephemeris
  // start time will be our reference frame to calculate seconds against.
  SecondsFrom convert( parse_time( eph.start_time ) );

  // I'm going make the assumption that EPH and ATT are sampled at the same rate and time.
  VW_ASSERT( eph.position_vec.size() == att.quat_vec.size(),
	     vw::MathErr() << "Ephemeris and attitude don't have the same number of samples." );
  VW_ASSERT( eph.start_time == att.start_time && eph.time_interval == att.time_interval,
	     vw::MathErr() << "Ephemeris and attitude don't seem to sample "
             << "with the same t0 or dt." );

  // Convert ephemeris to be position of camera. Change attitude to
  // be the rotation from camera frame to world frame. We also add an
  // additional rotation to the camera frame so X is the horizontal
  // direction to the picture and +Y points down the image (in the direction of flight).
  vw::Quat sensor_coordinate = vw::math::euler_xyz_to_quaternion
    (vw::Vector3(0,0,geo.detector_rotation - M_PI/2));
  for ( size_t i = 0; i < eph.position_vec.size(); i++ ) {
    eph.position_vec[i] += att.quat_vec[i].rotate( geo.perspective_center );
    att.quat_vec[i] = att.quat_vec[i] * geo.camera_attitude * sensor_coordinate;
  }

  // this must be independent on the observation
  // vw::vw_out() << "DG model load: sensor_coordinate = " << sensor_coordinate << std::endl;

  // Also eph and att must also be independent of the observation
  
  //geo.printDebugInfo(); // DEBUG INFO

  // Load up the time interpolation class. If the TLCList only has
  // one entry ... then we have to manually drop in the slope and offset.
  if ( img.tlc_vec.size() == 1 ) {
    double direction = 1;
    if ( boost::to_lower_copy( img.scan_direction ) != "forward" ) {
      direction = -1;
    }
    img.tlc_vec.push_back( std::make_pair(img.tlc_vec.front().first +
					  img.avg_line_rate, direction) );
  }

  // Build the TLCTimeInterpolation object and do a quick sanity check.
  vw::camera::TLCTimeInterpolation tlc_time_interpolation
    (img.tlc_vec, convert(parse_time(img.tlc_start_time)));
  VW_ASSERT(fabs(convert(parse_time(img.first_line_start_time)) - tlc_time_interpolation(0))
            < fabs(1.0 / (10.0 * img.avg_line_rate)),
	     vw::MathErr()
	     << "First Line Time and output from TLC lookup table "
	     << "do not agree of the ephemeris time for the first line of the image. "
	     << "If your XML camera files are not from the WorldView satellites, "
	     << "you may try the switch -t rpc to use the RPC camera model.\n"
	     << "The first image line ephemeris time is: "
  	     << convert( parse_time( img.first_line_start_time ) ) << ".\n"
	     << "The TLC look up table time is: " << tlc_time_interpolation( 0 ) << ".\n"
	     << "Maximum allowed difference is 1/10 of avg line rate, which is: "
	     << fabs( 1.0 / (10.0 * img.avg_line_rate ))
	     << ".\n"
    );
   
  vw::Vector2 final_detector_origin
    = subvector(inverse(sensor_coordinate).rotate(vw::Vector3(geo.detector_origin[0],
							      geo.detector_origin[1],
							      0)), 0, 2);

  // This must be independent of the observation
//   std::cout << "--final det origin is " << final_detector_origin << std::endl;

//   std::cout << "--geo principal distance: " << geo.principal_distance << std::endl;
  
  double et0 = convert( parse_time( eph.start_time ) );
  double at0 = convert( parse_time( att.start_time ) );
  double edt = eph.time_interval;
  double adt = att.time_interval;

//   std::cout << "--xx27 eph " << eph.position_vec[0] << ' ' << eph.velocity_vec[0] << ' '
//             << et0 << ' ' << edt << std::endl;

//   std::cout << "--att " << att.quat_vec[0] << ' ' << at0 << ' ' << adt << std::endl;

//   std::cout << "--geo principal distance " << geo.principal_distance << std::endl;

//   std::cout << "--mean ground elevation " << mean_ground_elevation << std::endl;

//   std::cout << "--image size " << img.image_size << std::endl;
  
  // This is where we could set the Earth radius if we have that info.

  typedef boost::shared_ptr<DGCameraModel> CameraModelPtr;
  return CameraModelPtr(new DGCameraModel
                        (vw::camera::PiecewiseAPositionInterpolation(eph.position_vec,
                                                                     eph.velocity_vec, et0, edt ),
                         vw::camera::LinearPiecewisePositionInterpolation(eph.velocity_vec,
                                                                          et0, edt),
                         vw::camera::SLERPPoseInterpolation(att.quat_vec, at0, adt),
                         tlc_time_interpolation, img.image_size,
                         final_detector_origin,
                         geo.principal_distance,
                         mean_ground_elevation,
                         !stereo_settings().disable_correct_velocity_aberration,
                         !stereo_settings().disable_correct_atmospheric_refraction, ms_offset,
                         coeffsx, coeffsy, coeffsz));
} // End function load_dg_camera_model()
  

// ----------------- Start here

class LinescanModelFreq: public vw::camera::CameraModel {
  
public:
  DGCameraModel * m_cam;
  std::vector<double> m_coeffsx, m_coeffsy, m_coeffsz; // Fourier coefficients
  double m_ms_offset;
  
public:
  LinescanModelFreq(DGCameraModel * cam): m_cam(cam), m_ms_offset(0) {}
  
  virtual ~LinescanModelFreq() {}
  virtual std::string type() const { return "LinescanModelFreq"; }
  
  // A time-dependent rotation around the x axis
  // There are two copies of this!
  vw::Quat get_adj(double t) const {
    return vw::camera::get_adj(t, m_cam->m_time_func(0), m_cam->m_time_func(number_of_lines() - 1),
                               m_coeffsx, m_coeffsy, m_coeffsz);
  }
  
  vw::Vector3 get_local_pixel_vector(vw::Vector2 const& pix) const {
    vw::Vector3 local_vec(pix[0] + m_cam->m_detector_origin[0],
                          m_cam->m_detector_origin[1] + m_ms_offset,
                          m_cam->m_focal_length);
    return normalize(local_vec);
  }
  
  /// Gives the camera position in world coordinates.
  virtual vw::Vector3 camera_center(vw::Vector2 const& pix) const {
    return m_cam->get_camera_center_at_time(this->get_time_at_line(pix.y()));
  }

  /// Gives a pose vector which represents the rotation from camera to world units
  virtual vw::Quat camera_pose(vw::Vector2 const& pix) const {
    if (!m_coeffsx.empty() && !m_cam->m_coeffsx.empty()) {
      vw_throw(vw::camera::PixelToRayErr() << "Cannot have both the unadjusted and "
               << "adjusted cameras have frequency modulation");
    }
    return this->get_camera_pose_at_time(this->get_time_at_line(pix.y()));
  }
  
  /// Returns the image size in pixels
  vw::Vector2i get_image_size() const { return m_cam->m_image_size; }
  
  int samples_per_line() const { return m_cam->m_image_size[0]; }
  int number_of_lines () const { return m_cam->m_image_size[1]; }
  
  /// Gives the camera velocity in world coordinates.
  vw::Vector3 camera_velocity(vw::Vector2 const& pix) const {
    return this->get_camera_velocity_at_time(this->get_time_at_line(pix.y()));
  }
  
  /// Gives the camera center at a time.
  virtual vw::Vector3 get_camera_center_at_time(double time) const {
    return m_cam->get_camera_center_at_time(time);
  }
  
  /// Gives the camera velocity at a time.
  virtual vw::Vector3 get_camera_velocity_at_time(double time) const {
    return m_cam->get_camera_velocity_at_time(time);
  }
  
  /// Get the pose at a time.
  virtual vw::Quat get_camera_pose_at_time(double time) const {
    return m_cam->get_camera_pose_at_time(time) * get_adj(time);
  }

  /// Return the computed time for a given line.
  virtual double get_time_at_line(double line) const {
    return m_cam->m_time_func(line - m_ms_offset);
  }

  class LinescanFreqLMA : public vw::math::LeastSquaresModelBase<LinescanFreqLMA> {
    const LinescanModelFreq* m_model;
    vw::Vector3 m_point;
    double m_ms_offset;
  public:
    typedef vw::Vector<double> result_type;   // 1D error on the optical plane.
    typedef result_type        domain_type;   // 1D linescan number
    typedef vw::Matrix<double> jacobian_type;
    
    LinescanFreqLMA(const LinescanModelFreq* model, const vw::Vector3& pt, double ms_offset):
      m_model(model), m_point(pt), m_ms_offset(ms_offset) {}
    
    inline result_type operator()(domain_type const& y) const {
      // Get point in camera's frame and rescale to pixel units
      double       t        = m_model->get_time_at_line(y[0]);
      vw::Quat     pose     = m_model->get_camera_pose_at_time(t); // camera to world
      vw::Vector3  position = m_model->get_camera_center_at_time(t);
      vw::Vector3 pt        = inverse(pose).rotate(m_point - position);
      pt                   *= m_model->m_cam->m_focal_length / pt.z();
      
      result_type result(1);
      // Error against the location of the detector.
      result[0] = pt.y() - m_model->m_cam->m_detector_origin[1] - m_ms_offset; 
      return result;
    }
  };

  vw::Vector2 point_to_pixel_uncorrected(vw::Vector3 const& point, double starty) const {
    
    // Solve for the correct line number to use
    LinescanFreqLMA model(this, point, m_ms_offset);
    int status;
    vw::Vector<double> objective(1), start(1);
    start[0] = m_cam->m_image_size.y()/2; 
    // Use a refined guess, if available, otherwise the center line.
    if (starty >= 0)
      start[0] = starty;
    
    // Run the solver
    const double ABS_TOL = 1e-16;
    const double REL_TOL = 1e-16;
    const int    MAX_ITERATIONS = 1e+5;
    vw::Vector<double> solution = vw::math::levenberg_marquardt(model, start, objective, status,
                                                                ABS_TOL, REL_TOL, MAX_ITERATIONS);
    
    VW_ASSERT( status > 0, vw::camera::PointToPixelErr()
               << "Unable to project point into LinescanDG model" );
    
    // Solve for sample location now that we know the correct line
    double       t        = this->get_time_at_line(solution[0]);
    vw::Quat     pose     = this->get_camera_pose_at_time(t); // camera to world
    vw::Vector3  position = this->get_camera_center_at_time(t);
    vw::Vector3 pt        = inverse(pose).rotate(point - position);
    pt                   *= this->m_cam->m_focal_length / pt.z();
    
    return vw::Vector2(pt.x() - m_cam->m_detector_origin[0], solution[0]);
  }

    
  // Here we use an initial guess for the line number
  virtual vw::Vector2 point_to_pixel(vw::Vector3 const& point, double starty) const {
    
    // Use the uncorrected function to get a fast but good starting seed.
    vw::camera::CameraGenericLMA model(this, point);
    int status;
    vw::Vector2 start = point_to_pixel_uncorrected(point, starty);
    
    // Run the solver
    vw::Vector3 objective(0, 0, 0);
    const double ABS_TOL = 1e-16;
    const double REL_TOL = 1e-16;
    const int    MAX_ITERATIONS = 1e+5;
    vw::Vector2 solution = vw::math::levenberg_marquardtFixed<vw::camera::CameraGenericLMA, 2,3>
      (model, start, objective, status, ABS_TOL, REL_TOL, MAX_ITERATIONS);
    VW_ASSERT( status > 0,
               vw::camera::PointToPixelErr() << "Unable to project point into LinescanDG model" );
    
    return solution;
  }
  
  /// Get the pixel that observes a point in world coordinates.
  virtual vw::Vector2 point_to_pixel(vw::Vector3 const& point) const {
    return point_to_pixel(point, -1); // Redirect to the function with no guess
  }
  
  virtual vw::Vector3 pixel_to_vector(vw::Vector2 const& pixel) const {
    try {
      // Compute local vector from the pixel out of the sensor
      // - m_detector_origin and m_focal_length have been converted into units of pixels
      vw::Vector3 local_vec = this->get_local_pixel_vector(pixel);
      // Put the local vector in world coordinates using the pose information.
      double       t        = this->get_time_at_line(pixel.y());
      vw::Quat     pose     = this->get_camera_pose_at_time(t); // camera to world
      vw::Vector3 output_vector = pose.rotate(local_vec);
      
      vw::Vector3 cam_ctr = this->camera_center(pixel);
      if (!m_cam->m_correct_atmospheric_refraction) 
        output_vector
          = vw::camera::apply_atmospheric_refraction_correction(cam_ctr,
                                                                m_cam->m_mean_earth_radius,
                                                                m_cam->m_mean_surface_elevation,
                                                                output_vector);
      
      if (!m_cam->m_correct_velocity_aberration) 
        return output_vector;
      else
        return vw::camera::apply_velocity_aberration_correction(cam_ctr,
                                                                this->camera_velocity(pixel),
                                                                m_cam->m_mean_earth_radius,
                                                                output_vector);
      
    } catch(const vw::Exception &e) {
      // Repackage any of our exceptions thrown below this point as a 
      //  pixel to ray exception that other code will be able to handle.
      vw_throw(vw::camera::PixelToRayErr() << e.what());
    }
  }
};
  
}      // namespace asp
  
  
#endif//__STEREO_CAMERA_LINESCAN_DG_MODEL_H__
