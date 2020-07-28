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

  class LinescanModelASP : public vw::camera::CameraModel {

  public:
    //------------------------------------------------------------------
    // Constructors / Destructors
    //------------------------------------------------------------------
    LinescanModelASP(Vector2i const& image_size,
                     bool            correct_velocity_aberration,
                     bool            correct_atmospheric_refraction,
                     double ms_offset) :
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
    virtual Quat camera_pose(Vector2 const& pix) const {
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
    virtual Quat    get_camera_pose_at_time    (double time) const = 0;
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
    
  protected:

    /// Image size in pixels: [num lines, num samples]
    Vector2i m_image_size;      

    double m_mean_earth_radius;
    double m_mean_surface_elevation;

    /// Set this flag to enable velocity aberration correction.
    /// - For satellites this makes a big difference, make sure it is set!
    bool m_correct_velocity_aberration;
    
    /// Set this flag to enable atmospheric refraction correction.
    bool m_correct_atmospheric_refraction;

    double m_ms_offset; // the multispectral sensor line offset relative to pan

  protected:

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
                    double ms_offset = 0.0):
      vw::camera::LinescanModelASP(image_size, correct_velocity, correct_atmosphere, ms_offset),
      m_position_func(position), m_velocity_func(velocity),
      m_pose_func(pose),         m_time_func(time),
      m_detector_origin(detector_origin),
      m_focal_length(focal_length) {
      m_mean_surface_elevation = mean_ground_elevation; // Set base class value
    } 
    virtual ~LinescanDGModel() {}
    virtual std::string type() const { return "LinescanDG"; }
    
    // -- This set of functions implements virtual functions from LinescanModelASP.h --
    
    // Implement the functions from the LinescanModelASP class using functors
    virtual vw::Vector3 get_camera_center_at_time  (double time) const { return m_position_func(time); }
    virtual vw::Vector3 get_camera_velocity_at_time(double time) const { return m_velocity_func(time); }
    virtual vw::Quat    get_camera_pose_at_time    (double time) const { return m_pose_func    (time); }
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
    PoseFuncT     const& get_pose_func    () const {return m_pose_func;    } ///< Access the pose     function
    vw::camera::TLCTimeInterpolation
    const& get_time_func    () const {return m_time_func;    } ///< Access the time     function
    
    
    virtual vw::Vector2 point_to_pixel(vw::Vector3 const& point, double starty) const;
    
    virtual vw::Vector3 pixel_to_vector(vw::Vector2 const& pixel) const;
    
    
  protected: // Functions
    
    /// Low accuracy function used by point_to_pixel to get a good solver starting seed.
    vw::Vector2 point_to_pixel_uncorrected(vw::Vector3 const& point, double starty) const;

  protected: // Variables
  
    // Extrinsics
    PositionFuncT                                    m_position_func; ///< Yields position at time T
    vw::camera::LinearPiecewisePositionInterpolation m_velocity_func; ///< Yields velocity at time T
    PoseFuncT                                        m_pose_func;     ///< Yields pose     at time T
    vw::camera::TLCTimeInterpolation                 m_time_func;     ///< Yields time at a given line.

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

}      // namespace asp

#include <asp/Camera/LinescanDGModel.tcc>


#endif//__STEREO_CAMERA_LINESCAN_DG_MODEL_H__
