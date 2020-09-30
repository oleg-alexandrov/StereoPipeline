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



#include <vw/Math/EulerAngles.h>
#include <vw/Camera/CameraSolve.h>
#include <asp/Core/StereoSettings.h>
#include <asp/Camera/RPCModel.h>
#include <asp/Camera/RPC_XML.h>
#include <boost/date_time/posix_time/posix_time.hpp>

namespace vw {
namespace camera {

}} // namespace vw::camera


namespace asp {

// -----------------------------------------------------------------
// LinescanDGModel class functions

template <class PositionFuncT, class PoseFuncT>
vw::camera::PinholeModel LinescanDGModel<PositionFuncT, PoseFuncT>
::linescan_to_pinhole(double y) const {

  double t = this->get_time_at_line(y);
  return vw::camera::PinholeModel(this->m_position_func(t),  this->get_camera_pose_at_time(t).rotation_matrix(),
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

// -----------------------------------------------------------------
// LinescanDGModel solver functions

// Function to minimize with the no-correction LMA optimizer.
template <class PositionFuncT, class PoseFuncT>
typename LinescanDGModel<PositionFuncT, PoseFuncT>::LinescanLMA::result_type
LinescanDGModel<PositionFuncT, PoseFuncT>::LinescanLMA::operator()( domain_type const& y ) const {

  // Get point in camera's frame and rescale to pixel units
  double       t        = m_model->get_time_at_line(y[0]);
  vw::Quat     pose     = m_model->get_camera_pose_at_time(t); // camera to world
  vw::Vector3  position = m_model->m_position_func(t);
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
  vw::Vector3  position = this->m_position_func(t);
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
  std::cout << "--band id is " << band_id << std::endl;
  if (band_id == "Multi") {
    char * ptr = getenv("MS_OFFSET");
    if (ptr == NULL) {
      vw::vw_throw(vw::ArgumentErr() << "Must set MS_OFFSET.");
    }
    ms_offset = atof(ptr);
  }
  std::cout << "--ms offset is " << ms_offset << std::endl;
  std::vector<double> coeffs;
  char * adj = getenv("ADJ");
  std::cout << "--adj is " << adj << std::endl;
  if (adj != NULL) {
    coeffs.clear();
    std::istringstream ifs(adj);
    double val;
    while (ifs >> val) {
      coeffs.push_back(val);
      std::cout << "--coeff is " << val << std::endl;
    }
  }
  
  // Get an estimate of the surface elevation from the corners specified in the file.
  // - Not every file has this information, in which case we will just use zero.
  double mean_ground_elevation = 0;
  vw::BBox3 bbox = rpc.get_lon_lat_height_box();
  if (!bbox.empty())
    mean_ground_elevation = (bbox.min()[2] + bbox.max()[2]) / 2.0;
  
  // Convert measurements in millimeters to pixels.
  geo.principal_distance /= geo.detector_pixel_pitch;
  geo.detector_origin    /= geo.detector_pixel_pitch;

  std::cout << "--det origin in pixels " << geo.detector_origin << std::endl;
  
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
  vw::vw_out() << "DG model load: sensor_coordinate = " << sensor_coordinate << std::endl;

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
  std::cout << "--final det origin is " << final_detector_origin << std::endl;

  std::cout << "--geo principal distance: " << geo.principal_distance << std::endl;
  
  double et0 = convert( parse_time( eph.start_time ) );
  double at0 = convert( parse_time( att.start_time ) );
  double edt = eph.time_interval;
  double adt = att.time_interval;

  std::cout << "--xx27 eph " << eph.position_vec[0] << ' ' << eph.velocity_vec[0] << ' '
            << et0 << ' ' << edt << std::endl;

  std::cout << "--att " << att.quat_vec[0] << ' ' << at0 << ' ' << adt << std::endl;

  std::cout << "--geo principal distance " << geo.principal_distance << std::endl;

  std::cout << "--mean ground elevation " << mean_ground_elevation << std::endl;

  std::cout << "--image size " << img.image_size << std::endl;
  
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
                         coeffs));
} // End function load_dg_camera_model()
  
  
} // end namespace asp

