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

/// \file jitter_solve.cc
///

#include <asp/Core/BundleAdjustUtils.h>
#include <asp/Core/InterestPointMatching.h>
#include <asp/Core/Macros.h>
#include <asp/Core/StereoSettings.h>
#include <asp/Sessions/StereoSession.h>
#include <asp/Sessions/StereoSessionFactory.h>
#include <asp/Camera/LinescanDGModel.h>
#include <vw/BundleAdjustment/ControlNetwork.h>
#include <vw/BundleAdjustment/ControlNetworkLoader.h>
#include <vw/Cartography/CameraBBox.h>
#include <vw/InterestPoint/InterestData.h>
#include <vw/Core/Stopwatch.h>
#include <xercesc/util/PlatformUtils.hpp>

#include <ceres/ceres.h>
#include <ceres/loss_function.h>

namespace po = boost::program_options;
namespace fs = boost::filesystem;
using namespace vw;

typedef boost::scoped_ptr<asp::StereoSession> SessionPtr;
typedef boost::shared_ptr<vw::camera::CameraModel> CamPtr;
typedef vw::ba::CameraRelationNetwork<vw::ba::JFeature> CRNJ;
typedef vw::ba::CameraNode<vw::ba::JFeature>::iterator crn_iter;

const int NUM_POINT_PARAMS = 3;
const int PIXEL_SIZE = 2;

struct Options : vw::cartography::GdalWriteOptions {
  std::string pan_image, b7_image, b8_image, pan_camera, b7_camera, b8_camera, out_prefix,
    images_to_use, mode;
  int num_frequencies; // TODO(oalexan1): Use instead the frequency, in Hz
  int num_iterations;
  double ms_offset, min_triangulation_angle, horiz_weight, vert_weight;
  bool read_xyz_on_dem, float_free_term, float_ms_offset;
};

/// A Ceres cost function. We pass in the observation and the model.
///  The result is the residual, the difference in the observation 
///  and the projection of the point into the camera.
struct ReprojectionError {
  ReprojectionError(int sign, Vector2 const& observation, CamPtr cam, 
                    std::vector<int> const& block_sizes):
    m_sign(sign),
    m_observation(observation),
    m_cam(cam),
    m_block_sizes(block_sizes){}

  // Call to work with ceres::DynamicCostFunctions. Takes array of
  // arrays. Each such array is a parameter block of a size that is
  // determined at run-time.
  bool operator()(double const * const * parameters, double * residuals) const {

    // Need to create here the adjusted camera from the original camera

    for (size_t it = 0; it < m_block_sizes.size(); it++) {
//        std::cout << "--block size is " << it << " " << m_block_sizes[it] << std::endl;
    }

    asp::DGCameraModel * dg_cam = dynamic_cast<asp::DGCameraModel*>(m_cam.get());
    if (dg_cam == NULL) 
      vw::vw_throw( vw::ArgumentErr() << "Expecting a Digital Globe camera.\n\n");

    // Create the camera
    asp::LinescanModelFreq freq_cam(dg_cam);

//     std::cout << "--old offset " << m_sign << ' ' << dg_cam->m_ms_offset << std::endl;
    freq_cam.m_ms_offset = m_sign * parameters[4][0];
//     std::cout << "--new offset " << m_icam << ' ' << freq_cam.m_ms_offset << std::endl;
//     std::cout << std::endl;

    int num_coeffsx = m_block_sizes[0];
    freq_cam.m_coeffsx.resize(num_coeffsx);
    for (int it = 0; it < num_coeffsx; it++) freq_cam.m_coeffsx[it] = parameters[0][it];
    
    int num_coeffsy = m_block_sizes[1];
    freq_cam.m_coeffsy.resize(num_coeffsy);
    for (int it = 0; it < num_coeffsy; it++) freq_cam.m_coeffsy[it] = parameters[1][it];

    int num_coeffsz = m_block_sizes[2];
    freq_cam.m_coeffsz.resize(num_coeffsz);
    for (int it = 0; it < num_coeffsz; it++) freq_cam.m_coeffsz[it] = parameters[2][it];
    
    // Create the point
    Vector3 xyz;
    for (int it = 0; it < m_block_sizes[3]; it++) {
      xyz[it] = parameters[3][it];
    }
     //std::cout << "--got xyz = " << xyz << std::endl;
    
    try {
      // Use the camera model wrapper to handle all of the parameter blocks.
      Vector2 prediction = freq_cam.point_to_pixel(xyz);

       //std::cout << "Got prediction " << prediction << std::endl;
       //std::cout << "--got observation " << m_observation << std::endl;
      
      // The error is the difference between the predicted and observed position.
      residuals[0] = prediction[0] - m_observation[0]; // Input units are pixels
      residuals[1] = prediction[1] - m_observation[1];

       //std::cout << "Got residuals " << residuals[0] << ", " << residuals[1] << std::endl;

    } catch (std::exception const& e) { // TODO: Catch only projection errors?
      // Failed to compute residuals

      //std::cout << "Caught exception!" << std::endl;

      // Not sure if to return false or true here
      residuals[0] = 0;
      residuals[1] = 0;
      return false;
    }
    
    return true;
  }

  // Factory to hide the construction of the CostFunction object from the client code.
  static ceres::CostFunction* Create(int sign, Vector2 const& observation, CamPtr cam,
                                     std::vector<int> const& block_sizes){
    ceres::DynamicNumericDiffCostFunction<ReprojectionError>* cost_function =
      new ceres::DynamicNumericDiffCostFunction<ReprojectionError>
      (new ReprojectionError(sign, observation, cam, block_sizes));
    
    // The residual size is always the same.
    cost_function->SetNumResiduals(PIXEL_SIZE);

    // Set the size of each optimization block
    for (size_t i = 0; i < block_sizes.size(); i++) {
      cost_function->AddParameterBlock(block_sizes[i]);
    }
    return cost_function;
  }

private:
  int m_sign;
  Vector2 m_observation; // The pixel observation for this camera/point pair
  CamPtr m_cam;
  std::vector<int> m_block_sizes;
}; // End class ReprojectionError

// TODO(oalexan1): Expose the threshold
ceres::LossFunction* get_jitter_loss_function(){
  return new ceres::CauchyLoss(0.5);
}

// A ceres cost function. The residual is the vertical component of
// the difference between the observed 3D point and the current
// (floating) 3D point.
struct XYZVertError {
  XYZVertError(Vector3 const& obs, double vert_weight):
    m_obs(obs), m_vert_weight(vert_weight) {}

  // Call to work with ceres::DynamicCostFunctions. Takes array of
  // arrays. Each such array is a parameter block of a size that is
  // determined at run-time.
  bool operator()(double const * const * parameters, double * residuals) const {
    
    Vector3 xyz;
    for (size_t p = 0; p < 3; p++) 
      xyz[p] = parameters[0][p];

    Vector3 err = m_obs * dot_prod(xyz - m_obs, m_obs) / dot_prod(m_obs, m_obs);
    err *= m_vert_weight;
    
    for (size_t p = 0; p < 3; p++)
      residuals[p] = err[p];
    
    return true;
  }

  // Factory to hide the construction of the CostFunction object from the client code.
  static ceres::CostFunction* Create(Vector3 const& obs, double vert_weight){
    ceres::DynamicNumericDiffCostFunction<XYZVertError>* cost_function =
      new ceres::DynamicNumericDiffCostFunction<XYZVertError>
      (new XYZVertError(obs, vert_weight));
    
    // The residual size is always the same.
    cost_function->SetNumResiduals(3);
    cost_function->AddParameterBlock(3);
    return cost_function;
  }

private:
  Vector3 m_obs; 
  double m_vert_weight;
}; // End class XYZVertError

// A ceres cost function. The residual is the horizical component of
// the difference between the observed 3D point and the current
// (floating) 3D point.
struct XYZHorizError {
  XYZHorizError(Vector3 const& obs, double horiz_weight):
    m_obs(obs), m_horiz_weight(horiz_weight) {}

  // Call to work with ceres::DynamicCostFunctions. Takes array of
  // arrays. Each such array is a parameter block of a size that is
  // determined at run-time.
  bool operator()(double const * const * parameters, double * residuals) const {

    
    Vector3 xyz;
    for (size_t p = 0; p < 3; p++) 
      xyz[p] = parameters[0][p];

     Vector3 err = m_obs * dot_prod(xyz - m_obs, m_obs) / dot_prod(m_obs, m_obs);
     err = xyz - m_obs - err;
     err *= m_horiz_weight;
    
    for (size_t p = 0; p < 3; p++)
      residuals[p] = err[p];
    
    return true;
  }

  // Factory to hide the construction of the CostFunction object from the client code.
  static ceres::CostFunction* Create(Vector3 const& obs, double horiz_weight){
    ceres::DynamicNumericDiffCostFunction<XYZHorizError>* cost_function =
      new ceres::DynamicNumericDiffCostFunction<XYZHorizError>
      (new XYZHorizError(obs, horiz_weight));
    
    // The residual size is always the same.
    cost_function->SetNumResiduals(3);
    cost_function->AddParameterBlock(3);
    return cost_function;
  }

private:
  Vector3 m_obs; 
  double m_horiz_weight;
}; // End class XYZHorizError

/// Compute the residuals
void compute_residuals(Options const& opt,
                       std::vector<int> const& cam_residual_counts,
                       int num_points,
                       ceres::Problem &problem,
                       // output
                       std::vector<double> & residuals) {
  
  double cost = 0.0;
  ceres::Problem::EvaluateOptions eval_options;
  eval_options.apply_loss_function = false;
  eval_options.num_threads = opt.num_threads;
  
  problem.Evaluate(eval_options, &cost, &residuals, 0, 0);
  int num_residuals = residuals.size();
  
  // Verify our residual calculations are correct
  int num_expected_residuals = 0;
  for (size_t i = 0; i < cam_residual_counts.size(); i++) 
    num_expected_residuals += cam_residual_counts[i]*PIXEL_SIZE;

  if (opt.horiz_weight > 0 || opt.vert_weight > 0) {
    num_expected_residuals += num_points * NUM_POINT_PARAMS; // horiz residuals
    num_expected_residuals += num_points * NUM_POINT_PARAMS; // vert residuals
  }
  
  if (num_expected_residuals != num_residuals)
    vw_throw( LogicErr() << "Expected " << num_expected_residuals
              << " residuals but instead got " << num_residuals);
}

/// Compute residual map by averaging all the reprojection error at a given point
void compute_mean_residuals_at_xyz(CRNJ & crn,
                                   int num_points, int num_cameras,
                                   std::vector<double> const& residuals,
                                   // outputs
                                   std::vector<double> & mean_residuals,
                                   std::vector< std::vector<double> > & indiv_residuals,
                                   std::vector<int>  & num_point_observations) {

  mean_residuals.resize(num_points);
  num_point_observations.resize(num_points);

  // Allocate memory for individual residuals. Each camera has x and y residuals.
  indiv_residuals.resize(num_points);
  for (int ipt = 0; ipt < num_points; ipt++) {
    indiv_residuals[ipt].resize(num_cameras * PIXEL_SIZE, 0);
  }
    
  // Observation residuals are stored at the beginning of the residual vector in the 
  //  same order they were originally added to Ceres.
  
  size_t residual_index = 0;
  // Double loop through cameras and crn entries will give us the correct order
  for (int icam = 0; icam < num_cameras; icam++) {
    for (crn_iter fiter = crn[icam].begin(); fiter != crn[icam].end(); fiter++){

      // The index of the 3D point
      int ipt = (**fiter).m_point_id;

      // Get the residual error for this observation
      double errorX         = residuals[residual_index + 0];
      double errorY         = residuals[residual_index + 1];
      double residual_error = (std::abs(errorX) + std::abs(errorY)) / 2;
      residual_index += PIXEL_SIZE;

      indiv_residuals[ipt][PIXEL_SIZE * icam + 0] = errorX;
      indiv_residuals[ipt][PIXEL_SIZE * icam + 1] = errorY;
      
      // Update information for this point
      num_point_observations[ipt] += 1;
      mean_residuals        [ipt] += residual_error;
    }
  } // End double loop through all the observations

  // Do the averaging
  for (int i = 0; i < num_points; i++) {
    mean_residuals[i] /= static_cast<double>(num_point_observations[i]);
  }
  
} // End function compute_mean_residuals_at_xyz

/// Write out a .csv file recording the residual error at each location on the ground
void write_residual_map(// Mean residual of each point
                        std::vector<double> const& mean_residuals,
                        // Individual x and y residuals per camera
                        std::vector< std::vector<double> > const& indiv_residuals,
                         // Num non-outlier pixels per point
                        std::vector<int> const& num_point_observations,
                        int num_points, int num_cameras,
                        std::vector<double> const& points_vec,
                        vw::ba::ControlNetwork const& cnet,
                        vw::cartography::Datum const& datum,
                        std::string const& output_path,
                        Options const& opt) {

  if (cnet.size() != num_points || num_points * NUM_POINT_PARAMS != int(points_vec.size())) {
    vw_throw( LogicErr()
              << "The number of points points "
              << "does not agree with number of points in cnet.\n");
  }
  
  // Open the output file and write the header
  vw_out() << "Writing: " << output_path << std::endl;
  
  std::ofstream file(output_path.c_str());
  file.precision(18);
  file << "# lon, lat, height_above_datum, mean_residual, num_observations, indiv residuals\n";
  file << "# " << datum << std::endl;
  
  // Now write all the points to the file
  for (int i = 0; i < num_points; i++) {

    int k = NUM_POINT_PARAMS * i;
    Vector3 xyz(points_vec[k], points_vec[k + 1], points_vec[k + 2]);
    
    Vector3 llh = datum.cartesian_to_geodetic(xyz);
    
    file << llh[0] <<", "<< llh[1] <<", "<< llh[2] <<", "<< mean_residuals[i] <<", "
         << num_point_observations[i];
    
    for (size_t j = 0; j < indiv_residuals[i].size(); j++)
      file << ", " << indiv_residuals[i][j];

    file << "\n";
  }
  file.close();
  
} // End function write_residual_map

void save_residuals(Options const& opt, CRNJ & crn,
                    vw::ba::ControlNetwork const& cnet,
                    std::vector<double> const& points_vec,
                    ceres::Problem &problem,
                    vw::cartography::Datum const& datum,
                    int num_points, int num_cameras,
                    std::vector<int> const& cam_residual_counts,
                    std::string const& output_path) {
    
    std::vector<double> residuals;
    std::vector<double> mean_residuals;
    std::vector< std::vector<double> > indiv_residuals;
    std::vector<int> num_point_observations;
    
    compute_residuals(opt, cam_residual_counts, num_points, problem, residuals);
    compute_mean_residuals_at_xyz(crn, num_points, num_cameras, residuals,  
                                   // outputs
                                  mean_residuals, indiv_residuals, num_point_observations);

    // Write out a .csv file recording the residual error at each location on the ground
    write_residual_map(mean_residuals, indiv_residuals,  
                       num_point_observations, num_points, num_cameras,
                       points_vec,
                       cnet, datum, output_path, opt);
}

  
// Read xyz for each col and row. Average values. Assume that we will
// read at most two such files.
void read_xyz(std::map< std::pair<float, float>, Vector3> & xyz, std::string const& filename) {
  std::ifstream ifs(filename.c_str());

  float col, row;
  double x, y, z;
  while (ifs >> col >> row >> x >> y >> z) {

    //std::cout << "--Read " << col << ' ' << row << ' ' << x << ' ' << y << ' ' << z << std::endl;
    std::pair<float, float> P(col, row);
    if (xyz.find(P) == xyz.end()) {
      // std::cout << "--no average" << std::endl;
      xyz[P] = Vector3(x, y, z);
    }else{
      // std::cout << "--yes average" << std::endl;
      xyz[P] = (xyz[P] + Vector3(x, y, z))/2;
    }
  }
  
}

void handle_arguments(int argc, char *argv[], Options& opt) {
  po::options_description general_options("");
  general_options.add_options()
    ("output-prefix,o",  po::value(&opt.out_prefix),
     "Prefix for output filenames.")
    ("num-frequencies",  po::value(&opt.num_frequencies)->default_value(0),
     "Number of jitter frequencies.")
    ("ms-offset",  po::value(&opt.ms_offset)->default_value(0.0),
     "Number of lines of offset between the multispectral and pan cameras.")
    ("num-iterations",       po::value(&opt.num_iterations)->default_value(1000),
     "Set the maximum number of iterations.") 
    ("min-triangulation-angle",      po::value(&opt.min_triangulation_angle)->default_value(1e-8),
     "The minimum angle, in degrees, at which rays must meet at a triangulated point to accept this point as valid. It must be a positive value.")
    ("horiz-weight",      po::value(&opt.horiz_weight)->default_value(-1),
     "The horizontal weight to use to constrain the triangulated points to the initial values.")
    ("vert-weight",      po::value(&opt.vert_weight)->default_value(-1),
     "The vertical weight to use to constrain the triangulated points to the initial values.")
    ("images-to-use",  po::value(&opt.images_to_use)->default_value(""),
     "Out of the provided images, use 'pan,b7', 'pan,b8', or 'pan,b7,b8'.")
    ("mode",  po::value(&opt.mode)->default_value(""),
     "Optimize the following frequencies: 'x', 'y', 'z', 'xy', 'xz', 'yz', 'xyz'.")
    ("read-xyz-on-dem",   po::bool_switch(&opt.read_xyz_on_dem)->default_value(false),
     "Read externally computed xyz values on the DEM ray intersections should constrain to.")
    ("float-free-term",   po::bool_switch(&opt.float_free_term)->default_value(false),
     "Float the free term in the sines and cosines representation.")
    ("float-ms-offset",   po::bool_switch(&opt.float_ms_offset)->default_value(false),
     "Float the MS offset.");

  general_options.add(vw::cartography::GdalWriteOptionsDescription(opt));
  
  po::options_description positional("");
  positional.add_options()
    ("pan-image", po::value(&opt.pan_image))
    ("b7-image", po::value(&opt.b7_image))
    ("b8-image", po::value(&opt.b8_image))
    ("pan-camera", po::value(&opt.pan_camera))
    ("b7-camera", po::value(&opt.b7_camera))
    ("b8-camera", po::value(&opt.b8_camera));
  
  po::positional_options_description positional_desc;
  positional_desc.add("pan-image", 1);
  positional_desc.add("b7-image", 1);
  positional_desc.add("b8-image", 1);
  positional_desc.add("pan-camera", 1);
  positional_desc.add("b7-camera", 1);
  positional_desc.add("b8-camera", 1);
  
  std::string usage("[options] <pan-image> <b7-image> <b8-image> <pan-camera> <b7-camera> <b8-camera>");
  bool allow_unregistered = false;
  std::vector<std::string> unregistered;
  po::variables_map vm =
    asp::check_command_line(argc, argv, opt, general_options, general_options,
                            positional, positional_desc, usage,
                             allow_unregistered, unregistered );
  
  if ( !vm.count("pan-image")  || !vm.count("b7-image")  || !vm.count("b8-image")  ||
       !vm.count("pan-camera") || !vm.count("b7-camera") || !vm.count("b8-camera") )
    vw::vw_throw( vw::ArgumentErr() << "Not all inputs were specified.\n\n"
                  << usage << general_options );
  
  if (opt.out_prefix == "") 
    vw::vw_throw( vw::ArgumentErr() << "The output prefix is required.\n\n"
                  << usage << general_options );

  if (opt.num_frequencies <= 0) 
    vw::vw_throw( vw::ArgumentErr() << "Expecting a positive number of frequencies.\n\n"
                  << usage << general_options );

  if (opt.ms_offset <= 0) 
    vw::vw_throw( vw::ArgumentErr() << "Expecting a positive MS offset.\n\n"
                  << usage << general_options );

  // NOTE(oalexan1): The reason min_triangulation_angle cannot be 0 is deep inside
  // StereoModel.cc. Better keep it this way than make too many changes there.
  if ( opt.min_triangulation_angle <= 0.0 )
    vw_throw( ArgumentErr() << "The minimum triangulation angle must be positive.\n"
                            << usage << general_options );

  if ( opt.horiz_weight < 0.0 || opt.vert_weight < 0 )
    vw_throw( ArgumentErr() << "The horizontal and vertical weights must be non-negative.\n"
                            << usage << general_options );
  
  // Create the directory in which the output image will be written.
  vw::create_out_dir(opt.out_prefix);
}


CamPtr load_camera(Options const& opt, std::string const& image, std::string const& camera) {
  
  std::string session_str = "dg", out_prefix = "out";
  SessionPtr session(asp::StereoSessionFactory::create(session_str, // may change
                                                       opt,
                                                       image, camera, image, camera,
                                                       out_prefix));
  return session->camera_model(image, camera);
}

int main(int argc, char* argv[]) {
  
  Options opt;
  try {
    xercesc::XMLPlatformUtils::Initialize();

    handle_arguments(argc, argv, opt);

    std::cout << "pan image and camera: " << opt.pan_image << ' ' << opt.pan_camera << std::endl;
    std::cout << "b7 image and camera: " << opt.b7_image << ' ' << opt.b7_camera << std::endl;
    std::cout << "b8 image and camera: " << opt.b8_image << ' ' << opt.b8_camera << std::endl;

    std::cout << "-num frequencies is " << opt.num_frequencies << std::endl;
    
    // Load the cameras
    CamPtr pan_cam = load_camera(opt, opt.pan_image, opt.pan_camera);
    CamPtr b7_cam = load_camera(opt, opt.b7_image, opt.b7_camera);
    CamPtr b8_cam = load_camera(opt, opt.b8_image, opt.b8_camera);

//     std::cout << "--done loading!" << std::endl;
//     // TODO(oalexan1): Need to load just one camera and the B8 offset

//     std::cout << "--will cast!" << std::endl;
    asp::DGCameraModel * pan_dg_cam = dynamic_cast<asp::DGCameraModel*>(pan_cam.get());
    asp::DGCameraModel * b7_dg_cam = dynamic_cast<asp::DGCameraModel*>(b7_cam.get());
    asp::DGCameraModel * b8_dg_cam = dynamic_cast<asp::DGCameraModel*>(b8_cam.get());

    if (pan_dg_cam == NULL || b7_dg_cam == NULL || b8_dg_cam == NULL) 
    vw::vw_throw( vw::ArgumentErr() << "Expecting Digital Globe cameras.\n\n");

    pan_dg_cam->m_ms_offset = 0;
    b7_dg_cam->m_ms_offset = -opt.ms_offset;
    b8_dg_cam->m_ms_offset = opt.ms_offset;

    pan_dg_cam->m_coeffsx = std::vector<double>();
    pan_dg_cam->m_coeffsy = std::vector<double>();
    pan_dg_cam->m_coeffsz = std::vector<double>();
    
    b7_dg_cam->m_coeffsx = std::vector<double>();
    b7_dg_cam->m_coeffsy = std::vector<double>();
    b7_dg_cam->m_coeffsz = std::vector<double>();
    
    b8_dg_cam->m_coeffsx = std::vector<double>();
    b8_dg_cam->m_coeffsy = std::vector<double>();
    b8_dg_cam->m_coeffsz = std::vector<double>();

    std::cout << "ms_offset = " << pan_dg_cam->m_ms_offset << " for " << opt.pan_image<< std::endl;
    std::cout << "ms_offset = " << b7_dg_cam->m_ms_offset << " for " << opt.b7_image<< std::endl;
    std::cout << "ms_offset = " << b8_dg_cam->m_ms_offset << " for " << opt.b8_image<< std::endl;

    std::cout << "--horiz weight = " << opt.horiz_weight << std::endl;
    std::cout << "--vert weight = " << opt.vert_weight << std::endl;
    
    std::vector<double> coeffsx(2 * opt.num_frequencies + 1, 0);
    std::vector<double> coeffsy(2 * opt.num_frequencies + 1, 0);
    std::vector<double> coeffsz(2 * opt.num_frequencies + 1, 0);

    std::cout << "export ADJX='";
    for (int it = 0; it < int(coeffsx.size()) - 1; it++) std::cout << coeffsx[it] << " ";
    std::cout << coeffsx.back() << "'\n";
    
    std::cout << "export ADJY='";
    for (int it = 0; it < int(coeffsy.size()) - 1; it++) std::cout << coeffsy[it] << " ";
    std::cout << coeffsy.back() << "'\n";

    std::cout << "export ADJZ='";
    for (int it = 0; it < int(coeffsz.size()) - 1; it++) std::cout << coeffsz[it] << " ";
    std::cout << coeffsz.back() << "'\n";

    // 3 cameras and images
    std::vector<std::string> images;
    std::vector<CamPtr> cameras;
    std::vector<int> sign;
    std::map<std::pair<int, int>, std::string> match_files;
    
    images.push_back(opt.pan_image);
    cameras.push_back(pan_cam);
    sign.push_back(0);
    
    if (opt.images_to_use == "pan,b7" || opt.images_to_use == "pan,b7,b8") {
      images.push_back(opt.b7_image);
      cameras.push_back(b7_cam);
      sign.push_back(-1);
    }
    
    if (opt.images_to_use == "pan,b8" || opt.images_to_use == "pan,b7,b8") {
      images.push_back(opt.b8_image);
      cameras.push_back(b8_cam);
      sign.push_back(1);
    }

    if (cameras.empty()) 
      vw_throw( ArgumentErr() << "Could not load images and cameras. Check your parameters.\n" );
    
    std::map< std::pair<float, float>, Vector3> xyz_map;
    
    if (opt.images_to_use == "pan,b7") {
      std::string match_file = ip::match_filename(opt.out_prefix, opt.pan_image, opt.b7_image);
      match_files[std::pair<int, int>(0, 1)] = match_file;
      if (opt.read_xyz_on_dem)
        read_xyz(xyz_map, match_file + ".xyz");
    }
    
    if (opt.images_to_use == "pan,b8") {
      std::string match_file = ip::match_filename(opt.out_prefix, opt.pan_image, opt.b8_image);
      match_files[std::pair<int, int>(0, 1)] = match_file;
      if (opt.read_xyz_on_dem)
        read_xyz(xyz_map, match_file + ".xyz");
    }
    
    if (opt.images_to_use == "pan,b7,b8") {
      std::string match_file = ip::match_filename(opt.out_prefix, opt.pan_image, opt.b7_image);
      match_files[std::pair<int, int>(0, 1)] = match_file;
      if (opt.read_xyz_on_dem)
        read_xyz(xyz_map, match_file + ".xyz");

      match_file = ip::match_filename(opt.out_prefix, opt.pan_image, opt.b8_image);
      match_files[std::pair<int, int>(0, 2)] = match_file;
      if (opt.read_xyz_on_dem)
        read_xyz(xyz_map, match_file + ".xyz");
    }
    
    int min_matches = 30;
    vw::ba::ControlNetwork cnet("JitterSolve");
    bool triangulate_control_points = true;
    double forced_triangulation_distance = -1;
    bool success = build_control_network(triangulate_control_points,
                                         cnet, // output
                                         cameras, images, match_files,
                                         min_matches,
                                         opt.min_triangulation_angle*M_PI/180.0,
                                         forced_triangulation_distance);
    
    if (!success)
      vw_throw( ArgumentErr() << "Insufficient number of matches to solve for jitter.\n" );

    // Overwrite xyz to force them to be on the DEM
    if (opt.horiz_weight > 0 || opt.vert_weight > 0) {
      CRNJ crn_tmp;
      crn_tmp.read_controlnetwork(cnet);
      int num_cameras = crn_tmp.size();
      
      for (int icam = 0; icam < num_cameras; icam++) { // Camera loop
        for (crn_iter fiter = crn_tmp[icam].begin(); fiter != crn_tmp[icam].end(); fiter++){
          
          // The index of the 3D point this IP is for.
          int ipt = (**fiter).m_point_id;
          
          Vector2 observation = (**fiter).m_location;
          std::pair<float, float> obs(float(observation.x()), float(observation.y()));

          if (icam == 0) {
            // Index each xyz by the ip position in the pan camera (icam == 0)
            if (xyz_map.find(obs) == xyz_map.end()) 
              vw::vw_throw( vw::ArgumentErr() << "Expecting observation to be in the xyz map.\n\n");

            cnet[ipt].set_position(xyz_map[obs]);
          }
          
        }
      }
    }
    
    int num_points = cnet.size();
    std::cout << "Number of points in the network " << num_points << std::endl;

    std::vector<double> points_vec(num_points*NUM_POINT_PARAMS, 0.0);
    for (int ipt = 0; ipt < num_points; ipt++){
      for (int q = 0; q < NUM_POINT_PARAMS; q++){
        points_vec[ipt*NUM_POINT_PARAMS + q] = cnet[ipt].position()[q];
      }
    }
    double* points = &points_vec[0];
    
    CRNJ crn;
    crn.read_controlnetwork(cnet);

    std::cout << "--crn size is " << crn.size() << std::endl;
    int num_cameras = crn.size();

    double ms_offset = opt.ms_offset;
    std::cout << "export MS_OFFSET=" << ms_offset << std::endl;

    // The size of each vector that we optimize
    std::vector<int> block_sizes;
    block_sizes.push_back(coeffsx.size());
    block_sizes.push_back(coeffsy.size());
    block_sizes.push_back(coeffsz.size());
    block_sizes.push_back(NUM_POINT_PARAMS);
    block_sizes.push_back(1); // ms offset

    // The ceres problem
    ceres::Problem problem;
    
    // TODO(oalexan1): Get the datum from the stereo session
    vw::cartography::Datum datum("WGS84");
    std::cout << "--datum is " << datum << std::endl;

//     std::vector< std::vector<Vector3> > LLH(2);
//     std::vector< std::vector<Vector2> > Diff(2);

    std::vector<int> cam_residual_counts(num_cameras);

    for (int icam = 0; icam < num_cameras; icam++) { // Camera loop

      cam_residual_counts[icam] = 0;

      for (crn_iter fiter = crn[icam].begin(); fiter != crn[icam].end(); fiter++){ // IP loop

        // The index of the 3D point this IP is for.
        int ipt = (**fiter).m_point_id;

        VW_ASSERT(int(icam) < num_cameras,
                  ArgumentErr() << "Out of bounds in the number of cameras");
        VW_ASSERT(int(ipt)  < num_points,
                  ArgumentErr() << "Out of bounds in the number of points");

        // The observed value for the projection of point with index ipt into
        // the camera with index icam.
        Vector2 observation = (**fiter).m_location;

        // Note that xyz and point below must have the same value at this stage
        Vector3 xyz = cnet[ipt].position();

        // Each observation corresponds to a pair of a camera and a point
        double * point = points + ipt * NUM_POINT_PARAMS;

        Vector2 pix = cameras[icam]->point_to_pixel(xyz);
        Vector2 diff = pix - observation;

        // Vector3 llh = datum.cartesian_to_geodetic(xyz);

        // Need to add here the pointer to the given DG camera
        
        ceres::CostFunction* cost_function =
          ReprojectionError::Create(sign[icam], observation, cameras[icam], block_sizes);

        ceres::LossFunction* loss_function = get_jitter_loss_function();
        
        problem.AddResidualBlock(cost_function, loss_function,
                                 &coeffsx[0], &coeffsy[0], &coeffsz[0], point, &ms_offset);

        cam_residual_counts[icam] += 1; // Track the number of residual blocks for each camera
        
//         LLH[icam].push_back(llh);
//         Diff[icam].push_back(diff);
        
      } // end iterating over points
    } // end iterating over cameras

    // Add the horizontal and vertical xyz constraints
    if (opt.horiz_weight > 0 || opt.vert_weight > 0) {

      for (int ipt = 0; ipt < num_points; ipt++){
        
        Vector3 observation = cnet[ipt].position();
        double * point = points + ipt * NUM_POINT_PARAMS;
        
        ceres::CostFunction* horiz_cost_function
          = XYZHorizError::Create(observation, opt.horiz_weight);
        ceres::LossFunction* horiz_loss_function = get_jitter_loss_function();
        problem.AddResidualBlock(horiz_cost_function, horiz_loss_function, point);
        
        ceres::CostFunction* vert_cost_function
          = XYZVertError::Create(observation, opt.vert_weight);
        ceres::LossFunction* vert_loss_function = get_jitter_loss_function();
        problem.AddResidualBlock(vert_cost_function, vert_loss_function, point);
        
      }
    }
    
    std::string output_path = opt.out_prefix + "-initial_point_log.csv";
    save_residuals(opt, crn, cnet,
                   points_vec,
                   problem, datum, num_points, num_cameras,  
                   cam_residual_counts, output_path);
    
    vw::vw_out() << "Solving for jitter" << std::endl;
    Stopwatch sw;
    sw.start();
    // Solve the problem
    ceres::Solver::Options options;
    options.gradient_tolerance = 1e-16;
    options.function_tolerance = 1e-16;
    options.max_num_iterations = opt.num_iterations;
    options.max_num_consecutive_invalid_steps = 100; // try hard
    options.minimizer_progress_to_stdout = true;

    options.num_threads = opt.num_threads;
    
    options.linear_solver_type = ceres::SPARSE_SCHUR;
    //options.ordering_type = ceres::SCHUR;
    //options.eta = 1e-3; // FLAGS_eta;
    //options->max_solver_time_in_seconds = FLAGS_max_solver_time;
    //options->use_nonmonotonic_steps = FLAGS_nonmonotonic_steps;
    //if (FLAGS_line_search) {
    //  options->minimizer_type = ceres::LINE_SEARCH;
    //}

    if (opt.mode == "x") {
      if (!opt.float_free_term) {
        // Freeze the free term, coeffsx[0].
        problem.SetParameterization(&coeffsx[0],
                                    new ceres::SubsetParameterization(coeffsx.size(), {0}));
      }
      // Fix all y and z coeffs
      problem.SetParameterBlockConstant(&coeffsy[0]);
      problem.SetParameterBlockConstant(&coeffsz[0]);
    } else if (opt.mode == "y") {
      if (!opt.float_free_term) {
        // Freeze the free term, coeffsy[0].
        problem.SetParameterization(&coeffsy[0],
                                    new ceres::SubsetParameterization(coeffsy.size(), {0}));
      }
      // Fix all x and z coeffs
      problem.SetParameterBlockConstant(&coeffsx[0]);
      problem.SetParameterBlockConstant(&coeffsz[0]);
      
    } else if (opt.mode == "z") {
      if (!opt.float_free_term) {
        // Freeze the free term, coeffsz[0].
        problem.SetParameterization(&coeffsz[0],
                                    new ceres::SubsetParameterization(coeffsz.size(), {0}));
      }
      // Fix all x and y coeffs
      problem.SetParameterBlockConstant(&coeffsx[0]);
      problem.SetParameterBlockConstant(&coeffsy[0]);
    } else if (opt.mode == "xy") {
      if (!opt.float_free_term) {
        // Freeze the free term, coeffsx[0].
        problem.SetParameterization(&coeffsx[0],
                                    new ceres::SubsetParameterization(coeffsx.size(), {0}));
        // Freeze the free term, coeffsy[0].
        problem.SetParameterization(&coeffsy[0],
                                    new ceres::SubsetParameterization(coeffsy.size(), {0}));
      }
      // Fix all z coeffs
      problem.SetParameterBlockConstant(&coeffsz[0]);
    } else if (opt.mode == "xz") {
      if (!opt.float_free_term) {
        // Freeze the free term, coeffsx[0].
        problem.SetParameterization(&coeffsx[0],
                                    new ceres::SubsetParameterization(coeffsx.size(), {0}));
        // Freeze the free term, coeffsz[0].
        problem.SetParameterization(&coeffsz[0],
                                    new ceres::SubsetParameterization(coeffsz.size(), {0}));
      }
      // Fix all y coeffs
      problem.SetParameterBlockConstant(&coeffsy[0]);
    } else if (opt.mode == "yz") {
      if (!opt.float_free_term) {
        // Freeze the free term, coeffsy[0].
        problem.SetParameterization(&coeffsy[0],
                                    new ceres::SubsetParameterization(coeffsy.size(), {0}));
        // Freeze the free term, coeffsz[0].
        problem.SetParameterization(&coeffsz[0],
                                    new ceres::SubsetParameterization(coeffsz.size(), {0}));
      }
      // Fix all x coeffs
      problem.SetParameterBlockConstant(&coeffsx[0]);
    } else if (opt.mode == "xyz") {
      if (!opt.float_free_term) {
        // Freeze the free term, coeffsx[0].
        problem.SetParameterization(&coeffsx[0],
                                    new ceres::SubsetParameterization(coeffsx.size(), {0}));
        // Freeze the free term, coeffsy[0].
        problem.SetParameterization(&coeffsy[0],
                                    new ceres::SubsetParameterization(coeffsy.size(), {0}));
        
        // Freeze the free term, coeffsz[0].
        problem.SetParameterization(&coeffsz[0],
                                    new ceres::SubsetParameterization(coeffsz.size(), {0}));
      }
    } else {
      vw::vw_throw( vw::ArgumentErr() << "Unknown mode: " << opt.mode << ".\n\n");
    }

    if (!opt.float_ms_offset) 
      problem.SetParameterBlockConstant(&ms_offset);
        
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    vw_out() << summary.FullReport() << "\n";
    if (summary.termination_type == ceres::NO_CONVERGENCE){
      // Print a clarifying message, so the user does not think that the algorithm failed.
      vw_out() << "Found a valid solution, but did not reach the actual minimum." << std::endl;
    }

    output_path = opt.out_prefix + "-final_point_log.csv";
    save_residuals(opt, crn, cnet,
                   points_vec,
                   problem, datum, num_points, num_cameras,  
                   cam_residual_counts, output_path);
    
    std::cout << "export ADJX='";
    for (int it = 0; it < int(coeffsx.size()) - 1; it++) std::cout << coeffsx[it] << " ";
    std::cout << coeffsx.back() << "'\n";
    
    std::cout << "export ADJY='";
    for (int it = 0; it < int(coeffsy.size()) - 1; it++) std::cout << coeffsy[it] << " ";
    std::cout << coeffsy.back() << "'\n";

    std::cout << "export ADJZ='";
    for (int it = 0; it < int(coeffsz.size()) - 1; it++) std::cout << coeffsz[it] << " ";
    std::cout << coeffsz.back() << "'\n";

    std::cout << "export MS_OFFSET=" << ms_offset << std::endl;
    
    sw.stop();
    vw::vw_out() << "Jitter solve elapsed time: " << sw.elapsed_seconds() << std::endl;
    
    xercesc::XMLPlatformUtils::Terminate();
  } ASP_STANDARD_CATCHES;
}
