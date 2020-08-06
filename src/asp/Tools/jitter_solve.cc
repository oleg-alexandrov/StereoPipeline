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

/// A Ceres cost function. We pass in the observation and the model.
///  The result is the residual, the difference in the observation 
///  and the projection of the point into the camera, normalized by pixel_sigma.
struct ReprojectionError {
  ReprojectionError(Vector2 const& observation, Vector2 const& pixel_sigma,
                    std::vector<int> const& block_sizes):
    m_observation(observation),
    m_pixel_sigma(pixel_sigma),
    m_block_sizes(block_sizes) {}

  // Call to work with ceres::DynamicCostFunctions.
  // - Takes array of arrays.
  bool operator()(double const * const * parameters, double * residuals) const {

    // Need to create here the adjusted camera from the original camera
    
    try {
//       // Unpack the parameter blocks
//       std::vector<double const*> param_blocks(m_num_param_blocks);
//       for (size_t i = 0; i < m_num_param_blocks; i++) {
//         param_blocks[i] = parameters[i];
//       }

      // Use the camera model wrapper to handle all of the parameter blocks.
      Vector2 prediction; // = m_camera_wrapper->evaluate(param_blocks);

      //std::cout << "Got prediction " << prediction << std::endl;

      // The error is the difference between the predicted and observed position,
      // normalized by sigma.
      residuals[0] = (prediction[0] - m_observation[0])/m_pixel_sigma[0]; // Input units are pixels
      residuals[1] = (prediction[1] - m_observation[1])/m_pixel_sigma[1];

      residuals[0] = 0;
      residuals[1] = 0;
      
      //std::cout << "Got residuals " << residuals[0] << ", " << residuals[1] << std::endl;

    } catch (std::exception const& e) { // TODO: Catch only projection errors?
      // Failed to compute residuals

      //std::cout << "Caught exception!" << std::endl;

      residuals[0] = 1e+20;
      residuals[1] = 1e+20;
      return false;
    }
    
    return true;
  }


  // Factory to hide the construction of the CostFunction object from the client code.
  static ceres::CostFunction* Create(Vector2 const& observation, Vector2 const& pixel_sigma,
                                     std::vector<int> const& block_sizes){
    const int NUM_RESIDUALS = 2;

    ceres::DynamicNumericDiffCostFunction<ReprojectionError>* cost_function =
      new ceres::DynamicNumericDiffCostFunction<ReprojectionError>
      (new ReprojectionError(observation, pixel_sigma, block_sizes));
    
    // The residual size is always the same.
    cost_function->SetNumResiduals(NUM_RESIDUALS);

    // The camera wrapper knows all of the block sizes to add.
    for (size_t i = 0; i < block_sizes.size(); i++) {
      cost_function->AddParameterBlock(block_sizes[i]);
    }
    return cost_function;
  }

private:
  Vector2 m_observation;     ///< The pixel observation for this camera/point pair.
  Vector2 m_pixel_sigma;
  std::vector<int> m_block_sizes;
}; // End class ReprojectionError


struct Options : vw::cartography::GdalWriteOptions {
  std::string pan_image, ms7_image, ms8_image, pan_camera, ms7_camera, ms8_camera, out_prefix;
  int num_frequencies; // TODO(oalexan1): Use instead the frequency, in Hz
  int ms_offset;
  double min_triangulation_angle;
};

void handle_arguments(int argc, char *argv[], Options& opt) {
  po::options_description general_options("");
  general_options.add_options()
    ("output-prefix,o",  po::value(&opt.out_prefix),
     "Prefix for output filenames.")
    ("num-frequencies",  po::value(&opt.num_frequencies)->default_value(0),
     "Number of jitter frequencies.")
    ("ms-offset",  po::value(&opt.ms_offset)->default_value(0),
     "Number of lines of offset between the multispectral and pan cameras.")
    ("min-triangulation-angle",      po::value(&opt.min_triangulation_angle)->default_value(1e-8),
     "The minimum angle, in degrees, at which rays must meet at a triangulated point to accept this point as valid. It must be a positive value.");

  general_options.add(vw::cartography::GdalWriteOptionsDescription(opt));
  
  po::options_description positional("");
  positional.add_options()
    ("pan-image", po::value(&opt.pan_image))
    ("ms7-image", po::value(&opt.ms7_image))
    ("ms8-image", po::value(&opt.ms8_image))
    ("pan-camera", po::value(&opt.pan_camera))
    ("ms7-camera", po::value(&opt.ms7_camera))
    ("ms8-camera", po::value(&opt.ms8_camera));
  
  po::positional_options_description positional_desc;
  positional_desc.add("pan-image", 1);
  positional_desc.add("ms7-image", 1);
  positional_desc.add("ms8-image", 1);
  positional_desc.add("pan-camera", 1);
  positional_desc.add("ms7-camera", 1);
  positional_desc.add("ms8-camera", 1);
  
  std::string usage("[options] <pan-image> <ms7-image> <ms8-image> <pan-camera> <ms7-camera> <ms8-camera>");
  bool allow_unregistered = false;
  std::vector<std::string> unregistered;
  po::variables_map vm =
    asp::check_command_line(argc, argv, opt, general_options, general_options,
                            positional, positional_desc, usage,
                             allow_unregistered, unregistered );
  
  if ( !vm.count("pan-image")  || !vm.count("ms7-image")  || !vm.count("ms8-image")  ||
       !vm.count("pan-camera") || !vm.count("ms7-camera") || !vm.count("ms8-camera") )
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

// TODO(oalexan1): Expose the threshold
ceres::LossFunction* get_jitter_loss_function(){
  return new ceres::CauchyLoss(0.5);
}

int main(int argc, char* argv[]) {
  
  Options opt;
  try {
    xercesc::XMLPlatformUtils::Initialize();

    handle_arguments(argc, argv, opt);

    std::cout << "pan image and camera: " << opt.pan_image << ' ' << opt.pan_camera << std::endl;
    std::cout << "ms7 image and camera: " << opt.ms7_image << ' ' << opt.ms7_camera << std::endl;
    std::cout << "ms8 image and camera: " << opt.ms8_image << ' ' << opt.ms8_camera << std::endl;

    std::cout << "-num frequencies is " << opt.num_frequencies << std::endl;
    std::cout << "-ms offset is " << opt.ms_offset << std::endl;
    
    // Load the cameras
    CamPtr pan_cam = load_camera(opt, opt.pan_image, opt.pan_camera);
    CamPtr ms7_cam = load_camera(opt, opt.ms7_image, opt.ms7_camera);
    CamPtr ms8_cam = load_camera(opt, opt.ms8_image, opt.ms8_camera);

    std::cout << "--done loading!" << std::endl;
    // TODO(oalexan1): Need to load just one camera and the MS8 offset

    std::cout << "--will cast!" << std::endl;
    asp::DGCameraModel * pan_dg_cam = dynamic_cast<asp::DGCameraModel*>(pan_cam.get());
    asp::DGCameraModel * ms7_dg_cam = dynamic_cast<asp::DGCameraModel*>(ms7_cam.get());
    asp::DGCameraModel * ms8_dg_cam = dynamic_cast<asp::DGCameraModel*>(ms8_cam.get());

    if (pan_dg_cam == NULL || ms7_dg_cam == NULL || ms8_dg_cam == NULL) 
    vw::vw_throw( vw::ArgumentErr() << "Expecting Digital Globe cameras.\n\n");

    pan_dg_cam->m_ms_offset = 0;
    ms7_dg_cam->m_ms_offset = -opt.ms_offset;
    ms8_dg_cam->m_ms_offset = opt.ms_offset;
    

#if 1
    int mode = atoi(getenv("MODE"));
    std::cout << "--mode is " << mode << std::endl;
    asp::LinescanModelFreq f_pan_cam(pan_dg_cam, mode);
    asp::LinescanModelFreq f_ms7_cam(ms7_dg_cam, mode);
    asp::LinescanModelFreq f_ms8_cam(ms8_dg_cam, mode);
    int num = 10;
    std::vector<double> coeffs(2*num+1, 0);
    coeffs.back() = 1e-6; // temporary!
    f_pan_cam.m_coeffs = coeffs;
    f_ms7_cam.m_coeffs = coeffs;
    f_ms8_cam.m_coeffs = coeffs;
#endif
    
//     int num =
#if 0
    int mode = 1;
    asp::LinescanModelFreq f_ms8_cam(ms8_dg_cam, mode);

    for (int col = 0; col < 100; col++) {
      for (int row = 0; row < 100; row++) {
        Vector2 pix(col, row);
        Vector3 ctr1 = ms8_dg_cam->camera_center(pix);
        Vector3 ctr2 = f_ms8_cam.camera_center(pix);
        std::cout << "-diff1 is " << ctr1 << ' ' << ctr2 << ' ' << norm_2(ctr1-ctr2) << std::endl;

        Vector3 dir1 = ms8_dg_cam->pixel_to_vector(pix);
        Vector3 dir2 = f_ms8_cam.pixel_to_vector(pix);
        std::cout << "-diff2 is " << dir1 << ' ' << dir2 << ' ' << norm_2(dir1-dir2) << std::endl;

        Vector3 xyz = ctr1 + 600000 * dir1;
        Vector2 pix1 = ms8_dg_cam->point_to_pixel(xyz);
        Vector2 pix2 = f_ms8_cam.point_to_pixel(xyz);
        std::cout << "--diff3 is " << pix1 << ' ' << pix2 << ' ' << norm_2(pix1 - pix2) << std::endl;
      }
    }
#endif
#if 0
    // Two cameras
    std::string match_filename = ip::match_filename(opt.out_prefix, opt.pan_image, opt.ms8_image);

    std::cout << "--reading binary match file " << match_filename << std::endl;
    
    std::vector<ip::InterestPoint> ip1, ip2;
    ip::read_binary_match_file(match_filename, ip1, ip2);
    
    std::cout << "--now here!" << std::endl;
    std::cout << "number of ip is " << ip1.size() << ' ' << ip2.size() << std::endl;

    bool least_squares = false;
    vw::stereo::StereoModel sm(pan_dg_cam, ms8_dg_cam, least_squares,
                               opt.min_triangulation_angle*M_PI/180.0 );

    std::string res_file = opt.out_prefix + "-residuals.csv";
    std::cout << "Writing: " << res_file << std::endl;
    std::ofstream ofs(res_file.c_str());

    // TODO(oalexan1): Get the datum from the stereo session
    vw::cartography::Datum datum("WGS84");
    std::cout << "--datum is " << datum << std::endl;
    
    ofs << "# lon, lat, height_above_datum, mean_residual, num_observations, indiv residuals\n";
    ofs.precision(18);
    
    for (size_t it = 0; it < ip1.size(); it++) {
      Vector3 err;
      Vector2 p1 = Vector2(ip1[it].x, ip1[it].y);
      Vector2 p2 = Vector2(ip2[it].x, ip2[it].y);
      Vector3 xyz = sm(p1, p2, err);
      
      Vector2 pix1 = pan_dg_cam->point_to_pixel(xyz);
      Vector2 pix2 = ms8_dg_cam->point_to_pixel(xyz);

      Vector2 diff1 = pix1 - p1;
      Vector2 diff2 = pix2 - p2;

      double mean_res = (std::abs(diff1[0]) + std::abs(diff1[1])
                         + std::abs(diff2[0]) + std::abs(diff2[1]))/4.0;
      
      Vector3 llh = datum.cartesian_to_geodetic(xyz);

      ofs << llh[0] << ", " << llh[1] << ", " << llh[2] << ", "
          << mean_res << ", " << 2 << ", "
          << diff1.x() << ", " << diff1.y() << ", "
          << diff2.x() << ", " << diff2.y() << std::endl;
    }
    ofs.close();
#endif

    // 3 cameras and images
    std::vector<std::string> images;
    images.push_back(opt.pan_image);
    //images.push_back(opt.ms7_image); // temporary!
    images.push_back(opt.ms8_image);
    std::vector<CamPtr> cameras;
    cameras.push_back(pan_cam);
    //cameras.push_back(ms7_cam); // temporary!
    cameras.push_back(ms8_cam);

    // Temporary!
    std::map< std::pair<int, int>, std::string> match_files;
    match_files[std::pair<int, int>(0, 1)]
      = ip::match_filename(opt.out_prefix, opt.pan_image, opt.ms8_image);
    // temporary!
//     match_files[std::pair<int, int>(0, 1)]
//       = ip::match_filename(opt.out_prefix, opt.pan_image, opt.ms7_image);
//     match_files[std::pair<int, int>(0, 2)]
//       = ip::match_filename(opt.out_prefix, opt.pan_image, opt.ms8_image);
    
    int min_matches = 30;
    ba::ControlNetwork cnet("JitterSolve");
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

    // The size of each vector that we optimize
    std::vector<int> block_sizes;
    block_sizes.push_back(coeffs.size());
    block_sizes.push_back(NUM_POINT_PARAMS);
    
    // The ceres problem
    ceres::Problem problem;
    
    // TODO(oalexan1): Get the datum from the stereo session
    vw::cartography::Datum datum("WGS84");
    std::cout << "--datum is " << datum << std::endl;

    std::vector< std::vector<Vector3> > LLH(2);
    std::vector< std::vector<Vector2> > Diff(2);
    
    for (int icam = 0; icam < num_cameras; icam++) { // Camera loop
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
        Vector2 pixel_sigma = (**fiter).m_scale;
        
        //std::cout << "--pixel sigma is " << pixel_sigma << std::endl;
        
        // This is a bugfix
        if (pixel_sigma != pixel_sigma) // nan check
          pixel_sigma = Vector2(1, 1);

        Vector3 xyz = cnet[ipt].position();

        Vector2 pix = cameras[icam]->point_to_pixel(xyz);
        Vector2 diff = pix - observation;

        Vector3 llh = datum.cartesian_to_geodetic(xyz);

        // Each observation corresponds to a pair of a camera and a point
        double * point = points + ipt * NUM_POINT_PARAMS;

        // Need to add here the pointer to the given DG camera
        
        ceres::CostFunction* cost_function =
          ReprojectionError::Create(observation, pixel_sigma, block_sizes);

        ceres::LossFunction* loss_function = get_jitter_loss_function();
        
        problem.AddResidualBlock(cost_function, loss_function, &coeffs[0], point);
        
        LLH[icam].push_back(llh);
        Diff[icam].push_back(diff);
        
        //         // Call function to add the appropriate Ceres residual block.
        //         add_reprojection_residual_block(observation, pixel_sigma, ipt, icam,
        //                                         is_gcp, param_storage, opt, problem);
        
      } // end iterating over points
    } // end iterating over cameras
    
    std::string res_file = opt.out_prefix + "-residuals.csv";
    std::cout << "Writing: " << res_file << std::endl;
    std::ofstream ofs(res_file.c_str());
    ofs << "# lon, lat, height_above_datum, mean_residual, num_observations, indiv residuals\n";
    ofs.precision(18);
    
    for (size_t it = 0; it < Diff[0].size(); it++) {
      
      Vector3 llh = LLH[0][it];
      Vector2 diff1 = Diff[0][it];
      Vector2 diff2 = Diff[1][it];

      double mean_res = (std::abs(diff1[0]) + std::abs(diff1[1])
                         + std::abs(diff2[0]) + std::abs(diff2[1]))/4.0;
      
      ofs << llh[0] << ", " << llh[1] << ", " << llh[2] << ", "
          << mean_res << ", " << 2 << ", "
          << diff1.x() << ", " << diff1.y() << ", "
          << diff2.x() << ", " << diff2.y() << std::endl;
    }
    ofs.close();

    vw::vw_out() << "Solving for jitter" << std::endl;
    Stopwatch sw;
    sw.start();
    // Solve the problem
    ceres::Solver::Options options;
    options.gradient_tolerance = 1e-16;
    options.function_tolerance = 1e-16;
    options.max_num_iterations = 1000; // TODO(oalexan1): Expose this
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
    
    // Use a callback function at every iteration.
    //   PiecewiseBaCallback callback;
    //   options.callbacks.push_back(&callback);
    //   options.update_state_every_iteration = true; // ensure we have the latest adjustments
    
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    vw_out() << summary.FullReport() << "\n";
    if (summary.termination_type == ceres::NO_CONVERGENCE){
      // Print a clarifying message, so the user does not think that the algorithm failed.
      vw_out() << "Found a valid solution, but did not reach the actual minimum." << std::endl;
    }
    
    sw.stop();
    vw::vw_out() << "Jitter solve elapsed time: " << sw.elapsed_seconds() << std::endl;
    
    xercesc::XMLPlatformUtils::Terminate();
  } ASP_STANDARD_CATCHES;
}
