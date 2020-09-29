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

#include <asp/Core/Macros.h>
#include <asp/Sessions/StereoSession.h>
#include <asp/Sessions/StereoSessionFactory.h>
#include <asp/Core/StereoSettings.h>
#include <vw/Cartography/CameraBBox.h>
#include <vw/InterestPoint/InterestData.h>
#include <vw/InterestPoint/Matcher.h>
#include <xercesc/util/PlatformUtils.hpp>
#include <asp/Camera/LinescanDGModel.h>

namespace po = boost::program_options;
namespace fs = boost::filesystem;
using namespace vw;

struct Options : vw::cartography::GdalWriteOptions {
  std::string pan_image, pan_camera, ms_camera, bundle_adjust_prefix, dem, pan_to_ms_disp,
    ms_image, out_prefix, band;
  int num_matches_x, num_matches_y;
  bool synthetic_jitter;
  double ms_offset;
};

void handle_arguments(int argc, char *argv[], Options& opt) {
  
  po::options_description general_options("");
  general_options.add_options()
    ("dem",  po::value(&opt.dem)->default_value(""),
     "The DEM to use to undo the processing that was applied to the MS image.")
    ("pan-to-ms-disp",  po::value(&opt.pan_to_ms_disp)->default_value(""),
     "The disparity from the PAN image to the processed MS image.")
    ("num-matches-x",  po::value(&opt.num_matches_x)->default_value(0),
     "How many matches to find between the MS and PAN image in the x direction.")
    ("num-matches-y",  po::value(&opt.num_matches_y)->default_value(0),
     "How many matches to find between the MS and PAN image in the y direction.")
    ("ms-offset",  po::value(&opt.ms_offset)->default_value(0.0),
     "Number of lines of offset between the multispectral and pan cameras.")
    ("synthetic-jitter",   po::bool_switch(&opt.synthetic_jitter)->default_value(false),
     "Form synthetic jitter by using ADJX, ADJY and ADJZ env vars rather than from disparity.")
    ("band",  po::value(&opt.band)->default_value(""),
     "The band for the image, b7 (NIR1) or b8 (NIR2).")
    ("output-prefix,o",  po::value(&opt.out_prefix),
     "Prefix for output filenames.");

  general_options.add(vw::cartography::GdalWriteOptionsDescription(opt));
  
  po::options_description positional("");
  positional.add_options()
    ("pan-image", po::value(&opt.pan_image))
    ("proc-ms-image", po::value(&opt.ms_image))
    ("pan-camera", po::value(&opt.pan_camera))
    ("ms-camera", po::value(&opt.ms_camera));
  
  po::positional_options_description positional_desc;
  positional_desc.add("pan-image", 1);
  positional_desc.add("proc-ms-image", 1);
  positional_desc.add("pan-camera", 1);
  positional_desc.add("ms-camera", 1);
  
  std::string usage("[options] <pan-image> <proc-ms-image> <pan-camera> <ms-camera>");
  bool allow_unregistered = false;
  std::vector<std::string> unregistered;
  po::variables_map vm =
    asp::check_command_line(argc, argv, opt, general_options, general_options,
                            positional, positional_desc, usage,
                             allow_unregistered, unregistered);
  
  if ( !vm.count("pan-image") || !vm.count("proc-ms-image") ||
       !vm.count("pan-camera") || !vm.count("ms-camera") )
    vw::vw_throw( vw::ArgumentErr() << "Not all inputs were specified.\n\n"
                  << usage << general_options);

//   if (opt.ms_offset <= 0) {
//     vw::vw_throw( vw::ArgumentErr() << "The MS offset must be positive.\n\n"
//                   << usage << general_options);
//   }
  
  if (opt.dem == "") 
    vw::vw_throw( vw::ArgumentErr() << "An input DEM is required.\n\n"
                  << usage << general_options);

  if (opt.pan_to_ms_disp == "") 
    vw::vw_throw( vw::ArgumentErr() << "The processed MS to PAN image disparity is required.\n\n"
                  << usage << general_options);

  if (opt.num_matches_x <= 0 || opt.num_matches_y <= 0) 
    vw::vw_throw( vw::ArgumentErr() << "Number of matches to compute must be positive.\n\n"
                  << usage << general_options);

  if (opt.out_prefix == "") 
    vw::vw_throw( vw::ArgumentErr() << "The output prefix is required.\n\n"
                  << usage << general_options);

  // Create the directory in which the output image will be written.
  vw::create_out_dir(opt.out_prefix);
  
}

// Take an MS pixel, travel to the DEM representing the ground,
// and then project into the PAN camera. This will give the pixel
// as processed by Digital Globe when they create synthetic (processed)
// MS images.
// TODO(oalexan1): Sometimes this function fails when it should not!
// Likely the code doing the intersection with the ground is not robust.
bool raw_to_proc_ms_pix(Vector2 const& ms_pix,
                        ImageViewRef< PixelMask<float> > dem,
                        vw::cartography::GeoReference const& dem_georef,
                        boost::shared_ptr<vw::camera::CameraModel> pan_cam, 
                        boost::shared_ptr<vw::camera::CameraModel> ms_cam,
                        Vector2 & pan_pix) {
  
  // Project from the MS camera into the ground
  bool   has_intersection     = false;
  bool   treat_nodata_as_zero = false;
  double height_error_tol     = 1e-5;   // error in DEM height
  double max_abs_tol          = 1e-14;  // abs cost function change b/w iterations
  double max_rel_tol          = 1e-14;
  int    num_max_iter         = 100;
  Vector3 xyz_guess           = Vector3();
  Vector3 camera_ctr          = ms_cam->camera_center(ms_pix);  // Get ray from this pixel
  Vector3 camera_vec          = ms_cam->pixel_to_vector(ms_pix);
         
  // Use iterative solver call to compute an intersection of the pixel with the DEM	
  Vector3 xyz
    = vw::cartography::camera_pixel_to_dem_xyz(camera_ctr, camera_vec,
                                               dem, dem_georef,
                                               treat_nodata_as_zero,
                                               has_intersection,
                                               height_error_tol, max_abs_tol, max_rel_tol,
                                               num_max_iter, xyz_guess);
         
  // Quit if we did not find an intersection
  if (!has_intersection)
    return false;

  pan_pix = pan_cam->point_to_pixel(xyz);
  
  return true;
}

// Take an MS pixel in the processed image, travel to the DEM representing the ground,
// via the PAN cam, and then project into the MS camera. This will give the raw
// MS pixel.
// TODO(oalexan1): Sometimes this function fails when it should not!
// Likely the code doing the intersection with the ground is not robust.
bool proc_ms_to_raw_ms_pix(Vector2 const& proc_ms_pix,
                           ImageViewRef< PixelMask<float> > dem,
                           vw::cartography::GeoReference const& dem_georef,
                           boost::shared_ptr<vw::camera::CameraModel> pan_cam, 
                           boost::shared_ptr<vw::camera::CameraModel> ms_cam,
                           Vector2 & raw_ms_pix, Vector3 & xyz) {
  
  // Project from the MS camera into the ground
  bool   has_intersection     = false;
  bool   treat_nodata_as_zero = false;
  double height_error_tol     = 1e-5;   // error in DEM height
  double max_abs_tol          = 1e-14;  // abs cost function change b/w iterations
  double max_rel_tol          = 1e-14;
  int    num_max_iter         = 100;
  Vector3 xyz_guess           = Vector3();
  Vector3 camera_ctr          = pan_cam->camera_center(proc_ms_pix);  // Get ray from this pixel
  Vector3 camera_vec          = pan_cam->pixel_to_vector(proc_ms_pix);
         
  // Use iterative solver call to compute an intersection of the pixel with the DEM	
  xyz = vw::cartography::camera_pixel_to_dem_xyz(camera_ctr, camera_vec,
                                                 dem, dem_georef,
                                                 treat_nodata_as_zero,
                                                 has_intersection,
                                                 height_error_tol, max_abs_tol, max_rel_tol,
                                                 num_max_iter, xyz_guess);
  
  // Quit if we did not find an intersection
  if (!has_intersection)
    return false;

  raw_ms_pix = ms_cam->point_to_pixel(xyz);
  
  return true;
}

typedef boost::scoped_ptr<asp::StereoSession> SessionPtr;

int main(int argc, char* argv[]) {

  vw::cartography::Datum datum("WGS84");
  
  Options opt;
  try {
    xercesc::XMLPlatformUtils::Initialize();

    handle_arguments(argc, argv, opt);

    // Initialize the session
    std::string session_str = "dg";
    std::string out_prefix = "out";
    SessionPtr session(asp::StereoSessionFactory::create(session_str, // may change
                                                         opt,
                                                         opt.pan_image,
                                                         opt.ms_image,
                                                         opt.pan_camera,
                                                         opt.ms_camera,
                                                         out_prefix));

    // Load the pan camera
    boost::shared_ptr<vw::camera::CameraModel> pan_cam = session->camera_model(opt.pan_image,
                                                                           opt.pan_camera);

    // Load the ms camera
    boost::shared_ptr<vw::camera::CameraModel> ms_cam = session->camera_model(opt.ms_image,
                                                                           opt.ms_camera);

    asp::DGCameraModel * pan_dg_cam = dynamic_cast<asp::DGCameraModel*>(pan_cam.get());
    asp::DGCameraModel * ms_dg_cam = dynamic_cast<asp::DGCameraModel*>(ms_cam.get());
    
    if (pan_dg_cam == NULL || ms_dg_cam == NULL) 
      vw::vw_throw( vw::ArgumentErr() << "Expecting Digital Globe cameras.\n\n");

    pan_dg_cam->m_ms_offset = 0;
    if (opt.band == "b8") {
      ms_dg_cam->m_ms_offset  = opt.ms_offset;
    } else if (opt.band == "b7") {
      ms_dg_cam->m_ms_offset  = -opt.ms_offset;
    }else{
      vw::vw_throw( vw::ArgumentErr() << "Invalid band: " << opt.band);
    }

    std::cout << "Setting PAN ms offset to " << pan_dg_cam->m_ms_offset << std::endl;
    std::cout << "Setting MS ms offset to "  << ms_dg_cam->m_ms_offset  << std::endl;
    
    // Read the image
    //std::cout << "Reading image: " << opt.pan_image << std::endl;
    bool has_image_nodata = true; // we will surely have no-data on output, even if not on input
    float image_nodata = -std::numeric_limits<float>::max();
    boost::shared_ptr<DiskImageResource> img_rsrc(new DiskImageResourceGDAL(opt.pan_image));
    if (img_rsrc->has_nodata_read())
      image_nodata = img_rsrc->nodata_read();
    ImageViewRef< PixelMask<float> > pan_image
      (create_mask(DiskImageView<float>(opt.pan_image), image_nodata));
    
    // Read the DEM
    //std::cout << "dem file is " << opt.dem << std::endl;
    float dem_nodata = -std::numeric_limits<float>::max();
    if (vw::read_nodata_val(opt.dem, dem_nodata))
      vw_out() << "Dem nodata: " << dem_nodata << std::endl;
    ImageViewRef< PixelMask<float> > dem(create_mask(DiskImageView<float>(opt.dem),
                                                     dem_nodata));
    //std::cout << "cols and rows are " << dem.cols() << ' ' << dem.rows() << std::endl;
    vw::cartography::GeoReference dem_georef;
    if (!read_georeference(dem_georef, opt.dem))
      vw_throw( ArgumentErr() << "The input DEM " << opt.dem << " has no georeference.\n");

    //std::cout << "--Reading " << opt.pan_to_ms_disp << std::endl;
    DiskImageView< PixelMask<Vector2f> > disparity(opt.pan_to_ms_disp);

    if (disparity.cols() != pan_image.cols() || disparity.rows() != pan_image.rows()) 
      vw_throw( ArgumentErr() << "The input image and disparity must have the same dimensions.\n");
      
    ImageViewRef< PixelMask<Vector2f> > interp_disp  
      = interpolate(disparity, BilinearInterpolation(), ConstantEdgeExtension());
        
    int num_cols = pan_image.cols();
    int num_rows = pan_image.rows();

    double ratio_x = double(opt.num_matches_x)/num_cols;
    double ratio_y = double(opt.num_matches_y)/num_rows;

    int num_sample_cols = std::min(std::max(1, int(round(num_cols * ratio_x))), num_cols);
    int num_sample_rows = std::min(std::max(1, int(round(num_rows * ratio_y))), num_rows);

//     std::cout << "num sample cols and rows " << num_sample_cols << ' ' << num_sample_rows
//               << std::endl;
//     std::cout << "--ratio x and ratio y " << ratio_x << ' ' << ratio_y << std::endl;

    std::vector< std::vector<double> > xyz_list;

    ImageView<float> dispx;
    ImageView<float> dispy;
    dispx.set_size(num_sample_cols, num_sample_rows);
    dispy.set_size(num_sample_cols, num_sample_rows);
    
    std::vector<double> coeffsx, coeffsy, coeffsz;

    if (opt.synthetic_jitter) {
      asp::parse_coeffs(coeffsx, "ADJX", "Multi");
      asp::parse_coeffs(coeffsy, "ADJY", "Multi");
      asp::parse_coeffs(coeffsz, "ADJZ", "Multi");
    }else{
      pan_dg_cam->m_coeffsx = std::vector<double>();
      pan_dg_cam->m_coeffsy = std::vector<double>();
      pan_dg_cam->m_coeffsz = std::vector<double>();
      
      ms_dg_cam->m_coeffsx = std::vector<double>();
      ms_dg_cam->m_coeffsy = std::vector<double>();
      ms_dg_cam->m_coeffsz = std::vector<double>();
    }

    // Matches from PAN image to raw MS image
    std::vector<vw::ip::InterestPoint> left_ip, right_ip;

    for (int samp_col = 0; samp_col < num_sample_cols; samp_col++) {
      for (int samp_row = 0; samp_row < num_sample_rows; samp_row++) {

        dispx(samp_col, samp_row) = 0;
        dispy(samp_col, samp_row) = 0;
        
        //std::cout << std::endl;

        Vector2 pan_pix = Vector2(samp_col/ratio_x, samp_row/ratio_y);

        // Add jitter to cameras
        if (opt.synthetic_jitter) {
          pan_dg_cam->m_coeffsx = coeffsx;
          pan_dg_cam->m_coeffsy = coeffsy;
          pan_dg_cam->m_coeffsz = coeffsz;
          
          ms_dg_cam->m_coeffsx = coeffsx;
          ms_dg_cam->m_coeffsy = coeffsy;
          ms_dg_cam->m_coeffsz = coeffsz;
        }

        // The MS pixel when projected in the PAN camera
        Vector2 tweaked_pix;
        if (opt.synthetic_jitter) {
          tweaked_pix = pan_pix;
        }else{
          // Get it from disparity
          PixelMask<Vector2f> disp_val = interp_disp(pan_pix.x(), pan_pix.y());
          if (!is_valid(disp_val)) 
            continue;
          tweaked_pix = pan_pix + disp_val.child();
        }

#if 0
        if (samp_col == 2 && samp_row == 54) {
          std::cout << "--temporary!" << std::endl;
          std::cout << "samp col and row is " << samp_col << ' ' << samp_row << std::endl;
          tweaked_pix = pan_pix;
          pan_dg_cam->m_coeffsx = std::vector<double>();
          pan_dg_cam->m_coeffsy = std::vector<double>();
          pan_dg_cam->m_coeffsz = std::vector<double>();
          
          ms_dg_cam->m_coeffsx = std::vector<double>();
          ms_dg_cam->m_coeffsy = std::vector<double>();
          ms_dg_cam->m_coeffsz = std::vector<double>();
        }
#endif
        
        
        // Project from the PAN camera into the ground
        bool   has_intersection     = false;
        bool   treat_nodata_as_zero = false;
        double height_error_tol     = 1e-5;   // error in DEM height
        double max_abs_tol          = 1e-14;  // abs cost function change b/w iterations
        double max_rel_tol          = 1e-14;
        int    num_max_iter         = 100;
        Vector3 xyz_guess           = Vector3();
         
        // Use iterative solver call to compute an intersection of the pixel with the DEM	
        Vector3 camera_ctr          = pan_dg_cam->camera_center(tweaked_pix); 
        Vector3 camera_vec          = pan_dg_cam->pixel_to_vector(tweaked_pix);
        Vector3 xyz
          = vw::cartography::camera_pixel_to_dem_xyz(camera_ctr, camera_vec,
                                                     dem, dem_georef,
                                                     treat_nodata_as_zero,
                                                     has_intersection,
                                                     height_error_tol, max_abs_tol, max_rel_tol,
                                                     num_max_iter, xyz_guess);
         
        // Quit if we did not find an intersection
        if (!has_intersection)
          continue;
        
        // std::cout << "-xyz is " << xyz << std::endl;
        // Vector2 ms_pix = m_ms_dg_cam->point_to_pixel(xyz);
        // std::cout << "--ms pix " << pix << ' ' << ms_pix << ' ' << norm_2(pix - ms_pix)
        // << std::endl;

        // The MS pixel, now in the MS camera
        Vector2 ms_pix = ms_dg_cam->point_to_pixel(xyz);

#if 0
        if (samp_col == 2 && samp_row == 54) {
          std::cout << std::endl;
          std::cout << std::endl;
          std::cout << "--x" << std::endl;
        
          double min_triangulation_angle = 1e-10;
          double angle_tol =
            vw::stereo::StereoModel::robust_1_minus_cos(min_triangulation_angle*M_PI/180);
          bool use_least_squares = false;
          //std::vector<const vw::camera::CameraModel *> camera_ptrs;
          //camera_ptrs.push_back(pan_dg_cam);
          //camera_ptrs.push_back(ms_dg_cam);
          //vw::stereo::StereoModel stereo_model(camera_ptrs, use_least_squares, angle_tol);

          vw::stereo::StereoModel sm(pan_dg_cam, ms_dg_cam, use_least_squares, angle_tol);
          Vector3 err;
          Vector3 xyz3 = sm(pan_pix, ms_pix, err);
          std::cout.precision(18);
          std::cout << "xyz is " << xyz << ' ' << xyz3 << ' ' << norm_2(xyz - xyz3) << std::endl;

          Vector2 ms_pix0 = ms_dg_cam->point_to_pixel(xyz3);
          Vector2 pan_pix0 = pan_dg_cam->point_to_pixel(xyz3);

          std::cout << "--difff 0 " << pan_pix0 << ' '
                    << ms_pix0 - ms_pix << ' '
                    << pan_pix0 - pan_pix << std::endl;

        
          pan_dg_cam->m_coeffsx = coeffsx;
          pan_dg_cam->m_coeffsy = coeffsy;
          pan_dg_cam->m_coeffsz = coeffsz;
          
          ms_dg_cam->m_coeffsx = coeffsx;
          ms_dg_cam->m_coeffsy = coeffsy;
          ms_dg_cam->m_coeffsz = coeffsz;
        
          ms_pix0 = ms_dg_cam->point_to_pixel(xyz3);
          pan_pix0 = pan_dg_cam->point_to_pixel(xyz3);

          std::cout << "--difff 1 " << pan_pix0 << ' '
                    << ms_pix0 - ms_pix << ' '
                    << pan_pix0 - pan_pix << std::endl;
          std::cout << std::endl;
          std::cout << std::endl;

          Vector3 xyz4 = sm(pan_pix, ms_pix, err);
          std::cout << "xyz4 is " << xyz << ' ' << xyz4 << ' ' << norm_2(xyz - xyz4) << std::endl;

          ms_pix0 = ms_dg_cam->point_to_pixel(xyz4);
          pan_pix0 = pan_dg_cam->point_to_pixel(xyz4);

          std::cout << "--difff 2 " << pan_pix0 << ' '
                    << ms_pix0 - ms_pix << ' '
                    << pan_pix0 - pan_pix << std::endl;
          std::cout << std::endl;
          std::cout << std::endl;

          
        }
#endif
        
        // Now compute the disparity due to the jitter
        // So we use pixels measured with jitter
        // but cameras without jitter.
        pan_dg_cam->m_coeffsx = std::vector<double>();
        pan_dg_cam->m_coeffsy = std::vector<double>();
        pan_dg_cam->m_coeffsz = std::vector<double>();
        
        ms_dg_cam->m_coeffsx = std::vector<double>();
        ms_dg_cam->m_coeffsy = std::vector<double>();
        ms_dg_cam->m_coeffsz = std::vector<double>();
        
         
        // Use iterative solver call to compute an intersection of the pixel with the DEM	
        camera_ctr = ms_dg_cam->camera_center(ms_pix); 
        camera_vec = ms_dg_cam->pixel_to_vector(ms_pix);
        Vector3 xyz2
          = vw::cartography::camera_pixel_to_dem_xyz(camera_ctr, camera_vec,
                                                     dem, dem_georef,
                                                     treat_nodata_as_zero,
                                                     has_intersection,
                                                     height_error_tol, max_abs_tol, max_rel_tol,
                                                     num_max_iter, xyz_guess);
        
        // Quit if we did not find an intersection
        if (!has_intersection)
          continue;
#if 0
        // Average xyz2 which came from the ms pix with the xyz3
        // coming from the pan pix
        camera_ctr = ms_dg_cam->camera_center(pan_pix); 
        camera_vec = ms_dg_cam->pixel_to_vector(pan_pix);
        Vector3 xyz3
          = vw::cartography::camera_pixel_to_dem_xyz(camera_ctr, camera_vec,
                                                     dem, dem_georef,
                                                     treat_nodata_as_zero,
                                                     has_intersection,
                                                     height_error_tol, max_abs_tol, max_rel_tol,
                                                     num_max_iter, xyz_guess);
#endif
        
        if (!has_intersection)
          continue;

        //Vector3 xyz_avg = (xyz2 + xyz3)/2.0;

//         std::cout << "xyz0 diff " << norm_2(xyz2 - xyz) << ' ' << xyz2 - xyz << std::endl;
//         // std::cout << "--xyz diff " << norm_2(xyz - xyz2) << std::endl;

//         Vector3 llh2 = datum.cartesian_to_geodetic(xyz2);
//         Vector3 llh3 = datum.cartesian_to_geodetic(xyz3);

//         std::cout << "llh diff " << llh2-llh3  << std::endl;
        
        Vector2 pan_pix2 = pan_dg_cam->point_to_pixel(xyz2);

        // The disparity from PAN to MS projected in the PAN plane
        dispx(samp_col, samp_row) = pan_pix2.x() - pan_pix.x();
        dispy(samp_col, samp_row) = pan_pix2.y() - pan_pix.y();

        // The matches from PAN to MS
        double sigma = 1.0;
        vw::ip::InterestPoint lip(pan_pix.x(), pan_pix.y(), sigma);
        vw::ip::InterestPoint rip(ms_pix.x(), ms_pix.y(), sigma);
        left_ip.push_back(lip);
        right_ip.push_back(rip);

        std::vector<double> vals;
        vals.push_back(pan_pix.x());
        vals.push_back(pan_pix.y());

        // Note how we use xyz_avg.
        // TODO(oalexan1): Need to get xyz3 as the projection
        // from the PAN pixel, and average that one with xyz2.

        // TODO(oalexan1): See the difference between xyz3 and xyz2!
        
        vals.push_back(xyz2.x());
        vals.push_back(xyz2.y());
        vals.push_back(xyz2.z());
        xyz_list.push_back(vals);
      }
    }

    std::string disp_x_file = opt.out_prefix + "-" + opt.band + "-disp_x.tif";
    std::string disp_y_file = opt.out_prefix + "-" + opt.band + "-disp_y.tif";
    bool has_georef = false, has_nodata = false;
    vw::cartography::GeoReference georef;
    double nodata = 0;
    
    TerminalProgressCallback tpc("asp", ": ");

    //std::cout << "imagesize " << dispx.cols() << ' ' << dispx.rows() << std::endl;
    std::cout << "Writing: " << disp_x_file << std::endl;
    block_write_gdal_image(disp_x_file, dispx, has_georef, georef, has_nodata, nodata, opt, tpc);

    std::cout << "Writing: " << disp_y_file << std::endl;
    block_write_gdal_image(disp_y_file, dispy, has_georef, georef, has_nodata, nodata, opt, tpc);

    std::string match_file = vw::ip::match_filename(opt.out_prefix, opt.pan_image, opt.ms_image);
    std::cout << "Writing: " << match_file << std::endl;
    ip::write_binary_match_file(match_file, left_ip, right_ip);

    // These intersections don't make sense for synthetic jitter
    if (!opt.synthetic_jitter) {
      std::string xyz_file = match_file + ".xyz";
      std::cout << "Writing: " << xyz_file << std::endl;
      std::ofstream ofs (xyz_file.c_str());
      ofs.precision(18);
      for (size_t row = 0; row < xyz_list.size(); row++) {
        for (size_t col = 0; col < xyz_list[row].size(); col++) {
          ofs << xyz_list[row][col] << ' ';
        }
        ofs << "\n";
      }
    }
    
    xercesc::XMLPlatformUtils::Terminate();
  } ASP_STANDARD_CATCHES;
}
