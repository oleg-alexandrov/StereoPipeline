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

/// \file jitter_correct.cc
///

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
    proc_ms_image, out_prefix;
  int num_matches;
};

void handle_arguments(int argc, char *argv[], Options& opt) {
  
  po::options_description general_options("");
  general_options.add_options()
    ("dem",  po::value(&opt.dem)->default_value(""),
     "The DEM to use to undo the processing that was applied to the MS image.")
    ("pan-to-ms-disp",  po::value(&opt.pan_to_ms_disp)->default_value(""),
     "The disparity from the PAN image to the processed MS image.")
    ("num-matches",  po::value(&opt.num_matches)->default_value(0),
     "How many matches to find between the MS and PAN image.")
    ("output-prefix,o",  po::value(&opt.out_prefix),
     "Prefix for output filenames.");

  general_options.add(vw::cartography::GdalWriteOptionsDescription(opt));
  
  po::options_description positional("");
  positional.add_options()
    ("pan-image", po::value(&opt.pan_image))
    ("proc-ms-image", po::value(&opt.proc_ms_image))
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
                             allow_unregistered, unregistered );
  
  if ( !vm.count("pan-image") || !vm.count("proc-ms-image") ||
       !vm.count("pan-camera") || !vm.count("ms-camera") )
    vw::vw_throw( vw::ArgumentErr() << "Not all inputs were specified.\n\n"
                  << usage << general_options );

  if (opt.dem == "") 
    vw::vw_throw( vw::ArgumentErr() << "An input DEM is required.\n\n"
                  << usage << general_options );

  if (opt.pan_to_ms_disp == "") 
    vw::vw_throw( vw::ArgumentErr() << "The processed MS to PAN image disparity is required.\n\n"
                  << usage << general_options );

  if (opt.num_matches <= 0) 
    vw::vw_throw( vw::ArgumentErr() << "Number of matches to compute must be positive.\n\n"
                  << usage << general_options );

  if (opt.out_prefix == "") 
    vw::vw_throw( vw::ArgumentErr() << "The output prefix is required.\n\n"
                  << usage << general_options );

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

  // std::cout << "-xyz is " << xyz << std::endl;
  // Vector2 ms_pix = m_ms_cam->point_to_pixel(xyz);
  // std::cout << "--ms pix " << pix << ' ' << ms_pix << ' ' << norm_2(pix - ms_pix)
  // << std::endl;
          
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
  
  Options opt;
  try {
    xercesc::XMLPlatformUtils::Initialize();

    handle_arguments(argc, argv, opt);

    std::cout << "image and camera: " << opt.pan_image << ' ' << opt.pan_camera << std::endl;
    
    // Initialize the session
    std::string session_str = "dg";
    std::string out_prefix = "out";
    SessionPtr session(asp::StereoSessionFactory::create(session_str, // may change
                                                         opt,
                                                         opt.pan_image,
                                                         opt.proc_ms_image,
                                                         opt.pan_camera,
                                                         opt.ms_camera,
                                                         out_prefix));

    // Load the pan camera
    boost::shared_ptr<vw::camera::CameraModel> pan_cam = session->camera_model(opt.pan_image,
                                                                           opt.pan_camera);

    // Load the ms camera
    boost::shared_ptr<vw::camera::CameraModel> ms_cam = session->camera_model(opt.proc_ms_image,
                                                                           opt.ms_camera);

    asp::DGCameraModel * pan_dg_cam = dynamic_cast<asp::DGCameraModel*>(pan_cam.get());
    asp::DGCameraModel * ms_dg_cam = dynamic_cast<asp::DGCameraModel*>(ms_cam.get());
    

    if (pan_dg_cam == NULL || ms_dg_cam == NULL) 
      vw::vw_throw( vw::ArgumentErr() << "Expecting Digital Globe cameras.\n\n");

    // Read the image
    std::cout << "Reading image: " << opt.pan_image << std::endl;
    bool has_image_nodata = true; // we will surely have no-data on output, even if not on input
    float image_nodata = -std::numeric_limits<float>::max();
    boost::shared_ptr<DiskImageResource> img_rsrc(new DiskImageResourceGDAL(opt.pan_image));
    if (img_rsrc->has_nodata_read())
      image_nodata = img_rsrc->nodata_read();
    ImageViewRef< PixelMask<float> > pan_image
      (create_mask(DiskImageView<float>(opt.pan_image), image_nodata));
    
    // Read the DEM
    std::cout << "dem file is " << opt.dem << std::endl;
    float dem_nodata = -std::numeric_limits<float>::max();
    if (vw::read_nodata_val(opt.dem, dem_nodata))
      vw_out() << "Dem nodata: " << dem_nodata << std::endl;
    ImageViewRef< PixelMask<float> > dem(create_mask(DiskImageView<float>(opt.dem),
                                                     dem_nodata));
    std::cout << "cols and rows are " << dem.cols() << ' ' << dem.rows() << std::endl;
    vw::cartography::GeoReference dem_georef;
    if (!read_georeference(dem_georef, opt.dem))
      vw_throw( ArgumentErr() << "The input DEM " << opt.dem << " has no georeference.\n" );

    std::cout << "--Reading " << opt.pan_to_ms_disp << std::endl;
    DiskImageView< PixelMask<Vector2f> > disparity(opt.pan_to_ms_disp);

    if (disparity.cols() != pan_image.cols() || disparity.rows() != pan_image.rows()) 
      vw_throw( ArgumentErr() << "The input image and disparity must have the same dimensions.\n");
      
    ImageViewRef< PixelMask<Vector2f> > interp_disp  
      = interpolate(disparity, BilinearInterpolation(), ConstantEdgeExtension());
        
    //std::cout << "---num matches is " << opt.num_matches << std::endl;

    int num_cols = pan_image.cols();
    int num_rows = pan_image.rows();

    double ratio = std::sqrt(double(opt.num_matches)/(num_rows*num_cols));

    //std::cout << "--ratio is " << ratio << std::endl;

    int num_sample_cols = std::min(std::max(1, int(round(num_cols * ratio))), num_cols);
    int num_sample_rows = std::min(std::max(1, int(round(num_rows * ratio))), num_rows);

    //std::cout << "num sample cols and rows " << num_sample_cols << ' ' << num_sample_rows
    //          << std::endl;
    

    std::vector< std::vector<double> > xyz_list;

    ImageView<float> dispx;
    ImageView<float> dispy;
    dispx.set_size(num_sample_cols, num_sample_rows);
    dispy.set_size(num_sample_cols, num_sample_rows);
    
    ImageView<float> dispx3;
    ImageView<float> dispy3;
    dispx3.set_size(num_sample_cols, num_sample_rows);
    dispy3.set_size(num_sample_cols, num_sample_rows);

    ImageView<float> dispx4;
    ImageView<float> dispy4;
    dispx4.set_size(num_sample_cols, num_sample_rows);
    dispy4.set_size(num_sample_cols, num_sample_rows);

    std::cout << "num sample cols and rows " << num_sample_cols << ' ' << num_sample_rows
              << std::endl;
    std::cout << "image size " << dispx.cols() << ' ' << dispx.rows() << std::endl;
    std::cout << "--imagesize 4 " << dispy.cols() << ' ' << dispy.rows() << std::endl;
    
    std::vector<double> coeffs;
    char * adj = getenv("ADJ_TMP");
    if (adj != NULL) {
      std::cout << "adj is not null!" << std::endl;
      coeffs.clear();
      std::istringstream ifs(adj);
      double val;
      while (ifs >> val) {
        coeffs.push_back(val);
        std::cout << "--coeff is " << val << std::endl;
      }
    }else{
      std::cout << "--adj is null!" << std::endl;
    }
    
    std::cout << "size of coeffs is " << coeffs.size() << std::endl;

    // Matches from PAN image to raw MS image
    std::vector<vw::ip::InterestPoint> left_ip, right_ip;

    for (int samp_col = 0; samp_col < num_sample_cols; samp_col++) {
      for (int samp_row = 0; samp_row < num_sample_rows; samp_row++) {

        dispx(samp_col, samp_row) = 0;
        dispy(samp_col, samp_row) = 0;
        
        //std::cout << std::endl;

        Vector2 pan_pix = Vector2(samp_col, samp_row)/ratio;

        //std::cout << "--pan pixel val is " << pan_pix << std::endl;

        // Add jitter to cameras
        pan_dg_cam->m_coeffs = coeffs;
        ms_dg_cam->m_coeffs  = coeffs;
        
        // Project from the PAN camera into the ground
        bool   has_intersection     = false;
        bool   treat_nodata_as_zero = false;
        double height_error_tol     = 1e-5;   // error in DEM height
        double max_abs_tol          = 1e-14;  // abs cost function change b/w iterations
        double max_rel_tol          = 1e-14;
        int    num_max_iter         = 100;
        Vector3 xyz_guess           = Vector3();
        Vector3 camera_ctr          = pan_dg_cam->camera_center(pan_pix); 
        Vector3 camera_vec          = pan_dg_cam->pixel_to_vector(pan_pix);
         
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
          continue;
        
        // std::cout << "-xyz is " << xyz << std::endl;
        // Vector2 ms_pix = m_ms_dg_cam->point_to_pixel(xyz);
        // std::cout << "--ms pix " << pix << ' ' << ms_pix << ' ' << norm_2(pix - ms_pix)
        // << std::endl;

        // The MS pix
        Vector2 ms_pix = ms_dg_cam->point_to_pixel(xyz);

        // Now compute the disparity due to the jitter
        // So we use pixels measured with jitter
        // but cameras without jitter.
        pan_dg_cam->m_coeffs = std::vector<double>();
        ms_dg_cam->m_coeffs  = std::vector<double>();
        
        camera_ctr = ms_dg_cam->camera_center(ms_pix); 
        camera_vec = ms_dg_cam->pixel_to_vector(ms_pix);
         
        // Use iterative solver call to compute an intersection of the pixel with the DEM	
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

//         std::cout << "--xyz diff " << xyz - xyz2 << ' ' << norm_2(xyz - xyz2) << std::endl;
        
        Vector2 pan_pix2 = pan_dg_cam->point_to_pixel(xyz2);

//         {
//           bool least_squares = false;
//           double min_triangulation_angle = 1e-8;
//           Vector3 err;
//           vw::stereo::StereoModel sm(pan_dg_cam, ms_dg_cam, least_squares,
//                                      min_triangulation_angle*M_PI/180.0);
//           Vector3 xyz3 = sm(pan_pix, ms_pix, err);

//           // Why is this biased in x???
//           Vector2 pan_pix3 = pan_dg_cam->point_to_pixel(xyz3);
//           Vector2 ms_pix3 = ms_dg_cam->point_to_pixel(xyz3);

//           //std::cout << std::endl;
          
//           //std::cout << "--pan 3 diff " << pan_pix3 << ' ' << pan_pix3 - pan_pix << std::endl;
//           //std::cout << "--ms 3 diff " << ms_pix3 << ' ' << ms_pix3 - ms_pix << std::endl;
          
//            dispx3(samp_col, samp_row) = pan_pix3.x() - pan_pix.x();
//            dispy3(samp_col, samp_row) = pan_pix3.y() - pan_pix.y();

//            Vector2 pan_pix4 = pan_dg_cam->point_to_pixel(xyz);
//            Vector2 ms_pix4 = ms_dg_cam->point_to_pixel(xyz);
           
//            dispx4(samp_col, samp_row) = pan_pix4.x() - pan_pix.x();
//            dispy4(samp_col, samp_row) = pan_pix4.y() - pan_pix.y();

//            //std::cout << "--pan 4 diff " << pan_pix4 << ' ' << pan_pix4 - pan_pix << std::endl;
//            //std::cout << "--ms 4 diff " << ms_pix4 << ' ' << ms_pix4 - ms_pix << std::endl;

//            //std::cout << std::endl;
//         }
        
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

        // Note how we use xyz2!!!!!! Temporary!!!!
        vals.push_back(xyz2.x());
        vals.push_back(xyz2.y());
        vals.push_back(xyz2.z());
        xyz_list.push_back(vals);
        
//         std::cout << "-2pan pix " << pan_pix << ' ' << pan_pix2 << ' ' << pan_pix - pan_pix2
//                   << std::endl;

        //std::cout << std::endl;
      }
    }

    std::string disp_x3_file = "disp_x3.tif";
    std::string disp_y3_file = "disp_y3.tif";

    std::string disp_x4_file = "disp_x4.tif";
    std::string disp_y4_file = "disp_y4.tif";

    std::string disp_x_file = "disp_x.tif";
    std::string disp_y_file = "disp_y.tif";
    bool has_georef = false, has_nodata = false;
    vw::cartography::GeoReference georef;
    double nodata = 0;
    
    TerminalProgressCallback tpc("asp", ": ");

    std::cout << "imagesize " << dispx.cols() << ' ' << dispx.rows() << std::endl;
    
    block_write_gdal_image(disp_x_file, dispx, has_georef, georef, has_nodata, nodata, opt, tpc);
    block_write_gdal_image(disp_y_file, dispy, has_georef, georef, has_nodata, nodata, opt, tpc);

//     std::cout << "Writing " << disp_x3_file << ' ' << disp_y3_file << std::endl;
//     block_write_gdal_image(disp_x3_file, dispx3, has_georef, georef, has_nodata, nodata, opt, tpc);
//     block_write_gdal_image(disp_y3_file, dispy3, has_georef, georef, has_nodata, nodata, opt, tpc);

//     std::cout << "Writing " << disp_x4_file << ' ' << disp_y4_file << std::endl;
//     block_write_gdal_image(disp_x4_file, dispx4, has_georef, georef, has_nodata, nodata, opt, tpc);
//     block_write_gdal_image(disp_y4_file, dispy4, has_georef, georef, has_nodata, nodata, opt, tpc);

    std::string match_file = vw::ip::match_filename(opt.out_prefix, opt.pan_image,
                                                    opt.proc_ms_image);
    std::cout << "Writing: " << match_file << std::endl;
    ip::write_binary_match_file(match_file, left_ip, right_ip);

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
    
    xercesc::XMLPlatformUtils::Terminate();
  } ASP_STANDARD_CATCHES;
}
