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
#include <xercesc/util/PlatformUtils.hpp>

namespace po = boost::program_options;
namespace fs = boost::filesystem;
using namespace vw;

struct Options : vw::cartography::GdalWriteOptions {
  std::string proc_ms_image, pan_camera, ms_camera, bundle_adjust_prefix, dem, ms_to_pan_disp,
    out_raw_ms_image, out_match_file;
  int num_matches;
};

void handle_arguments(int argc, char *argv[], Options& opt) {
  
  po::options_description general_options("");
  general_options.add_options()
    ("dem",  po::value(&opt.dem)->default_value(""),
     "The DEM to use to undo the processing that was applied to the MS image.")
    ("ms-to-pan-disp",  po::value(&opt.ms_to_pan_disp)->default_value(""),
     "The disparity from the processed MS image to the PAN image.")
    ("num-matches",  po::value(&opt.num_matches)->default_value(0),
     "How many matches to find between the MS and PAN image.")
    ("out-match-file",  po::value(&opt.out_match_file)->default_value(""),
     "Save here a set of matches from the raw MS image to the PAN image.")
    ("bundle-adjust-prefix", po::value(&opt.bundle_adjust_prefix),
     "Adjustment to apply to the cameras.");
    
  general_options.add(vw::cartography::GdalWriteOptionsDescription(opt));
  
  po::options_description positional("");
  positional.add_options()
    ("proc-ms-image", po::value(&opt.proc_ms_image))
    ("pan-camera", po::value(&opt.pan_camera))
    ("ms-camera", po::value(&opt.ms_camera))
    ("raw-ms-image", po::value(&opt.out_raw_ms_image));
  
  po::positional_options_description positional_desc;
  positional_desc.add("proc-ms-image", 1);
  positional_desc.add("pan-camera", 1);
  positional_desc.add("ms-camera", 1);
  positional_desc.add("raw-ms-image", 1);
  
  std::string usage("[options] <proc-ms-image> <pan-camera> <ms-camera> <raw-ms-image>");
  bool allow_unregistered = false;
  std::vector<std::string> unregistered;
  po::variables_map vm =
    asp::check_command_line(argc, argv, opt, general_options, general_options,
                            positional, positional_desc, usage,
                             allow_unregistered, unregistered );
  
  if ( !vm.count("proc-ms-image") || !vm.count("pan-camera") || !vm.count("ms-camera") ||
       !vm.count("raw-ms-image") )
    vw::vw_throw( vw::ArgumentErr() << "Not all inputs were specified.\n\n"
                  << usage << general_options );

  if (opt.dem == "") 
    vw::vw_throw( vw::ArgumentErr() << "An input DEM is required.\n\n"
                  << usage << general_options );

  if (opt.ms_to_pan_disp == "") 
    vw::vw_throw( vw::ArgumentErr() << "The processed MS to PAN image disparity is required.\n\n"
                  << usage << general_options );

  if (opt.out_match_file == "") 
    vw::vw_throw( vw::ArgumentErr() << "The output match file is required.\n\n"
                  << usage << general_options );

  if (opt.num_matches <= 0) 
    vw::vw_throw( vw::ArgumentErr() << "Number of matches to compute must be positive.\n\n"
                  << usage << general_options );
  
  // Need this to be able to load adjusted camera models. That will happen
  // in the stereo session.
  asp::stereo_settings().bundle_adjust_prefix = opt.bundle_adjust_prefix;

  // Create the directory in which the output image will be written.
  vw::create_out_dir(opt.out_raw_ms_image);
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


// Class to make an image seen from the PAN sensor perspective into an
// image seen from an MS sensor perspective. For that, one has to
// travel from an MS camera pixel to the DEM representing the ground,
// and then project into the PAN camera.
class PanToMsCamView: public ImageViewBase<PanToMsCamView> {

  ImageViewRef< PixelMask<float> > m_pan_image, m_dem;
  vw::cartography::GeoReference m_dem_georef;
  boost::shared_ptr<vw::camera::CameraModel> m_pan_cam, m_ms_cam;
  
  typedef PixelMask<float> PixelT;

public:
  PanToMsCamView(ImageViewRef< PixelMask<float> > pan_image, 
                 ImageViewRef< PixelMask<float> > dem, 
                 vw::cartography::GeoReference const& dem_georef,
                 boost::shared_ptr<vw::camera::CameraModel> pan_cam,
                 boost::shared_ptr<vw::camera::CameraModel> ms_cam):
    m_pan_image(pan_image), m_dem(dem), m_dem_georef(dem_georef),
    m_pan_cam(pan_cam), m_ms_cam(ms_cam) {}
  
  typedef PixelT pixel_type;
  typedef PixelT result_type;
  typedef ProceduralPixelAccessor<PanToMsCamView> pixel_accessor;

  inline int32 cols() const { return m_pan_image.cols(); }
  inline int32 rows() const { return m_pan_image.rows(); }
  inline int32 planes() const { return 1; }

  inline pixel_accessor origin() const { return pixel_accessor( *this, 0, 0 ); }

  inline pixel_type operator()( double/*i*/, double/*j*/, int32/*p*/ = 0 ) const {
    vw_throw(NoImplErr() << "PanToMsCamView::operator()(...) is not implemented");
    return pixel_type();
  }

  typedef CropView<ImageView<pixel_type> > prerasterize_type;
  inline prerasterize_type prerasterize(BBox2i const& bbox) const {

    ImageView<result_type> tile(bbox.width(), bbox.height());
    
    ImageViewRef< PixelMask<float> > interp_pan_img 
      = interpolate(m_pan_image, BicubicInterpolation(), ConstantEdgeExtension());
    
     for (int col = bbox.min().x(); col < bbox.max().x(); col++) {
       for (int row = bbox.min().y(); row < bbox.max().y(); row++) {

         // Initialize the output
         tile(col - bbox.min().x(), row - bbox.min().y()) = 0;
         tile(col - bbox.min().x(), row - bbox.min().y()).invalidate();

         // Project from the MS camera into the ground
         Vector2 ms_pix(col, row);
         Vector2 pan_pix;
         bool ans = raw_to_proc_ms_pix(ms_pix, m_dem, m_dem_georef, m_pan_cam, m_ms_cam, pan_pix);
         if (!ans) 
           continue;
         
         if (pan_pix.x() < 1 || pan_pix.x() > m_pan_image.cols() - 2) continue;
         if (pan_pix.y() < 1 || pan_pix.y() > m_pan_image.rows() - 2) continue;

         // Do bicubic interpolation 
         tile(col - bbox.min().x(), row - bbox.min().y())
           = interp_pan_img(pan_pix.x(), pan_pix.y());
       }
     }
     
    return prerasterize_type(tile, -bbox.min().x(), -bbox.min().y(), cols(), rows());
  }

  template <class DestT>
  inline void rasterize(DestT const& dest, BBox2i bbox) const {
    vw::rasterize(prerasterize(bbox), dest, bbox);
  }
};

PanToMsCamView raw_to_proc_ms_image(ImageViewRef< PixelMask<float> > pan_image, 
                             ImageViewRef< PixelMask<float> > dem, 
                             vw::cartography::GeoReference const& dem_georef,
                             boost::shared_ptr<vw::camera::CameraModel> pan_cam,
                             boost::shared_ptr<vw::camera::CameraModel> ms_cam) {
  return PanToMsCamView(pan_image, dem, dem_georef, pan_cam, ms_cam);
}

typedef boost::scoped_ptr<asp::StereoSession> SessionPtr;

int main(int argc, char* argv[]) {
  
  Options opt;
  try {
    xercesc::XMLPlatformUtils::Initialize();

    handle_arguments(argc, argv, opt);

    std::cout << "image and camera: " << opt.proc_ms_image << ' ' << opt.pan_camera << std::endl;
    
    // Initialize the session
    std::string session_str = "dg";
    std::string out_prefix = "out";
    SessionPtr session(asp::StereoSessionFactory::create(session_str, // may change
                                                         opt,
                                                         opt.proc_ms_image,
                                                         opt.out_raw_ms_image,
                                                         opt.pan_camera,
                                                         opt.ms_camera,
                                                         out_prefix));

    // Load the pan camera
    boost::shared_ptr<vw::camera::CameraModel> pan_cam = session->camera_model(opt.proc_ms_image,
                                                                           opt.pan_camera);

    // Load the ms camera
    boost::shared_ptr<vw::camera::CameraModel> ms_cam = session->camera_model(opt.out_raw_ms_image,
                                                                           opt.ms_camera);

    // Read the image
    std::cout << "Reading image: " << opt.proc_ms_image << std::endl;
    bool has_image_nodata = true; // we will surely have no-data on output, even if not on input
    float image_nodata = -std::numeric_limits<float>::max();
    boost::shared_ptr<DiskImageResource> img_rsrc(new DiskImageResourceGDAL(opt.proc_ms_image));
    if (img_rsrc->has_nodata_read())
      image_nodata = img_rsrc->nodata_read();
    ImageViewRef< PixelMask<float> > pan_image
      (create_mask(DiskImageView<float>(opt.proc_ms_image), image_nodata));
    
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

#if 0
    // To actually save a transformed image rather than use the transform at selected
    // pixels
    ImageViewRef< PixelMask<float> > ms_image = raw_to_proc_ms_image(pan_image, dem, dem_georef,
                                                                   pan_cam,  ms_cam);
    
    vw_out() << "Writing: " << opt.out_raw_ms_image << std::endl;
    bool has_image_georef = false;
    vw::cartography::GeoReference image_georef;
    
    TerminalProgressCallback tpc("asp", "\t--> ");
    vw::cartography::block_write_gdal_image(opt.out_raw_ms_image,
                                            apply_mask(ms_image, image_nodata),
                                            has_image_georef, image_georef,
                                            has_image_nodata, image_nodata,
                                            opt, tpc);
#endif
    
    std::cout << "--Reading " << opt.ms_to_pan_disp << std::endl;
    DiskImageView< PixelMask<Vector2f> > disparity(opt.ms_to_pan_disp);

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
    
    // Matches from raw MS image to the processed PAN image
    std::vector<vw::ip::InterestPoint> left_ip, right_ip;

    for (int samp_col = 0; samp_col < num_sample_cols; samp_col++) {
      for (int samp_row = 0; samp_row < num_sample_rows; samp_row++) {
        
        //std::cout << std::endl;

        Vector2 raw_ms_pix = Vector2(samp_col, samp_row)/ratio;
        //std::cout << "--raw ms pix is " << raw_ms_pix << std::endl;
        
        Vector2 proc_ms_pix;
        bool ans = raw_to_proc_ms_pix(raw_ms_pix, dem, dem_georef, pan_cam, ms_cam, proc_ms_pix);
        if (!ans) 
          continue;

        //std::cout << "--proc_ms_pix " << proc_ms_pix << std::endl;

        if (proc_ms_pix.x() < 0 || proc_ms_pix.x() > interp_disp.cols() - 1)
          continue;
        if (proc_ms_pix.y() < 0 || proc_ms_pix.y() > interp_disp.rows() - 1)
          continue;

        PixelMask<Vector2f> disp_val = interp_disp(proc_ms_pix.x(), proc_ms_pix.y());
        if (!is_valid(disp_val)) 
          continue;

        //std::cout << "--disp val is " << disp_val << std::endl;
        
        Vector2 pan_pix = proc_ms_pix + disp_val.child();

        //std::cout << "--match " << raw_ms_pix << ' ' << pan_pix << std::endl;
        
        double sigma = 1.0;
        vw::ip::InterestPoint lip(raw_ms_pix.x(), raw_ms_pix.y(), sigma);
        vw::ip::InterestPoint rip(pan_pix.x(), pan_pix.y(), sigma);
        left_ip.push_back(lip);
        right_ip.push_back(rip);
      }
    }

    std::cout << "Writing: " << opt.out_match_file << std::endl;
    ip::write_binary_match_file(opt.out_match_file, left_ip, right_ip);
    
    xercesc::XMLPlatformUtils::Terminate();
  } ASP_STANDARD_CATCHES;
}
