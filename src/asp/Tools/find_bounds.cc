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
#include <vw/Camera/CameraModel.h>
#include <vw/Core/Stopwatch.h>
#include <xercesc/util/PlatformUtils.hpp>

#include <ceres/ceres.h>
#include <ceres/loss_function.h>

namespace po = boost::program_options;
namespace fs = boost::filesystem;
using namespace vw;
using namespace vw::camera;

typedef boost::scoped_ptr<asp::StereoSession> SessionPtr;
typedef boost::shared_ptr<vw::camera::CameraModel> CamPtr;
typedef vw::ba::CameraRelationNetwork<vw::ba::JFeature> CRNJ;
typedef vw::ba::CameraNode<vw::ba::JFeature>::iterator crn_iter;

const int NUM_POINT_PARAMS = 3;
const int PIXEL_SIZE = 2;

struct Options : vw::cartography::GdalWriteOptions {
  std::string left_image, right_image, left_camera, right_camera, out_prefix,
    dem, bundle_adjust_prefix;
  int num_frequencies; // TODO(oalexan1): Use instead the frequency, in Hz
  int num_iterations, freeze_num;
  double ms_offset, min_triangulation_angle, horiz_weight, vert_weight, xyz_weight;
  bool read_xyz_on_dem, float_free_term, float_ms_offset;
};

void handle_arguments(int argc, char *argv[], Options& opt) {
  po::options_description general_options("");
  general_options.add_options()
    ("output-prefix,o",  po::value(&opt.out_prefix),
     "Prefix for output filenames.")
    ("dem",  po::value(&opt.dem)->default_value(""),
     "DEM to use.")
    ("bundle-adjust-prefix", po::value(&opt.bundle_adjust_prefix),
     "Adjustment to apply to the cameras.");

  general_options.add(vw::cartography::GdalWriteOptionsDescription(opt));
  
  po::options_description positional("");
  positional.add_options()
    ("left-image", po::value(&opt.left_image))
    ("right-image", po::value(&opt.right_image))
    ("left-camera", po::value(&opt.left_camera))
    ("right-camera", po::value(&opt.right_camera));
  
  po::positional_options_description positional_desc;
  positional_desc.add("left-image", 1);
  positional_desc.add("right-image", 1);
  positional_desc.add("left-camera", 1);
  positional_desc.add("right-camera", 1);
  
  std::string usage("[options] <left-image> <right-image> <left-camera> <right-camera>");
  bool allow_unregistered = false;
  std::vector<std::string> unregistered;
  po::variables_map vm =
    asp::check_command_line(argc, argv, opt, general_options, general_options,
                            positional, positional_desc, usage,
                             allow_unregistered, unregistered );
  
  if ( !vm.count("left-image")  || !vm.count("right-image")  ||
       !vm.count("left-camera") || !vm.count("right-camera") )
    vw::vw_throw( vw::ArgumentErr() << "Not all inputs were specified.\n\n"
                  << usage << general_options );
  
    // Need this to be able to load adjusted camera models. That will happen
  // in the stereo session.
  asp::stereo_settings().bundle_adjust_prefix = opt.bundle_adjust_prefix;

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

    std::cout << "left image and camera: " << opt.left_image << ' ' << opt.left_camera << "\n";
    std::cout << "right image and camera: " << opt.right_image << ' ' << opt.right_camera << "\n";

    std::cout << "dem is " << opt.dem << std::endl;
    float dem_nodata_val = -std::numeric_limits<float>::max();
    if (vw::read_nodata_val(opt.dem, dem_nodata_val)){
      vw_out() << "Dem nodata: " << dem_nodata_val << std::endl;
    }

    // Load the cameras
    CamPtr left_cam = load_camera(opt, opt.left_image, opt.left_camera);
    CamPtr right_cam = load_camera(opt, opt.right_image, opt.right_camera);

    AdjustedCameraModel * left_adj_cam = dynamic_cast<AdjustedCameraModel*>(left_cam.get());
    AdjustedCameraModel * right_adj_cam = dynamic_cast<AdjustedCameraModel*>(right_cam.get());
    if (left_adj_cam == NULL || right_adj_cam == NULL) 
      vw_throw( ArgumentErr() << "get_isis_cam: Expecting adjusted camera models.\n" );

    ImageView< PixelMask<float> > dem (create_mask(DiskImageView<float>(opt.dem), dem_nodata_val));
    std::cout << "cols and rows are " << dem.cols() << ' ' << dem.rows() << std::endl;
    vw::cartography::GeoReference georef;
    if (!read_georeference(georef, opt.dem))
      vw_throw( ArgumentErr() << "The input DEM " << opt.dem << " has no georeference.\n" );

    BBox2i left_box, right_box;
    for (int col = 0; col < dem.cols(); col += 10) {
      for (int row = 0; row < dem.rows(); row += 10) {
        PixelMask<float> val = dem(col, row);
        if (!is_valid(val))
          continue;
        Vector3 llh; 
        subvector(llh, 0, 2) = georef.pixel_to_lonlat(Vector2(col, row));
        llh[2] = val.child();
        Vector3 xyz = georef.datum().geodetic_to_cartesian(llh);

        Vector2 left_pix = left_cam->point_to_pixel(xyz);
        left_box.grow(left_pix);
        
        Vector2 right_pix = right_cam->point_to_pixel(xyz);
        right_box.grow(right_pix);
      }
    }

    std::cout << "left box " << left_box.min().x() << ' ' << left_box.min().y() << ' '
              << left_box.width() << ' ' << left_box.height() << std::endl;
    std::cout << "right box " << right_box.min().x() << ' ' << right_box.min().y() << ' '
              << right_box.width() << ' ' << right_box.height() << std::endl;
    
    xercesc::XMLPlatformUtils::Terminate();
  } ASP_STANDARD_CATCHES;

  return 0;
}
