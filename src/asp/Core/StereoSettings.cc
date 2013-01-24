// __BEGIN_LICENSE__
//  Copyright (c) 2009-2012, United States Government as represented by the
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


/// \file StereoSettings.h
///
#include <fstream>

#include <vw/Core/Thread.h>
#include <vw/Core/Log.h>

#include <asp/Core/Common.h>
#include <asp/Core/StereoSettings.h>

namespace po = boost::program_options;
using namespace vw;

// Create a single instance of the StereoSettings
// ---------------------------------------------------
namespace {
  vw::RunOnce stereo_settings_once = VW_RUNONCE_INIT;
  boost::shared_ptr<asp::StereoSettings> stereo_settings_ptr;
  void init_stereo_settings() {
    stereo_settings_ptr = boost::shared_ptr<asp::StereoSettings>(new asp::StereoSettings());
  }
}

namespace asp {
  StereoSettings& stereo_settings() {
    stereo_settings_once.run( init_stereo_settings );
    return *stereo_settings_ptr;
  }

  // Define our options that are available
  // ----------------------------------------------------------
  PreProcessingDescription::PreProcessingDescription() : po::options_description("Preprocessing Options") {

    StereoSettings& global = stereo_settings();

    double nan = std::numeric_limits<double>::quiet_NaN();
    global.nodata_threshold = nan;
    global.nodata_percentage = nan;
    global.nodata_optimal_threshold_factor = nan;

    (*this).add_options()
      ("alignment-method", po::value(&global.alignment_method),
       "Rough alignment for input images. [Homography, Epipolar, None]")
      ("individually-normalize", po::bool_switch(&global.individually_normalize),
       "Individually normalize the input images between 0.0-1.0 using +- 2.5 sigmas about their mean values.")
      ("force-use-entire-range", po::bool_switch(&global.force_max_min),
       "Normalize images based on the global min and max values from both images. Don't use this option if you are using normalized cross correlation.")
      ("nodata-threshold", po::value(&global.nodata_threshold),
       "Pixels with values less than this number times the maximum image value are treated as no-data.")
      ("nodata-percentage", po::value(&global.nodata_percentage),
       "The percentage of (low-value) pixels treated as no-data.")
      ("nodata-optimal-threshold-factor", po::value(&global.nodata_optimal_threshold_factor),
       "Pixels with values less than this factor times the optimal Otsu threshold are treated as no-data. Suggested value: 0.1 to 0.2.");
  }

  CorrelationDescription::CorrelationDescription() : po::options_description("Correlation Options") {
    StereoSettings& global = stereo_settings();
    (*this).add_options()
      ("prefilter-kernel-width", po::value(&global.slogW)->default_value(1.5),
       "Sigma value for gaussian kernel used in prefilter for correlator.")
      ("prefilter-mode", po::value(&global.pre_filter_mode)->default_value(2),
       "Preprocessing filter mode. [0 None, 1 Gaussian, 2 LoG, 3 Sign of LoG]")
      ("corr-seed-mode", po::value(&global.seed_mode)->default_value(1),
       "Correlation seed strategy. [0 None, 1 Use Low-Res Disparity]")
      ("corr-sub-seed-percent", po::value(&global.seed_percent_pad)->default_value(.25),
       "Percent fudge factor for disparity seed's search range")
      ("cost-mode", po::value(&global.cost_mode)->default_value(2),
       "Correlation cost metric. [0 Absolute, 1 Squared, 2 Normalized Cross Correlation]")
      ("xcorr-threshold", po::value(&global.xcorr_threshold)->default_value(2),
       "L-R vs R-L agreement threshold in pixels.")
      ("corr-kernel", po::value(&global.kernel)->default_value(Vector2i(21,21),"21 21"),
       "Kernel size used for integer correlator.")
      ("corr-search", po::value(&global.search_range)->default_value(BBox2i(0,0,0,0), "auto"),
       "Disparity search range. Specify in format: hmin vmin hmax vmax.")
      ("corr-max-levels", po::value(&global.corr_max_levels)->default_value(5),
       "Max pyramid levels to process when using the integer correlator. (0 is just a single level).")
      ("edge-density-threshold", po::value(&global.edge_density_threshold)->default_value(0),
       "Minimum acceptable edge density at a given pyramid level.");
    
    po::options_description backwards_compat_options("Aliased backwards compatibility options");
    backwards_compat_options.add_options()
      ("h-kernel", po::value(&global.kernel[0]), "Correlation kernel width")
      ("v-kernel", po::value(&global.kernel[1]), "Correlation kernel height")
      ("h-corr-min", po::value(&global.search_range.min()[0]), "Correlation window size min x")
      ("h-corr-max", po::value(&global.search_range.max()[0]), "Correlation window size max x")
      ("v-corr-min", po::value(&global.search_range.min()[1]), "Correlation window size min y")
      ("v-corr-max", po::value(&global.search_range.max()[1]), "Correlation window size max y");
    (*this).add( backwards_compat_options );
  }

  SubpixelDescription::SubpixelDescription() : po::options_description("Subpixel Options") {
    StereoSettings& global = stereo_settings();
    (*this).add_options()
      ("subpixel-mode",  po::value(&global.subpixel_mode)->default_value(2),
       "Subpixel algorithm. [0 None, 1 Parabola, 2 Bayes EM]")
      ("subpixel-kernel", po::value(&global.subpixel_kernel)->default_value(Vector2i(35,35), "35 35"),
       "Kernel size used for subpixel method.")
      ("disable-h-subpixel", po::bool_switch(&global.disable_h_subpixel),
       "Disable calculation of subpixel in horizontal direction.")
      ("disable-v-subpixel", po::bool_switch(&global.disable_v_subpixel),
       "Disable calculation of subpixel in vertical direction.")
      ("subpixel-max-levels", po::value(&global.subpixel_max_levels)->default_value(2),
       "Max pyramid levels to process when using the BayesEM refinement. (0 is just a single level).");

    po::options_description experimental_subpixel_options("Experimental Subpixel Options");
    experimental_subpixel_options.add_options()
      ("subpixel-em-iter", po::value(&global.subpixel_em_iter)->default_value(15),
       "Maximum number of EM iterations for EMSubpixelCorrelator")
      ("subpixel-affine-iter", po::value(&global.subpixel_affine_iter)->default_value(5),
       "Maximum number of affine optimization iterations for EMSubpixelCorrelator")
      ("subpixel-pyramid-levels", po::value(&global.subpixel_pyramid_levels)->default_value(3),
       "Number of pyramid levels for EMSubpixelCorrelator");
    (*this).add( experimental_subpixel_options );

    po::options_description backwards_compat_options("Aliased backwards compatibility options");
    backwards_compat_options.add_options()
      ("subpixel-h-kernel", po::value(&global.subpixel_kernel[0]), "Subpixel kernel width")
      ("subpixel-v-kernel", po::value(&global.subpixel_kernel[1]), "Subpixel kernel height");
    (*this).add( backwards_compat_options );
  }

  FilteringDescription::FilteringDescription() : po::options_description("Filtering Options") {
    StereoSettings& global = stereo_settings();
    (*this).add_options()
      ("rm-half-kernel", po::value(&global.rm_half_kernel)->default_value(Vector2i(5,5), "5 5"),
       "Low confidence pixel removal kernel (half sized)")
      ("rm-min-matches", po::value(&global.rm_min_matches)->default_value(60),
       "Minimum number of pixels to be matched to keep sample")
      ("rm-threshold", po::value(&global.rm_threshold)->default_value(3),
       "Maximum distance between samples to be considered still matched")
      ("rm-cleanup-passes", po::value(&global.rm_cleanup_passes)->default_value(1),
       "Number of passes for cleanup during the post-processing phase")
      ("erode-max-size", po::value(&global.erode_max_size)->default_value(1000),
       "Max size of islands that should be removed")
      ("disable-fill-holes", po::bool_switch(&global.disable_fill_holes),
       "Disable filling of holes using an inpainting method")
      ("fill-holes-max-size", po::value(&global.fill_hole_max_size)->default_value(100000),
       "Max size in pixels of holes that can be filled in.")
      ("mask-flatfield", po::bool_switch(&global.mask_flatfield),
       "Mask dust found on the sensor or film. (For use with Apollo Metric Cameras only!)");

    po::options_description backwards_compat_options("Aliased backwards compatibility options");
    backwards_compat_options.add_options()
      ("rm-h-half-kern", po::value(&global.rm_half_kernel[0]), "Filter kernel half width")
      ("rm-v-half-kern", po::value(&global.rm_half_kernel[1]), "Filter kernel half height");
    (*this).add( backwards_compat_options );
  }

  TriangulationDescription::TriangulationDescription() : po::options_description("Triangulation Options") {
    StereoSettings& global = stereo_settings();
    (*this).add_options()
      ("universe-center", po::value(&global.universe_center)->default_value("None"),
       "Center for radius measurement thresholding. [Camera, Zero, None]")
      ("near-universe-radius", po::value(&global.near_universe_radius)->default_value(0.0),
       "Radius of inner boundary of universe in meters.")
      ("far-universe-radius", po::value(&global.far_universe_radius)->default_value(0.0),
       "Radius of outer boundary of universe in meters.")
      ("use-least-squares", po::bool_switch(&global.use_least_squares)->default_value(false)->implicit_value(true),
       "Use rigorous least squares triangulation process. This is slow for ISIS processes.")
      ("compute-error-vector", po::bool_switch(&global.compute_error_vector)->default_value(false)->implicit_value(true),
       "Compute the triangulation error vector, not just its length.")
      ;
  }

  DGDescription::DGDescription() : po::options_description("DG Options") {
    StereoSettings& global = stereo_settings();
    (*this).add_options()
      ("disable-correct-velocity-aberration", po::bool_switch(&global.disable_correct_velocity_aberration)->default_value(false)->implicit_value(true),
       "Apply the velocity aberration correction for Digital Globe cameras.");
  }

  po::options_description
  generate_config_file_options( asp::BaseOptions& opt ) {
    po::options_description cfg_options;
    cfg_options.add( asp::BaseOptionsDescription( opt ) );
    cfg_options.add( PreProcessingDescription() );
    cfg_options.add( CorrelationDescription() );
    cfg_options.add( SubpixelDescription() );
    cfg_options.add( FilteringDescription() );
    cfg_options.add( TriangulationDescription() );
    cfg_options.add( DGDescription() );

    return cfg_options;
  }

  void StereoSettings::validate() {
    using namespace boost::algorithm;
    to_lower( alignment_method );
    trim( alignment_method );
    VW_ASSERT( alignment_method == "none" || alignment_method == "homography" ||
               alignment_method == "epipolar",
               ArgumentErr() << "\"" <<  alignment_method
               << "\" is not a valid option for ALIGNMENT_METHOD." );

    to_lower( universe_center );
    trim( universe_center );
    VW_ASSERT( universe_center == "camera" || universe_center == "zero" ||
               universe_center == "none",
               ArgumentErr() << "\"" << universe_center
               << "\" is not a valid option for UNIVERSE_CENTER." );
  }

  void StereoSettings::write_copy( int argc, char *argv[],
                                   std::string const& input_file,
                                   std::string const& output_file ) const {
    using namespace std;
    ifstream in( input_file.c_str() );
    ofstream out( output_file.c_str() );

    // Write some log information
    out << "# ASP Stereo Configuration Copy:" << endl;
    out << "# " << current_posix_time_string() << endl;
    out << "# > ";
    for ( int i = 0; i < argc; i++ ) {
      if ( i )
        out << " ";
      out << string(argv[i]);
    }
    out << endl << endl;

    out << in.rdbuf();
    in.close();
    out.close();
  }

  bool StereoSettings::is_search_defined() const {
    return !( search_range.min() == Vector2i() &&
              search_range.max() == Vector2i() );
  }

  bool asp_config_file_iterator::getline( std::string& s ) {
    std::string ws;

    if ( !std::getline( *is, ws, '\n') )
      return false;

    // Remove any comments that might be one the line
    size_t n = ws.find('#');
    if ( n != std::string::npos ) {
      ws = ws.substr(0,n);
      boost::trim(ws);
    }

    // Handle empty lines .. just pass them on through
    if ( ws.empty() ) {
      s = ws;
      return true;
    }

    // If there is not an equal sign, the first space is turned to
    // equal ... or it's just appended. Also, use lower case for the key.
    n = ws.find('=');
    if ( n  == std::string::npos ) {
      n = ws.find(' ');

      if ( n == std::string::npos ) {
        ws += "=";
        boost::to_lower(ws);
      } else {
        ws[n] = '=';
        std::string lowered_key = boost::to_lower_copy( ws.substr(0,n) );
        ws.replace( 0, n, lowered_key );
      }
    } else {
      std::string lowered_key = boost::to_lower_copy( ws.substr(0,n) );
      ws.replace( 0, n, lowered_key );
    }
    s = ws;
    return true;
  }

  asp_config_file_iterator::asp_config_file_iterator( std::basic_istream<char>& is,
                                                      const std::set<std::string>& allowed_options,
                                                      bool allow_unregistered ) :
    po::detail::common_config_file_iterator(allowed_options, allow_unregistered) {
    this->is.reset(&is, po::detail::null_deleter());
    get();
  }

  po::basic_parsed_options<char>
  parse_asp_config_file( std::basic_istream<char>& is,
                         const po::options_description& desc,
                         bool allow_unregistered ) {
    std::set<std::string> allowed_options;

    const std::vector<boost::shared_ptr<po::option_description> >& options = desc.options();
    for (size_t i = 0; i < options.size(); ++i) {
      const po::option_description& d = *options[i];

        if (d.long_name().empty())
          boost::throw_exception(
                                 po::error("long name required for config file"));

        allowed_options.insert(d.long_name());
      }

    // Parser return char strings
    po::parsed_options result(&desc);
    std::copy(asp_config_file_iterator(is, allowed_options, allow_unregistered),
              asp_config_file_iterator(),
              std::back_inserter(result.options));
    // Convert char strings into desired type.
    return po::basic_parsed_options<char>(result);
  }

  po::basic_parsed_options<char>
  parse_asp_config_file( std::string const& filename,
                         const po::options_description& desc,
                         bool allow_unregistered ) {
    std::basic_ifstream<char> strm(filename.c_str());
    if (!strm) {
      boost::throw_exception(po::reading_file(filename.c_str()));
    }
    return parse_asp_config_file(strm, desc, allow_unregistered);
  }
} // end namespace asp
