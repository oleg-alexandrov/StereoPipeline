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

#ifndef __ASP_CORE_STEREO_SETTINGS_H__
#define __ASP_CORE_STEREO_SETTINGS_H__

#include <asp/Core/Common.h>
#include <boost/program_options.hpp>
#include <boost/program_options/detail/config_file.hpp>

namespace asp {

  // Program Options for each executable/step
  struct PreProcessingDescription : public boost::program_options::options_description {
    PreProcessingDescription();
  };
  struct CorrelationDescription : public boost::program_options::options_description {
    CorrelationDescription();
  };
  struct SubpixelDescription : public boost::program_options::options_description {
    SubpixelDescription();
  };
  struct FilteringDescription : public boost::program_options::options_description {
    FilteringDescription();
  };
  struct TriangulationDescription : public boost::program_options::options_description {
    TriangulationDescription();
  };
  struct DGDescription : public boost::program_options::options_description {
    DGDescription();
  };

  boost::program_options::options_description
  generate_config_file_options( asp::BaseOptions& opt );

  // Structure holding variables
  class StereoSettings {
  public:
    void validate();
    void write_copy( int argc, char *argv[],
                     std::string const& input_file,
                     std::string const& output_file ) const;
    bool is_search_defined() const;

    // ----------------
    // Public variables
    // ----------------

    // Preprocessing options
    std::string alignment_method;     // Valid options are: [Homography, Epipolar, None]
    bool individually_normalize;      // if > 1, normalize the images
                                      //         individually with their
                                      //         own hi's and lo's
    bool force_max_min;               // Use entire dynamic range of image.
    double nodata_threshold;          // Pixels with value less than this are treated as no-data
    double nodata_percentage;         // the percentage of low-value pixels treated as no-data
    double nodata_optimal_threshold_factor; // Pixels with values less than this factor times the optimal Otsu threshold
                                      // are treated as no-data
    
    // Correlation Options
    float slogW;                      // Preprocessing filter width
    vw::uint16 pre_filter_mode;       // 0 = None
                                      // 1 = Gaussian Blur
                                      // 2 = Log Filter
                                      // 3 = SLog Filter  
    vw::uint16  seed_mode;            // 0 = User global search for each tile
                                      // 1 = Narrow search for each tile to low
                                      //     res disparity seed
                                      // 2 = Affine Transform and narrow search
                                      //     based on disparity seed 
    float seed_percent_pad;           // Pad amound towards the IP found
    float xcorr_threshold;
    vw::uint16 cost_mode;             // 0 = absolute difference
                                      // 1 = squared difference
                                      // 2 = normalized cross correlation 
    vw::uint16 corr_max_levels;       // Max pyramid levels to process. 0 hits only once.
    double edge_density_threshold;    // Minimum acceptable edge density at a given pyramid level

                                      // search range, used to define the
                                      // search range of D_sub.
    vw::Vector2i kernel;              // Correlation kernel
    vw::Vector2i subpixel_kernel;     // Subpixel correlation kernel
    vw::BBox2i search_range;          // Correlation search range
    bool disable_h_subpixel, disable_v_subpixel;
    vw::uint16 subpixel_mode;         // 0 = parabola fitting
                                      // 1 = affine, robust weighting
                                      // 2 = affine, bayes weighting
                                      // 3 = affine, bayes EM weighting 
    vw::uint16 subpixel_max_levels;   // Max pyramid levels to process. 0 hits only once.

    // EMSubpixelCorrelator Options (mode 3 only)
    int subpixel_affine_iter;
    int subpixel_em_iter;
    int subpixel_pyramid_levels;

    // Filtering Options
    vw::Vector2i rm_half_kernel;      // Low confidence pixel removal kernel size 
    int rm_min_matches;               // Min # of pxl to be matched to keep pxl 
    int rm_threshold;                 // rm_treshold < disp[n]-disp[m] reject pxl 
    int rm_cleanup_passes;            // Number of times to perform cleanup
                                      // in the post-processing phase 
    int erode_max_size;               // Max island size in pixels that it'll remove
    bool disable_fill_holes;
    int fill_hole_max_size;           // Maximum hole size in pixels that we'll attempt
                                      // to fill 
    bool mask_flatfield;              // Masks pixels in the input images that are less
                                      // than 0 (for use with Apollo Metric Camera)

    // Triangulation Options
    std::string universe_center;      // Center for the radius clipping   
    float near_universe_radius;       // Radius of the universe in meters 
    float far_universe_radius;        // Radius of the universe in meters 
    bool use_least_squares;           // Use a more rigorous triangulation
    bool compute_error_vector;        // Compute the triangulation error vector, not just its length

    // DG Options
    bool disable_correct_velocity_aberration;
  };

  /// Return the singleton instance of the stereo setting structure.
  /// The stereo settings struct is created the first time this method
  /// is invoked.  You should *always* access the stereo settings
  /// through this function.
  StereoSettings& stereo_settings();

  // Custom readers for Boost Program Options
  class asp_config_file_iterator : public boost::program_options::detail::common_config_file_iterator {
    boost::shared_ptr<std::basic_istream<char> > is;
  private:
    bool getline(std::string& s); // Used to precondition string before reading
  public:
    asp_config_file_iterator() {
      found_eof();
    }

    // Creates a config file parser for the specified stream.
    asp_config_file_iterator(std::basic_istream<char>& is,
                             const std::set<std::string>& allowed_options,
                             bool allow_unregistered = false);
  };

  // Custom Parsers for ASP's stereo.default files
  boost::program_options::basic_parsed_options<char>
  parse_asp_config_file( std::basic_istream<char>&,
                         const boost::program_options::options_description&,
                         bool allow_unregistered = false );

  boost::program_options::basic_parsed_options<char>
  parse_asp_config_file( std::string const&,
                         const boost::program_options::options_description&,
                         bool allow_unregistered = false );

}

#endif//__ASP_CORE_STEREO_SETTINGS_H__
