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


/// \file point2las.cc
///

#include <fstream>
#include <iostream>
#include <string>
#include <boost/program_options.hpp>
#include <liblas/liblas.hpp>

#include <asp/Core/Macros.h>
#include <asp/Core/Common.h>
#include <asp/Core/PointUtils.h>

#include <vw/FileIO.h>
#include <vw/Image.h>
#include <vw/Math.h>
#include <vw/Cartography/PointImageManipulation.h>

using namespace vw;
namespace po = boost::program_options;

// Allows FileIO to correctly read/write these pixel types
namespace vw {
  typedef Vector<float64,6> Vector6;
  template<> struct PixelFormatID<Vector3>   { static const PixelFormatEnum value = VW_PIXEL_GENERIC_3_CHANNEL; };
  template<> struct PixelFormatID<Vector4>   { static const PixelFormatEnum value = VW_PIXEL_GENERIC_4_CHANNEL; };
  template<> struct PixelFormatID<Vector6>   { static const PixelFormatEnum value = VW_PIXEL_GENERIC_6_CHANNEL; };
}

struct Options : asp::BaseOptions {
  // Input
  std::string reference_spheroid;
  std::string pointcloud_file;
  std::string target_srs_string;
  bool compressed;
  // Output
  std::string out_prefix;
  Options() : compressed(false){}
};

void handle_arguments( int argc, char *argv[], Options& opt ) {

  po::options_description general_options("General Options");
  general_options.add_options()
    ("compressed,c", po::bool_switch(&opt.compressed)->default_value(false)->implicit_value(true),
     "Compress using laszip.")
    ("output-prefix,o", po::value(&opt.out_prefix), "Specify the output prefix.")
    ("reference-spheroid,r", po::value(&opt.reference_spheroid),"Set the reference spheroid [Earth, Moon, Mars]. This will create a geo-referenced LAS file in respect to the spheroid. For Earth, the WGS84 datum is used.")
    ("t_srs", po::value(&opt.target_srs_string)->default_value(""),
     "Specify the output projection (PROJ.4 string).");

  general_options.add( asp::BaseOptionsDescription(opt) );

  po::options_description positional("");
  positional.add_options()
    ("input-file", po::value(&opt.pointcloud_file), "Input Point Cloud");

  po::positional_options_description positional_desc;
  positional_desc.add("input-file", 1);

  std::string usage("[options] <point-cloud>");
  bool allow_unregistered = false;
  std::vector<std::string> unregistered;
  po::variables_map vm =
    asp::check_command_line( argc, argv, opt, general_options, general_options,
                             positional, positional_desc, usage,
                             allow_unregistered, unregistered );

  if ( opt.pointcloud_file.empty() )
    vw_throw( ArgumentErr() << "Missing point cloud.\n"
              << usage << general_options );

  if ( opt.out_prefix.empty() )
    opt.out_prefix =
      asp::prefix_from_filename( opt.pointcloud_file );

  boost::to_lower( opt.reference_spheroid );

  // Create the output directory
  asp::create_out_dir(opt.out_prefix);

  // Turn on logging to file
  asp::log_to_file(argc, argv, "", opt.out_prefix);

}

int main( int argc, char *argv[] ) {

  // To do: need to understand what is the optimal strategy for
  // traversing the input point cloud file to minimize the reading
  // time.

  Options opt;
  try {
    handle_arguments( argc, argv, opt );

    // Save the las file in respect to a reference spheroid if provided
    // by the user.
    liblas::Header header;
    cartography::Datum datum;
    bool have_user_datum
      = asp::read_user_datum(0, 0, opt.reference_spheroid, datum);

    bool is_geodetic = false;
    if (have_user_datum || !opt.target_srs_string.empty()){

      cartography::GeoReference georef;
      // Set the srs string into georef.
      asp::set_srs_string(opt.target_srs_string,
                          have_user_datum, datum, georef);
      liblas::SpatialReference ref;
      std::string target_srs = georef.proj4_str();
      target_srs += " " + georef.datum().proj4_str();

      ref.SetFromUserInput(target_srs);
      vw_out() << "Using projection string: '" << target_srs << "'"<< std::endl;
      header.SetSRS(ref);

      is_geodetic = true;
      datum = georef.datum();
    }

    ImageViewRef<Vector3> point_image = asp::read_cloud<3>(opt.pointcloud_file);
    if (is_geodetic)
      point_image = cartesian_to_geodetic(point_image, datum);

    BBox3 cloud_bbox = asp::pointcloud_bbox(point_image, is_geodetic);

    // The las format stores the values as 32 bit integers. So, for a
    // given point, we store round((point-offset)/scale), as well as
    // the offset and scale values. Here we decide the values for
    // offset and scale to lose minimum amount of precision. We make
    // the scale almost as large as it can be without causing integer
    // overflow.
    Vector3 offset = (cloud_bbox.min() + cloud_bbox.max())/2.0;
    double maxInt = std::numeric_limits<int32>::max();
    maxInt *= 0.95; // Just in case stay a bit away
    Vector3 scale = cloud_bbox.size()/(2.0*maxInt);
    for (size_t i = 0; i < scale.size(); i++){
      if (scale[i] <= 0.0) scale[i] = 1.0e-16; // avoid degeneracy
    }

    // The line below causes trouble with compression in libLAS-1.7.0.
    //header.SetDataFormatId(liblas::ePointFormat1);
    header.SetScale(scale[0], scale[1], scale[2]);
    header.SetOffset(offset[0], offset[1], offset[2]);

    std::string lasFile;
    header.SetCompressed(opt.compressed);
    if (opt.compressed)
      lasFile = opt.out_prefix + ".laz";
    else
      lasFile = opt.out_prefix + ".las";

    vw_out() << "Writing LAS file: " << lasFile + "\n";
    std::ofstream ofs;
    ofs.open(lasFile.c_str(), std::ios::out | std::ios::binary);
    liblas::Writer writer(ofs, header);

    TerminalProgressCallback tpc("asp", "\t--> ");
    for (int row = 0; row < point_image.rows(); row++){
      tpc.report_fractional_progress(row, point_image.rows());
      for (int col = 0; col < point_image.cols(); col++){

        Vector3 point = point_image(col, row);

        // Skip no-data points
        bool is_good = ( (!is_geodetic && point != vw::Vector3()) ||
                         (is_geodetic  && !boost::math::isnan(point.z())) );
        if (!is_good) continue;

#if 0
        // For comparison later with las2txt.
        std::cout.precision(16);
        std::cout << "\npoint " << point[0] << ' ' << point[1] << ' '
                  << point[2] << std::endl;
#endif

        liblas::Point las_point(&header);
        las_point.SetCoordinates(point[0], point[1], point[2]);
        writer.WritePoint(las_point);

      }
    }
    tpc.report_finished();

  } ASP_STANDARD_CATCHES;

  return 0;
}
