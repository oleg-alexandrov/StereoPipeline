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

// Find the median disparity in each row

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
  std::string input_disp, output_disp; 
};

void handle_arguments(int argc, char *argv[], Options& opt) {
  
  po::options_description general_options("");
  general_options.add_options()
    ("input-disp",  po::value(&opt.input_disp)->default_value(""),
     "The input disparity.")
    ("output-disp",  po::value(&opt.output_disp)->default_value(""),
     "The output disparity.");

  general_options.add(vw::cartography::GdalWriteOptionsDescription(opt));
  
  po::options_description positional("");
  po::positional_options_description positional_desc;
  
  std::string usage("[options]");
  bool allow_unregistered = false;
  std::vector<std::string> unregistered;
  po::variables_map vm =
    asp::check_command_line(argc, argv, opt, general_options, general_options,
                            positional, positional_desc, usage,
                             allow_unregistered, unregistered);
  
  if (opt.input_disp == "") 
    vw::vw_throw( vw::ArgumentErr() << "The input disp is required.\n\n"
                  << usage << general_options);
  
  if (opt.output_disp == "") 
    vw::vw_throw( vw::ArgumentErr() << "The output disp is required.\n\n"
                  << usage << general_options);

  // Create the directory in which the output image will be written.
  vw::create_out_dir(opt.output_disp);
  
}

int main(int argc, char* argv[]) {
  
  Options opt;
  try {
    xercesc::XMLPlatformUtils::Initialize();

    handle_arguments(argc, argv, opt);

    std::cout << "Reading " << opt.input_disp << std::endl;
    DiskImageView< PixelMask<Vector2f> > in_disp(opt.input_disp);

    int cols = in_disp.cols();
    int rows = in_disp.rows();
    
    ImageView< PixelMask<Vector2f> > out_disp;
    out_disp.set_size(cols, rows);

    std::cout << "Num of cols and rows " << in_disp.cols() << ' ' << in_disp.rows() << std::endl;

    vw::TerminalProgressCallback tpc("asp", "\t--> ");
    double inc_amount = 1.0 / double(rows);
    tpc.report_progress(0);

    for (int row = 0; row < rows; row++) {

      std::vector<double> dx, dy;

      for (int col = 0; col < cols; col++) {
        PixelMask<Vector2f> pix = in_disp(col, row);
        if (!is_valid(pix)) continue;
        dx.push_back(pix.child().x());
        dy.push_back(pix.child().y());
      }


      if (dx.empty() || dy.empty()) {
        for (int col = 0; col < cols; col++) {
          out_disp(col, row) = Vector2f(0, 0);
          out_disp(col, row).invalidate();
        }
      }else{
        double x = dx[dx.size()/2];
        double y = dy[dy.size()/2];
        
        for (int col = 0; col < cols; col++) {
          out_disp(col, row) = Vector2f(x, y);
          out_disp(col, row).validate();
        }
      }
      
      tpc.report_incremental_progress(inc_amount);
    }
    tpc.report_finished();
    
    std::cout << "Writing: " << opt.output_disp << std::endl;

    bool has_georef = false, has_nodata = false;
    vw::cartography::GeoReference georef;
    double nodata = 0;
    
    block_write_gdal_image(opt.output_disp, out_disp, has_georef, georef, has_nodata, nodata,
                           opt, tpc);
    
    xercesc::XMLPlatformUtils::Terminate();
  } ASP_STANDARD_CATCHES;
}
