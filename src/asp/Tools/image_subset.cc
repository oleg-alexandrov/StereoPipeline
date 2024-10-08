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

// Given a set of overlapping georeferenced images, find a small
// subset with almost the same coverage. This is used in SfS.

#include <asp/Core/Common.h>
#include <asp/Core/Macros.h>

#include <vw/Cartography/GeoTransform.h>
#include <vw/Core/Stopwatch.h>

#include <vector>
#include <string>
#include <iostream>
#include <limits>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

namespace fs = boost::filesystem;
namespace po = boost::program_options;

struct Options : public vw::GdalWriteOptions {
  std::string image_list_file, out_list; 
  double threshold;
};

void handle_arguments(int argc, char *argv[], Options& opt) {
  
 double nan = std::numeric_limits<double>::quiet_NaN();

  po::options_description general_options("");
  general_options.add_options()
   ("image-list", po::value(&opt.image_list_file)->default_value(""),
    "The list of input images.")
   ("output-list,o", po::value(&opt.out_list)->default_value(""),
    "The file having the produced image subset, and for each image the number "
    "of contributing pixels (sorted in decreasing order of contribution).")
   ("threshold", po::value(&opt.threshold)->default_value(nan),
    "The image threshold. Pixels no less than this will contribute to the coverage.")
  ;

  general_options.add(vw::GdalWriteOptionsDescription(opt));
  po::options_description positional("");
  po::positional_options_description positional_desc;

  std::string usage("[options]");
  bool allow_unregistered = false;
  std::vector<std::string> unregistered;
  po::variables_map vm =
    asp::check_command_line(argc, argv, opt, general_options, general_options,
                            positional, positional_desc, usage,
                            allow_unregistered, unregistered);

  // Validate inputs
  if (opt.image_list_file == "")
    vw::vw_throw(vw::ArgumentErr() << "Missing image list.\n");
  if (opt.out_list == "")
    vw::vw_throw(vw::ArgumentErr() << "Missing output list.\n");
  if (std::isnan(opt.threshold))
    vw::vw_throw(vw::ArgumentErr() << "The threshold must be set.\n");
    
  // Create the output directory
  vw::create_out_dir(opt.out_list);

  // Turn on logging to file
  asp::log_to_file(argc, argv, "", opt.out_list);
    
  return;
}

// Calc the score of adding a new image. Apply the new value to the output image if
// requested.
double calc_score_apply(std::string const& image_file, 
                        double threshold, bool apply, 
                        vw::cartography::GeoReference const& out_georef,
                        // Output
                        vw::ImageView<vw::PixelMask<float>> & out_img) {
  
  // Read the image
  vw::DiskImageResourceGDAL in_rsrc(image_file);
  
  // Read the georef
  vw::cartography::GeoReference georef;
  bool is_good = read_georeference(georef, in_rsrc);
  if (!is_good)
    vw::vw_throw(vw::ArgumentErr() << "No georeference found in " << image_file << ".\n");
    
  // Read the no-data
  double nodata = -std::numeric_limits<float>::max();
  if (in_rsrc.has_nodata_read())
    nodata = in_rsrc.nodata_read();
  
  // Read the image and mask the no-data
  vw::ImageViewRef<vw::PixelMask<float>> img 
    = create_mask(vw::DiskImageView<float>(in_rsrc), nodata);
  
  // Bounds of the output and input images
  vw::BBox2 out_box = vw::bounding_box(out_img);
  vw::BBox2 img_box = vw::bounding_box(img);
  
  // Form the geo transform from the output to the input image
  vw::cartography::GeoTransform geotrans(out_georef, georef,
                                         out_box, img_box);
  
  // Find the current image bounding box in the output image coordinates
  vw::BBox2 trans_box = geotrans.reverse_bbox(vw::bounding_box(img));
  // Grow to int
  trans_box = vw::grow_bbox_to_int(trans_box);
  
  // trans_box must be contained within out_box, as that's how out_box was formed
  if (!out_box.contains(trans_box))
    vw::vw_throw(vw::ArgumentErr()
                  << "An image bounding box is not contained in output box.\n"); 

  // Prepare the image for interpolation   
  vw::PixelMask<float> no_data;
  no_data.invalidate();
  typedef vw::ValueEdgeExtension<vw::PixelMask<float>> NoDataType;
  NoDataType no_data_ext(no_data);
  vw::InterpolationView<vw::EdgeExtensionView<vw::ImageViewRef<vw::PixelMask<float>>, NoDataType>, vw::BilinearInterpolation> interp_img 
    = interpolate(img, vw::BilinearInterpolation(), no_data_ext);
 
  // Iterate over the transformed box and calculate the score or apply the
  // image.
  // Note: Multi-threading slowed things down. Likely because ASP does not
  // likely reading pixels from same large images in multiple threads.
  double score = 0.0;
  for (int col = trans_box.min().x(); col < trans_box.max().x(); col++) {
    for (int row = trans_box.min().y(); row < trans_box.max().y(); row++) {
      
      vw::Vector2 out_pix(col, row);
      
      // If the value of out_pix is valid and no less than the threshold,
      // the new image will not add anything to the coverage
      vw::PixelMask<float> out_val = out_img(col, row);
      if (is_valid(out_val) && out_val.child() >= threshold)
        continue;
        
      // Find the interpolated pixel value
      vw::Vector2 img_pix = geotrans.forward(out_pix);
      vw::PixelMask<float> val = interp_img(img_pix.x(), img_pix.y());
      
      // Skip invalid values and values less than the threshold
      if (!is_valid(val) || val.child() < threshold)
        continue;

      // This is a good value, increment the score
      score++;
      
      // Apply the value if requested
      if (apply)
        out_img(col, row) = val;
    }  
  }
   
  return score;
}
                           
void run_image_subset(Options const& opt) {
  
  vw::vw_out() << "Reading: " << opt.image_list_file << "\n";
  std::vector<std::string> image_files;
  asp::read_list(opt.image_list_file, image_files);
  
  vw::vw_out() << "Reading input georeferences.\n";
  vw::cartography::GeoReference out_georef;

  // Loop through all input images
  int num_images = image_files.size();  
  vw::BBox2 out_box;
  for (int image_iter = 0; image_iter < num_images; image_iter++) {

    std::string image_file = image_files[image_iter]; 
    vw::DiskImageResourceGDAL in_rsrc(image_file);
    vw::DiskImageView<float> img(in_rsrc);
    vw::cartography::GeoReference georef;
    bool is_good = read_georeference(georef, in_rsrc);
    if (!is_good)
      vw::vw_throw(vw::ArgumentErr() << "No georeference found in " << image_file << ".\n");
    
    // Input image bounding box
    vw::BBox2 img_box = vw::bounding_box(img);

    // Borrow the geo info from the first image
    if (image_iter == 0) {
      out_georef = georef;
      out_box = img_box;
    }
    
    // The georef transform from current to output image
    vw::cartography::GeoTransform geotrans(out_georef, georef,
                                           out_box, img_box);

    // Convert the bounding box of current image to output pixel coordinates
    vw::BBox2 trans_img_box = geotrans.reverse_bbox(img_box);
    out_box.grow(trans_img_box);
  } // End loop through DEM files

  // grow this to int
  out_box = vw::grow_bbox_to_int(out_box);
  
  // Crop the georef to the output box. Adjust the box.
  out_georef = crop(out_georef, out_box);
  out_box.max() -= out_box.min();
  out_box.min() = vw::Vector2(0, 0);
  
  // Create the output image as a pixel mask, with all pixels being float and invalid
  vw::ImageView<vw::PixelMask<float>> out_img(out_box.width(), out_box.height());
  // invalidate all pixels
  for (int col = 0; col < out_img.cols(); col++) {
    for (int row = 0; row < out_img.rows(); row++) {
      out_img(col, row).invalidate();
    }
  }
  
  // Now process the images  
  vw::vw_out() << "Processing the images.\n";
  vw::Stopwatch sw;
  sw.start();
  vw::TerminalProgressCallback tpc("", "\t--> ");
  tpc.report_progress(0);
  double inc_amount = 2.0 / (double(num_images) * double(num_images));
  
  // Do as many passes as images
  std::map<std::string, double> inspected;
  for (int pass = 0; pass < num_images; pass++) {
        
    // Skip inspected images
    if (inspected.find(image_files[pass]) != inspected.end())
      continue;

    double best_score = -1.0;
    int best_index = 0;
    
    // Inner iteration over all images. Skip the inspected ones.
    for (int inner = 0; inner < num_images; inner++) {
     
      if (inspected.find(image_files[inner]) != inspected.end())
        continue;

      bool apply = false;
      double score = calc_score_apply(image_files[inner], opt.threshold, apply, 
                                      out_georef, out_img);
      
      // TODO(oalexan1): If the score is zero, declare this image as not 
      // useful, as its score will continue to be 0 for future out_img,
      // since those only grow. This must be tested though.
        
      // Update the best score
      if (score > best_score) {
        best_score = score;
        best_index = inner;
      }
      tpc.report_incremental_progress(inc_amount);
    }

    // If the score is 0, we are done, as more passes won't help
    if (best_score == 0.0) {
      break;
    }
    
    // Add the contribution of the best image
    bool apply = true;
    calc_score_apply(image_files[best_index], opt.threshold, apply, 
                     out_georef, out_img);
    
    // Flag the current image as inspected and store its score
    inspected[image_files[best_index]] = best_score;
  }  

  tpc.report_finished();
  sw.stop();
  vw::vw_out() << "Elapsed time: " << sw.elapsed_seconds() << " s.\n";
  
  // Save in decreasing order of score
  std::vector<std::pair<double, std::string>> sorted;
  for (auto const& p: inspected)
    sorted.push_back(std::make_pair(p.second, p.first));
  std::sort(sorted.begin(), sorted.end());
  vw::vw_out() << "Writing: " << opt.out_list << "\n";
  std::ofstream ofs(opt.out_list.c_str());
  for (int i = sorted.size()-1; i >= 0; i--)
    ofs << sorted[i].second << " " << sorted[i].first << "\n";
  ofs.close();
  
} // End function run_image_subset()

int main(int argc, char *argv[]) {

  Options opt;
  try {
    handle_arguments(argc, argv, opt);
    run_image_subset(opt);
  } ASP_STANDARD_CATCHES;

  return 0;
}
