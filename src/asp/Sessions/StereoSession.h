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


/// \file StereoSession.h
///

#ifndef __STEREO_SESSION_H__
#define __STEREO_SESSION_H__

#include <vw/Image/ImageViewBase.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Camera/CameraModel.h>

#include <vw/Math/Functors.h>
#include <vw/Math/Geometry.h>
#include <vw/InterestPoint.h>

#include <boost/shared_ptr.hpp>
#include <boost/filesystem/operations.hpp>

#include <asp/Core/Common.h>
#include <asp/Core/InterestPointMatching.h>

namespace asp {

  // Specialize CompositionTransform that allows passing of BBox so
  // Map2CamTrans can cache itself.
  template <class Tx1T, class Tx2T>
  class CompositionTransformPassBBox : public vw::TransformBase<CompositionTransformPassBBox<Tx1T,Tx2T> > {
  public:
    CompositionTransformPassBBox( Tx1T const& tx1, Tx2T const& tx2 ) : tx1(tx1), tx2(tx2) {}

    Tx1T tx1; // Be sure to copy!
    Tx2T tx2; // public so that we can invoke caching manually for Map2CamTrans

    inline vw::Vector2 forward( vw::Vector2 const& p ) const { return tx1.forward( tx2.forward( p ) ); }
    inline vw::Vector2 reverse( vw::Vector2 const& p ) const { return tx2.reverse( tx1.reverse( p ) ); }

    inline vw::BBox2i reverse_bbox( vw::BBox2i const& bbox ) const {
      return this->tx2.reverse_bbox(this->tx1.reverse_bbox(bbox));
    }
  };

  // Stereo Sessions define for different missions or satellites how to:
  //   * Initialize, normalize, and align the input imagery
  //   * Extract the camera model
  //   * Custom code needed for correlation, filtering, and triangulation.
  //
  // All Stereo Session childs must define the following which are not
  // available in the the parent:
  //   typedef VWTransform left_tx_type;
  //   left_tx_type tx_left( void ) const;
  //   typedef VWTransform right_tx_type;
  //   right_tx_type tx_right( void ) const;
  //   typedef VWStereoModel stereo_model_type;

  class StereoSession {
  protected:
    asp::BaseOptions m_options;
    std::string m_left_image_file, m_right_image_file,
      m_left_camera_file, m_right_camera_file, m_out_prefix;
    std::string m_input_dem, m_extra_argument1,
      m_extra_argument2, m_extra_argument3;

    // This should probably be factored into the greater StereoSession
    // as ISIS Session does something very similar to this.
    template <class ViewT>
    vw::Vector4f gather_stats( vw::ImageViewBase<ViewT> const& view_base,
                               std::string const& tag) {
      using namespace vw;
      vw_out(InfoMessage) << "\t--> Computing statistics for the "+tag+" image\n";
      ViewT image = view_base.impl();
      int stat_scale = int(ceil(sqrt(float(image.cols())*float(image.rows()) / 1000000)));

      ChannelAccumulator<vw::math::CDFAccumulator<float> > accumulator;
      for_each_pixel( subsample( edge_extend(image, ConstantEdgeExtension()),
                                 stat_scale ), accumulator );
      Vector4f result( accumulator.quantile(0),
                       accumulator.quantile(1),
                       accumulator.approximate_mean(),
                       accumulator.approximate_stddev() );
      vw_out(InfoMessage) << "\t  " << tag << ": [ lo: " << result[0] << " hi: " << result[1]
                          << " m: " << result[2] << " s: " << result[3] << " ]\n";
      return result;
    }

    virtual void initialize (BaseOptions const& options,
                             std::string const& left_image_file,
                             std::string const& right_image_file,
                             std::string const& left_camera_file,
                             std::string const& right_camera_file,
                             std::string const& out_prefix,
                             std::string const& input_dem,
                             std::string const& extra_argument1,
                             std::string const& extra_argument2,
                             std::string const& extra_argument3);

  public:
    virtual ~StereoSession() {}

    // Methods for registering and creating stereo sessions.
    static StereoSession* create( std::string & session_type, // in-out variable
                                  BaseOptions const& options,
                                  std::string const& left_image_file = "",
                                  std::string const& right_image_file = "",
                                  std::string const& left_camera_file = "",
                                  std::string const& right_camera_file = "",
                                  std::string const& out_prefix = "",
                                  std::string const& input_dem = "",
                                  std::string const& extra_argument1 = "",
                                  std::string const& extra_argument2 = "",
                                  std::string const& extra_argument3 = "");
    typedef StereoSession* (*construct_func)();
    static void register_session_type( std::string const& id, construct_func func);

    // Helper function that retrieves both cameras.
    virtual void camera_models(boost::shared_ptr<vw::camera::CameraModel> &cam1,
                               boost::shared_ptr<vw::camera::CameraModel> &cam2);

    // Method that produces a Camera Model.
    virtual boost::shared_ptr<vw::camera::CameraModel>
    camera_model(std::string const& image_file,
                 std::string const& camera_file = "") = 0;

    // Method to help determine what session we actually have
    virtual std::string name() const = 0;

    virtual bool ip_matching( std::string const& match_filename,
                              double left_nodata_value,
                              double right_nodata_value );
    
    // Stage 1: Preprocessing
    //
    // Pre file is a pair of images.            ( ImageView<PixelT> )
    // Post file is a grayscale images.         ( ImageView<PixelGray<float> > )
    virtual void pre_preprocessing_hook(std::string const& input_file1,
                                        std::string const& input_file2,
                                        std::string &output_file1,
                                        std::string &output_file2);
    virtual void post_preprocessing_hook(std::string const& input_file1,
                                         std::string const& input_file2,
                                         std::string &output_file1,
                                         std::string &output_file2);

    // Stage 2: Correlation
    //
    // Pre file is a pair of grayscale images.  ( ImageView<PixelGray<float> > )
    // Post file is a disparity map.            ( ImageView<PixelDisparity> )
    virtual void pre_correlation_hook(std::string const& input_file1,
                                      std::string const& input_file2,
                                      std::string &output_file1,
                                      std::string &output_file2);
    virtual void post_correlation_hook(std::string const& input_file,
                                       std::string & output_file);

    // Stage 3: Filtering
    //
    // Pre file is a disparity map.  ( ImageView<PixelDisparity<float> > )
    // Post file is a disparity map. ( ImageView<PixelDisparity<float> > )
    virtual void pre_filtering_hook(std::string const& input_file,
                                    std::string & output_file);
    virtual void post_filtering_hook(std::string const& input_file,
                                     std::string & output_file);

    // Stage 4: Point cloud generation
    //
    // Pre file is a disparity map.  ( ImageView<PixelDisparity<float> > )
    // Post file is point image.     ( ImageView<Vector3> )
    virtual vw::ImageViewRef<vw::PixelMask<vw::Vector2f> >
    pre_pointcloud_hook(std::string const& input_file);
    virtual void post_pointcloud_hook(std::string const& input_file,
                                      std::string & output_file);

    void get_nodata_values(boost::shared_ptr<vw::DiskImageResource> left_rsrc,
                           boost::shared_ptr<vw::DiskImageResource> right_rsrc,
                           float & left_nodata_value,
                           float & right_nodata_value);
  };

} // end namespace asp

#endif // __STEREO_SESSION_H__
