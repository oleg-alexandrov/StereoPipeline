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


/// \file stereo_tri.cc
///
//#define USE_GRAPHICS

#include <asp/Tools/stereo.h>
#include <asp/Sessions/RPC/RPCModel.h>
#include <asp/Sessions/RPC/RPCStereoModel.h>
#include <vw/Cartography.h>
#include <vw/Camera/CameraModel.h>
#include <vw/Camera/Extrinsics.h>
#include <vw/Stereo/StereoView.h>
#include <vw/Camera/PinholeModel.h>
#include <asp/Sessions/DG/LinescanDGModel.h>
#include <asp/Sessions/DG/StereoSessionDG.h>

using namespace vw;
using namespace asp;

vw::camera::CameraModel * g_cam1;
vw::camera::CameraModel * g_cam2;

namespace vw {
  template<> struct PixelFormatID<PixelMask<Vector<float, 5> > >   { static const PixelFormatEnum value = VW_PIXEL_GENERIC_6_CHANNEL; };
  template<> struct PixelFormatID<Vector<double, 4> >   { static const PixelFormatEnum value = VW_PIXEL_GENERIC_4_CHANNEL; };
  template<> struct PixelFormatID<Vector<float, 2> > { static const PixelFormatEnum value = VW_PIXEL_GENERIC_2_CHANNEL; };
}

// Class definition
template <class DisparityImageT, class StereoModelT>
class StereoAndErrorView : public ImageViewBase<StereoAndErrorView<DisparityImageT, StereoModelT> >
{

  vw::camera::CameraModel const* m_camera_model1;
  vw::camera::CameraModel const* m_camera_model2;
  DisparityImageT m_disparity_map;
  StereoModelT m_stereo_model;
  typedef typename DisparityImageT::pixel_type dpixel_type;

  template <class PixelT>
  struct NotSingleChannel {
    static const bool value = (1 != CompoundNumChannels<typename UnmaskedPixelType<PixelT>::type>::value);
  };

  template <class T>
  inline typename boost::enable_if<IsScalar<T>,Vector3>::type
  StereoModelHelper( StereoModelT const& model, Vector2 const& index,
                     T const& disparity, double& error ) const {
    std::cout << "helper1" << std::endl;
    return model( index, Vector2( index[0] + disparity, index[1] ), error );
  }

  template <class T>
  inline typename boost::enable_if_c<IsCompound<T>::value && (CompoundNumChannels<typename UnmaskedPixelType<T>::type>::value == 1),Vector3>::type
  StereoModelHelper( StereoModelT const& model, Vector2 const& index,
                     T const& disparity, double& error ) const {
    std::cout << "helper2" << std::endl;
    return model( index, Vector2( index[0] + disparity, index[1] ), error );
  }

  template <class T>
  inline typename boost::enable_if_c<IsCompound<T>::value && (CompoundNumChannels<typename UnmaskedPixelType<T>::type>::value != 1),Vector3>::type
  StereoModelHelper( StereoModelT const& model, Vector2 const& index,
                     T const& disparity, double& error ) const {

    Vector3 res = model( index, Vector2( index[0] + disparity[0],
                                         index[1] + disparity[1] ), error );
    
    if (index[0] == 701 && index[1] == 701 ){

      std::cout.precision(20);
      std::cout << "\n\n\n----------pixel is " << index[0] << ' ' << index[1]  << ' ' << index[0] + disparity[0] << ' ' << index[1] + disparity[1]  << "\n\n\n" << std::endl;

      std::cout << "result is " << res << std::endl;
      std::cout << "error is: " << error << std::endl;

      std::cout << "point to pixel1: " << g_cam1->point_to_pixel(res) << std::endl;
      std::cout << "point to pixel2: " << g_cam2->point_to_pixel(res) << std::endl;
      
#if 0
      typedef LinescanDGModel<camera::PiecewiseAPositionInterpolation, camera::SLERPPoseInterpolation, camera::TLCTimeInterpolation> camera_type;

      camera_type const* cam1 = dynamic_cast<camera_type const*>(g_cam1);
      camera_type const* cam2 = dynamic_cast<camera_type const*>(g_cam2);
      std::cout << "\nbefore casting pointers: " << g_cam1 << ' ' << g_cam2 << std::endl;
      std::cout << "after casting pointers: " << cam1 << ' ' << cam2 << std::endl;

      Vector3 res1 = model( index, Vector2( index[0] + disparity[0],
                                            index[1] + disparity[1] ), error );
          
      camera::PinholeModel pinhole1 = cam1->linescan_to_pinhole(index[1]);
      camera::PinholeModel pinhole2 = cam2->linescan_to_pinhole(index[1] + disparity[1]);
          
      bool least_squares_refine = false;
      stereo::StereoModel model2(&pinhole1, &pinhole2, least_squares_refine);
          
      Vector3 res2 = model2( index, Vector2( index[0] + disparity[0],
                                             index[1] + disparity[1] ), error );
          
      std::cout << "result 1 is " << res1 << std::endl;
      std::cout << "result 2 is " << res2 << std::endl;

      std::string cam1File = "cam1.pinhole";
      std::cout << "Writing " << cam1File << std::endl;
      pinhole1.write(cam1File);
      
      std::string cam2File = "cam2.pinhole";
      std::cout << "Writing " << cam2File << std::endl;
      pinhole2.write(cam2File);

      std::cout << "pinhole camera1 is:\n\n" << pinhole1 << std::endl << std::endl;
      std::cout << "pinhole camera1 matrix is " << pinhole1.camera_matrix() << std::endl;
      std::cout << "pinhole camera2 is:\n\n" << pinhole2 << std::endl << std::endl;
      std::cout << "pinhole camera2 matrix is " << pinhole2.camera_matrix() << std::endl;
#endif

    }    
    return res;
  }
  
public:

  typedef Vector4 pixel_type;
  typedef const Vector4 result_type;
  typedef ProceduralPixelAccessor<StereoAndErrorView> pixel_accessor;

  StereoAndErrorView( DisparityImageT const& disparity_map,
                      vw::camera::CameraModel const* camera_model1,
                      vw::camera::CameraModel const* camera_model2,
                      bool least_squares_refine = false) :
    m_camera_model1(camera_model1),
    m_camera_model2(camera_model2),
    m_disparity_map(disparity_map),
    m_stereo_model(camera_model1, camera_model2, least_squares_refine) {}

  StereoAndErrorView( DisparityImageT const& disparity_map,
                      StereoModelT const& stereo_model) :
    m_disparity_map(disparity_map),
    m_stereo_model(stereo_model) {}

  inline int32 cols() const { return m_disparity_map.cols(); }
  inline int32 rows() const { return m_disparity_map.rows(); }
  inline int32 planes() const { return 1; }

  inline pixel_accessor origin() const { return pixel_accessor(*this); }

  inline result_type operator()( size_t i, size_t j, size_t p=0 ) const {
    if ( is_valid(m_disparity_map(i,j,p)) ) {
      pixel_type result;
      //std::cout << "--- in operator=" << std::endl;
      subvector(result,0,3) = StereoModelHelper( m_stereo_model, Vector2(i,j),
                                                 m_disparity_map(i,j,p), result[3] );
      return result;
    }
    // For missing pixels in the disparity map, we return a null 3D position.
    return pixel_type();
  }

  /// \cond INTERNAL
  typedef StereoAndErrorView<typename DisparityImageT::prerasterize_type, StereoModelT> prerasterize_type;
  inline prerasterize_type prerasterize( BBox2i const& bbox ) const { return prerasterize_type( m_disparity_map.prerasterize(bbox), m_stereo_model ); }
  template <class DestT> inline void rasterize( DestT const& dest, BBox2i const& bbox ) const { vw::rasterize( prerasterize(bbox), dest, bbox ); }
  /// \endcond
};

// Variant that uses LUT tables
template <class DisparityImageT, class LUTImage1T, class LUTImage2T, class StereoModelT>
class StereoLUTAndErrorView : public ImageViewBase<StereoLUTAndErrorView<DisparityImageT, LUTImage1T, LUTImage2T, StereoModelT> >
{
  vw::camera::CameraModel const* m_camera_model1;
  vw::camera::CameraModel const* m_camera_model2;
  DisparityImageT m_disparity_map;
  LUTImage1T m_lut_image1;
  InterpolationView< EdgeExtensionView<LUTImage2T, ConstantEdgeExtension>, BilinearInterpolation>  m_lut_image2;
  LUTImage2T m_lut_image2_org;
  StereoModelT m_stereo_model;
  typedef typename DisparityImageT::pixel_type dpixel_type;
  
  template <class PixelT>
  struct NotSingleChannel {
    static const bool value = (1 != CompoundNumChannels<typename UnmaskedPixelType<PixelT>::type>::value);
  };

  template <class T>
  inline typename boost::enable_if<IsScalar<T>,Vector3>::type
  StereoModelHelper( size_t i, size_t j, T const& disparity, double& error ) const {
    std::cout << "helper4" << std::endl;

    return m_stereo_model( m_lut_image1(i,j),
                           m_lut_image2( T(i) + disparity, j ), error );
  }

  template <class T>
  inline typename boost::enable_if_c<IsCompound<T>::value && (CompoundNumChannels<typename UnmaskedPixelType<T>::type>::value == 1),Vector3>::type
  StereoModelHelper( size_t i, size_t j, T const& disparity, double& error ) const {
    std::cout << "helper5" << std::endl;
    return m_stereo_model( m_lut_image1(i,j),
                           m_lut_image2( float(i) + disparity, j ),  error );
  }

  template <class T>
  inline typename boost::enable_if_c<IsCompound<T>::value && (CompoundNumChannels<typename UnmaskedPixelType<T>::type>::value != 1),Vector3>::type
    StereoModelHelper( size_t i, size_t j, T const& disparity, double& error ) const {

    float i2 = float(i) + disparity[0];
    float j2 = float(j) + disparity[1];
    if ( i2 < 0 || i2 >= m_lut_image2.cols() ||
         j2 < 0 || j2 >= m_lut_image2.rows() ||
         !is_valid( disparity ) ){
      return Vector3(); // out of bounds
    }

    return m_stereo_model( m_lut_image1(i,j), m_lut_image2(i2, j2), error );
  }

public:

  typedef Vector4 pixel_type;
  typedef const Vector4 result_type;
  typedef ProceduralPixelAccessor<StereoLUTAndErrorView> pixel_accessor;

  StereoLUTAndErrorView( ImageViewBase<DisparityImageT> const& disparity_map,
                         ImageViewBase<LUTImage1T> const& lut_image1,
                         ImageViewBase<LUTImage2T> const& lut_image2,
                         vw::camera::CameraModel const* camera_model1,
                         vw::camera::CameraModel const* camera_model2,
                         bool least_squares_refine = false) :
    m_camera_model1(camera_model1),
    m_camera_model2(camera_model2),
    m_disparity_map(disparity_map.impl()), m_lut_image1( lut_image1.impl() ),
    m_lut_image2(interpolate(lut_image2.impl())),
    m_lut_image2_org( lut_image2.impl() ),
    m_stereo_model(camera_model1, camera_model2, least_squares_refine) {
    std::cout << "pointers to camera are " << m_camera_model1 << ' ' << m_camera_model2 << std::endl;
  }

  StereoLUTAndErrorView( ImageViewBase<DisparityImageT> const& disparity_map,
                         ImageViewBase<LUTImage1T> const& lut_image1,
                         ImageViewBase<LUTImage2T> const& lut_image2,
                         StereoModelT const& stereo_model) :
    m_disparity_map(disparity_map.impl()), m_lut_image1(lut_image1.impl()),
    m_lut_image2(interpolate(lut_image2.impl())),
    m_lut_image2_org( lut_image2.impl() ),
    m_stereo_model(stereo_model) {}

  inline int32 cols() const { return m_disparity_map.cols(); }
  inline int32 rows() const { return m_disparity_map.rows(); }
  inline int32 planes() const { return 1; }

  inline pixel_accessor origin() const { return pixel_accessor(*this); }

  inline result_type operator()( size_t i, size_t j, size_t p=0 ) const {
    if ( is_valid(m_disparity_map(i,j,p)) ) {
      pixel_type result;
      subvector(result,0,3) = StereoModelHelper( i, j, m_disparity_map(i,j,p), result[3] );
      return result;
    }
    // For missing pixels in the disparity map, we return a null 3D position.
    return pixel_type();
  }

  /// \cond INTERNAL
  typedef StereoLUTAndErrorView<CropView<ImageView<typename DisparityImageT::pixel_type> >,
                                typename LUTImage1T::prerasterize_type,
                                typename LUTImage2T::prerasterize_type,
                                StereoModelT> prerasterize_type;
  inline prerasterize_type prerasterize( BBox2i const& bbox ) const {
    typedef typename DisparityImageT::pixel_type DPixelT;
    CropView<ImageView<DPixelT> > disparity_preraster =
      crop( ImageView<DPixelT>( crop( m_disparity_map, bbox ) ),
            -bbox.min().x(), -bbox.min().y(), cols(), rows() );

    // Calculate the range of disparities in our BBox to determine the
    // crop size needed for LUT 2.
    typedef typename UnmaskedPixelType<DPixelT>::type accum_t;
    PixelAccumulator<EWMinMaxAccumulator<accum_t> > accumulator;
    for_each_pixel( disparity_preraster.child(), accumulator );

    BBox2i preraster(0,0,0,0);
    if ( accumulator.is_valid() ){
      accum_t input_min = accumulator.minimum();
      accum_t input_max = accumulator.maximum();
      // Bugfix: expand the preraster window by 1 pixel to avoid segfaults.
      // This is needed for interpolation.
      preraster = BBox2i(bbox.min() + floor(Vector2f(input_min[0]-1,input_min[1]-1)),
                         bbox.max() + ceil(Vector2(input_max[0]+1,input_max[1]+1)) );
    }
    
    return prerasterize_type( disparity_preraster,
                              m_lut_image1.prerasterize(bbox),
                              m_lut_image2_org.prerasterize(preraster),
                              m_stereo_model );
  }
  template <class DestT> inline void rasterize( DestT const& dest, BBox2i const& bbox ) const { vw::rasterize( prerasterize(bbox), dest, bbox ); }
  /// \endcond
};

template <class ImageT, class StereoModelT>
StereoAndErrorView<ImageT, StereoModelT>
stereo_error_triangulate( ImageViewBase<ImageT> const& v,
                          vw::camera::CameraModel const* camera1,
                          vw::camera::CameraModel const* camera2 ) {
  std::cout << "now in triangulate!" << std::endl;
  return StereoAndErrorView<ImageT, StereoModelT>( v.impl(), camera1, camera2 );
}

template <class ImageT, class StereoModelT>
StereoAndErrorView<ImageT, StereoModelT>
lsq_stereo_error_triangulate( ImageViewBase<ImageT> const& v,
                              vw::camera::CameraModel const* camera1,
                              vw::camera::CameraModel const* camera2 ) {
  return StereoAndErrorView<ImageT, StereoModelT>( v.impl(), camera1, camera2, true );
}

template <class DisparityT, class LUT1T, class LUT2T, class StereoModelT>
StereoLUTAndErrorView<DisparityT, LUT1T, LUT2T, StereoModelT>
stereo_error_triangulate( ImageViewBase<DisparityT> const& disparity,
                          ImageViewBase<LUT1T> const& lut1,
                          ImageViewBase<LUT2T> const& lut2,
                          vw::camera::CameraModel const* camera1,
                          vw::camera::CameraModel const* camera2 ) {
  typedef StereoLUTAndErrorView<DisparityT, LUT1T, LUT2T, StereoModelT> result_type;
  return result_type( disparity.impl(), lut1.impl(), lut2.impl(),
                      camera1, camera2 );
}

template <class DisparityT, class LUT1T, class LUT2T, class StereoModelT>
StereoLUTAndErrorView<DisparityT, LUT1T, LUT2T, StereoModelT>
lsq_stereo_error_triangulate( ImageViewBase<DisparityT> const& disparity,
                              ImageViewBase<LUT1T> const& lut1,
                              ImageViewBase<LUT2T> const& lut2,
                              vw::camera::CameraModel const* camera1,
                              vw::camera::CameraModel const* camera2 ) {
  typedef StereoLUTAndErrorView<DisparityT,LUT1T,LUT2T, StereoModelT> result_type;
  return result_type( disparity.impl(), lut1.impl(), lut2.impl(),
                      camera1, camera2, true );
}

template <class StereoModelT>
void stereo_triangulation( Options const& opt ) {
  vw_out() << "\n[ " << current_posix_time_string()
           << " ] : Stage 4 --> TRIANGULATION \n";

  typedef ImageViewRef<PixelMask<Vector2f> > PVImageT;
  typedef ImageViewRef<Vector2f>             VImageT;
  try {
    PVImageT disparity_map =
      opt.session->pre_pointcloud_hook(opt.out_prefix+"-F.tif");

    boost::shared_ptr<camera::CameraModel> camera_model1, camera_model2;
    opt.session->camera_models(camera_model1, camera_model2);
    
#if HAVE_PKG_VW_BUNDLEADJUSTMENT
    // If the user has generated a set of position and pose
    // corrections using the bundle_adjust program, we read them in
    // here and incorporate them into our camera model.
    Vector3 position_correction;
    Quaternion<double> pose_correction;
    if (fs::exists(fs::path(in_file1).replace_extension("adjust"))) {
      read_adjustments(fs::path(in_file1).replace_extension("adjust").string(),
                       position_correction, pose_correction);
      camera_model1 =
        boost::shared_ptr<CameraModel>(new AdjustedCameraModel(camera_model1,
                                                               position_correction,
                                                               pose_correction));
    }
    if (fs::exists(fs::path(in_file2).replace_extension("adjust"))) {
      read_adjustments(fs::path(in_file2).replace_extension("adjust").string(),
                       position_correction, pose_correction);
      camera_model2 =
        boost::shared_ptr<CameraModel>(new AdjustedCameraModel(camera_model2,
                                                               position_correction,
                                                               pose_correction));
    }
#endif

    
#if 1
    g_cam1 = (vw::camera::CameraModel*) camera_model1.get();
    g_cam2 = (vw::camera::CameraModel*) camera_model2.get();
    std::cout << "----camera 1 and 2 are " << g_cam1 << ' ' << g_cam2 << std::endl;
#endif
    
    // If the distance from the left camera center to a point is
    // greater than the universe radius, we remove that pixel and
    // replace it with a zero vector, which is the missing pixel value
    // in the point_image.
    //
    // We apply the universe radius here and then write the result
    // directly to a file on disk.
    stereo::UniverseRadiusFunc universe_radius_func(Vector3(),0,0);
    if ( stereo_settings().universe_center == "camera" ) {

      if (opt.stereo_session_string == "rpc")
        vw_throw(InputErr() << "Stereo with RPC cameras cannot have the camera as the universe center.\n");
      
      universe_radius_func =
        stereo::UniverseRadiusFunc(camera_model1->camera_center(Vector2()),
                                   stereo_settings().near_universe_radius,
                                   stereo_settings().far_universe_radius);
      
    } else if ( stereo_settings().universe_center == "zero" ) {
      universe_radius_func =
        stereo::UniverseRadiusFunc(Vector3(),
                                   stereo_settings().near_universe_radius,
                                   stereo_settings().far_universe_radius);
    }

    
    // Apply radius function and stereo model in one go
    vw_out() << "\t--> Generating a 3D point cloud.   " << std::endl;
    ImageViewRef<Vector4> point_cloud;
    if ( opt.session->has_lut_images() ) {
      if ( stereo_settings().use_least_squares )
        point_cloud =
          per_pixel_filter(lsq_stereo_error_triangulate<PVImageT, VImageT, VImageT, StereoModelT>
                           ( disparity_map,
                             opt.session->lut_image_left(),
                             opt.session->lut_image_right(),
                             camera_model1.get(),
                             camera_model2.get() ),
                           universe_radius_func);
      else
        point_cloud =
          per_pixel_filter(stereo_error_triangulate<PVImageT, VImageT, VImageT, StereoModelT>
                           ( disparity_map,
                             opt.session->lut_image_left(),
                             opt.session->lut_image_right(),
                             camera_model1.get(),
                             camera_model2.get() ),
                           universe_radius_func);
    } else {
      if ( stereo_settings().use_least_squares )
        point_cloud =
          per_pixel_filter(lsq_stereo_error_triangulate<PVImageT, StereoModelT>( disparity_map,
                                                                                 camera_model1.get(),
                                                                                 camera_model2.get() ),
                           universe_radius_func);
      else
        point_cloud =
          per_pixel_filter(stereo_error_triangulate<PVImageT, StereoModelT>( disparity_map,
                                                                             camera_model1.get(),
                                                                             camera_model2.get() ),
                           universe_radius_func);
    }

    vw_out(VerboseDebugMessage,"asp") << "Writing Point Cloud: "
                                      << opt.out_prefix + "-PC.tif\n";

    DiskImageResource* rsrc =
      asp::build_gdal_rsrc( opt.out_prefix + "-PC.tif",
                            point_cloud, opt );
    if ( opt.stereo_session_string == "isis" )
      write_image(*rsrc, point_cloud,
                  TerminalProgressCallback("asp", "\t--> Triangulating: "));
    else
      block_write_image(*rsrc, point_cloud,
                        TerminalProgressCallback("asp", "\t--> Triangulating: "));
    delete rsrc;
    vw_out() << "\t--> " << universe_radius_func;

  } catch (IOErr const& e) {
    vw_throw( ArgumentErr() << "\nUnable to start at point cloud stage -- could not read input files.\n"
              << e.what() << "\nExiting.\n\n" );
  }
}

int main( int argc, char* argv[] ) {

  stereo_register_sessions();
  
  Options opt;
  try {
    handle_arguments( argc, argv, opt,
                      TriangulationDescription() );

    // Internal Processes
    //---------------------------------------------------------
    
    if (opt.stereo_session_string != "rpc"){
      stereo_triangulation<stereo::StereoModel>( opt );
    }else{
      // The RPC camera model does not fit in the existing framework
      // as its method of triangulation is quite different.
      stereo_triangulation<asp::RPCStereoModel>( opt );
    }

    vw_out() << "\n[ " << current_posix_time_string()
             << " ] : TRIANGULATION FINISHED \n";

  } ASP_STANDARD_CATCHES;

  return 0;
}
