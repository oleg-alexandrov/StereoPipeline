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


/// \file bundle_adjust.cc
///
// TODO(oalexan1):
// Make the camera_to_ground_dist a parameter.
// Create honest cameras using essential matrix.
// Study why the focal length turns out negative.
// Must make sure to use the Find3DAffineTransform only with >= 3 points!
// Must add unit test for Find3DAffineTransform()
// Must move a lot of code from .h to .cc
// Why convert to CameraRelationNetwork, use only the cnet, more clear that
// way.

#include <asp/Core/Macros.h>
#include <asp/Tools/bundle_adjust.h>
#include <asp/Sessions.h>
#include <ceres/ceres.h>
#include <ceres/loss_function.h>

namespace po = boost::program_options;
namespace fs = boost::filesystem;

using namespace vw;
using namespace vw::camera;
using namespace vw::ba;

std::string UNSPECIFIED_DATUM = "unspecified_datum";

Eigen::Affine3d Find3DAffineTransform(Eigen::Matrix3Xd in, Eigen::Matrix3Xd out) {
  // Default output
  Eigen::Affine3d A;
  A.linear() = Eigen::Matrix3d::Identity(3, 3);
  A.translation() = Eigen::Vector3d::Zero();

  if (in.cols() != out.cols())
    throw "Find3DAffineTransform(): input data mis-match";

  // First find the scale, by finding the ratio of sums of some distances,
  // then bring the datasets to the same scale.
  double dist_in = 0, dist_out = 0;
  for (int col = 0; col < in.cols()-1; col++) {
    dist_in  += (in.col(col+1) - in.col(col)).norm();
    dist_out += (out.col(col+1) - out.col(col)).norm();
  }
  if (dist_in <= 0 || dist_out <= 0)
    return A;
  double scale = dist_out/dist_in;
  out /= scale;

  // Find the centroids then shift to the origin
  Eigen::Vector3d in_ctr = Eigen::Vector3d::Zero();
  Eigen::Vector3d out_ctr = Eigen::Vector3d::Zero();
  for (int col = 0; col < in.cols(); col++) {
    in_ctr  += in.col(col);
    out_ctr += out.col(col);
  }
  in_ctr /= in.cols();
  out_ctr /= out.cols();
  for (int col = 0; col < in.cols(); col++) {
    in.col(col)  -= in_ctr;
    out.col(col) -= out_ctr;
  }

  // SVD
  Eigen::MatrixXd Cov = in * out.transpose();
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(Cov, Eigen::ComputeThinU | Eigen::ComputeThinV);

  // Find the rotation
  double d = (svd.matrixV() * svd.matrixU().transpose()).determinant();
  if (d > 0)
    d = 1.0;
  else
    d = -1.0;
  Eigen::Matrix3d I = Eigen::Matrix3d::Identity(3, 3);
  I(2, 2) = d;
  Eigen::Matrix3d R = svd.matrixV() * I * svd.matrixU().transpose();

  // The final transform
  A.linear() = scale * R;
  A.translation() = scale*(out_ctr - R*in_ctr);

  return A;
}

namespace asp{

  // Given a vector of strings, identify and store separately the the
  // list of GCPs. This should be useful for those programs who accept
  // their data in a mass input vector.
  std::vector<std::string>
  extract_gcps( std::vector<std::string>& image_files ) {
    std::vector<std::string> gcp_files;
    std::vector<std::string>::iterator it = image_files.begin();
    while ( it != image_files.end() ) {
      if ( boost::iends_with(boost::to_lower_copy(*it), ".gcp") ){
        gcp_files.push_back( *it );
        it = image_files.erase( it );
      } else
        it++;
    }

    return gcp_files;
  }

  bool images_are_cubes(std::vector<std::string>& image_files){
    bool are_cubes = true;
    for (int i = 0; i < (int)image_files.size(); i++){
      if ( ! boost::iends_with(boost::to_lower_copy(image_files[i]), ".cub") )
        are_cubes = false;
    }
    return are_cubes;
  }

} // end namespace asp

struct Options : public asp::BaseOptions {
  std::vector<std::string> image_files, camera_files, gcp_files;
  std::string cnet_file, out_prefix, stereo_session_string, cost_function, ba_type;

  double lambda, camera_weight, robust_threshold;
  int report_level, min_matches, max_iterations, overlap_limit;

  bool save_iteration, have_input_cams;
  std::string datum_str;
  double semi_major, semi_minor;

  boost::shared_ptr<ControlNetwork> cnet;
  std::vector<boost::shared_ptr<CameraModel> > camera_models;
  cartography::Datum datum;

  // Make sure all values are initialized, even though they will be
  // over-written later.
  Options():lambda(-1.0), camera_weight(0), robust_threshold(0), report_level(0), min_matches(0),
            max_iterations(0), overlap_limit(0), save_iteration(false), have_input_cams(true),
            semi_major(0), semi_minor(0),
            datum(cartography::Datum(UNSPECIFIED_DATUM, "User Specified Spheroid",
                                     "Reference Meridian", 1, 1, 0)){}
};

template<class ModelT>
typename boost::disable_if<boost::is_same<ModelT,BAPinholeModel>,void>::type
init_optimization_vars(ModelT & ba_model, Options & opt,
                          std::vector<double> & cameras_vec,
                          std::vector<double> & intrinsics_vec){
  // Do nothing, need this as part of the interface
}

template<class ModelT>
typename boost::enable_if<boost::is_same<ModelT,BAPinholeModel>,void>::type
init_optimization_vars(ModelT & ba_model, Options & opt,
                       std::vector<double> & cameras_vec,
                       std::vector<double> & intrinsics_vec){

  size_t num_camera_params    = ModelT::camera_params_n;
  size_t num_intrinsic_params = ModelT::intrinsic_params_n;

  // Pull info from the initial guess input cameras
  // we have created by now.
  for (int icam = 0; icam < (int)opt.camera_models.size(); icam++){
    vw::camera::PinholeModel2 * pincam
      = dynamic_cast<vw::camera::PinholeModel2*>(opt.camera_models[icam].get());
    VW_ASSERT(pincam != NULL,
              vw::ArgumentErr() << "A pinhole camera expected.\n");


    std::cout << "pinhole model2 is " << *pincam << std::endl;

    typename ModelT::camera_intr_vector_t b;
    ba_model.put_pinhole_model(*pincam, b);

    for (size_t q = 0; q < num_camera_params; q++)
      cameras_vec[icam*num_camera_params + q] = b[q];
    for (size_t q = 0; q < num_intrinsic_params; q++)
      intrinsics_vec[q] = b[num_camera_params + q];
  }
}

// Create pinhole models from scratch, and create initial guesses
// for triangulated points.
void update_cnet_and_init_cams(Options & opt, ControlNetwork& cnet){

  VW_ASSERT(!opt.have_input_cams,
            vw::ArgumentErr() << "Have input cameras, must not create new ones.\n");

#if 1

  for (int ipt = 0; ipt < (int)cnet.size(); ipt++){
    if (cnet[ipt].type() == ControlPoint::GroundControlPoint) {
      std::cout << "--gcp is " << cnet[ipt].position() << std::endl;
      continue;
    }
  }

  vw::Matrix3x3 R1, R2;
  vw::Vector3 t1, t2;
  double f, cx, cy;
  std::ifstream ifs("pinhole.txt");
  for (int row = 0; row < 3; row++) {
    for (int col = 0; col < 3; col++) {
      ifs >> R1(row, col);
    }
  }
  for (int row = 0; row < 3; row++) {
    ifs >> t1[row];
  }
  for (int row = 0; row < 3; row++) {
    for (int col = 0; col < 3; col++) {
      ifs >> R2(row, col);
    }
  }
  for (int row = 0; row < 3; row++) {
    ifs >> t2[row];
  }

  double scale = atof(getenv("S"));
  std::cout << "---scale is " << scale << std::endl;
  ifs >> f >> cx >> cy;

  f  *= scale;
  cx *= scale;
  cy *= scale;

  std::cout << "Left: " << R1 << "\n" << t1 << std::endl;
//   std::cout << "Right: " << R2 << "\n" << t2 << std::endl;
//   std::cout << "focal and cx, cy: " << f << ' ' << cx << ' ' << cy << std::endl;

//   vw::Matrix3x3 Q;
//   for (int i = 0; i < 3; i++) {
//     for (int j = 0; j < 3; j++) {
//       Q(i, j)=0;
//     }
//   }
//   Q(0, 0) = 1;
//   Q(1, 1) = -1;
//   Q(2, 2) = -1;
//   std::cout << "---temporary hack 4" << std::endl;

  for (int icam = 0; icam < (int)opt.camera_models.size(); icam++){

    vw::camera::PinholeModel2 * pincam
      = dynamic_cast<vw::camera::PinholeModel2*>(opt.camera_models[icam].get());
    VW_ASSERT(pincam != NULL,
              vw::ArgumentErr() << "A pinhole camera expected.\n");

    if (icam == 0)
      *pincam = vw::camera::PinholeModel2(-inverse(R1)*t1,
                                         inverse(R1),
                                         f, f,  // focal lengths. Why it turns out negative?
                                         cx, cy); // pixel offsets
    else
      *pincam = vw::camera::PinholeModel2(-inverse(R2)*t2,
                                         inverse(R2),
                                         f, f,  // focal lengths. Why it turns out negative?
                                         cx, cy); // pixel offsets

    std::cout << "pincam is " << *pincam << std::endl;

  }

  std::cout << "---test the logic here with gcp!!!!!!!" << std::endl;
  for (int ipt = 0; ipt < (int)cnet.size(); ipt++){
    if (cnet[ipt].type() == ControlPoint::GroundControlPoint) continue;

    double minimum_angle = 0;
    vw::ba::triangulate_control_point(cnet[ipt],
                                      opt.camera_models,
                                      minimum_angle);

    //std::cout << "triangulated: " << cnet[ipt].position() << std::endl;

    for ( size_t j = 0; j < cnet[ipt].size(); j++) {
      size_t j_cam_id = cnet[ipt][j].image_id();

      vw::camera::PinholeModel2 * pincam
        = dynamic_cast<vw::camera::PinholeModel2*>(opt.camera_models[j_cam_id].get());
      //      Vector2 pix = cnet[ipt][j].position();
//       std::cout << "id and pix is " << j_cam_id << ' ' << pix << std::endl;
//       std::cout << "---position is " << cnet[ipt].position() << std::endl;
//       Vector2 pix2 = pincam->point_to_pixel(cnet[ipt].position());
//       std::cout << "meas and calc: " << pix << ' ' << pix2 << ' '
//                 << norm_2(pix-pix2) << std::endl;

      Vector3 point = cnet[ipt].position();
      Matrix<double,3,4> M = pincam->camera_matrix();
      double denominator = M(2,0)*point(0) + M(2,1)*point(1) +
        M(2,2)*point(2) + M(2,3);
      //std::cout << "--den " << j_cam_id << ' ' << denominator << std::endl;

    }
  }

#else
  // If gcp were given, initialize all non-gcp control points to the
  // average of the gcp. Otherwise, just use an arbitrary point on the
  // equator.
  int num_gcp = 0;
  Vector3 sum;
  for (int ipt = 0; ipt < (int)cnet.size(); ipt++){
    if (cnet[ipt].type() != ControlPoint::GroundControlPoint) continue;
    sum += cnet[ipt].position();
    num_gcp++;
  }
  Vector3 controlPt;
  if (num_gcp > 0){
    controlPt = sum/num_gcp;
  }else{
    if (opt.datum.name() == UNSPECIFIED_DATUM)
      vw_throw( ArgumentErr() << "No datum was specified.\n");
    controlPt = opt.datum.geodetic_to_cartesian(Vector3(0, 0, 0));
  }
  for (int ipt = 0; ipt < (int)cnet.size(); ipt++){
    if (cnet[ipt].type() == ControlPoint::GroundControlPoint) continue;
    cnet[ipt].set_position(controlPt);
  }

  std::cout << "---temporary!!!" << std::endl;
  controlPt = opt.datum.semi_major_axis() * Vector3(0, 0, -1);

  std::cout << "---control point is " << controlPt << std::endl;

  // Init the pinhole cameras. We put them above ground, close to each
  // other, and pointing down.

  double cameras_to_ground_dist = 100000; // 100 km, hard-coded for now
  cameras_to_ground_dist = 1e+2*1e+3;
  Vector3 highPt = controlPt
    + cameras_to_ground_dist*controlPt/norm_2(controlPt);

  std::cout << "highPt is " << highPt << std::endl;

#if 0
  Vector3 offset(highPt[1], -highPt[0], 0); // perpendicular to highPt
#else
  Vector3 offset(0, -highPt[2], highPt[1]); // perpendicular to highPt
#endif

  offset = (cameras_to_ground_dist/100)*offset/norm_2(offset);
  std::cout << "222---offset is " << offset << std::endl;

  for (int icam = 0; icam < (int)opt.camera_models.size(); icam++){

    // Camera centers
    Vector3 camCtr = highPt + (2*icam-1)*offset;

    // Camera pose
    Vector3 llh = opt.datum.cartesian_to_geodetic(camCtr);
    Matrix3x3 M = opt.datum.lonlat_to_ned_matrix(subvector(llh, 0, 2));
    Quaternion<double> camPose(transpose(M)); // why transpose?

    //std::cout << "---temporary2" << std::endl;
    //M.set_identity();

    std::cout << "matrix is \n" << M << std::endl;

    // Init the models
    vw::camera::PinholeModel2 * pincam
      = dynamic_cast<vw::camera::PinholeModel2*>(opt.camera_models[icam].get());
    VW_ASSERT(pincam != NULL,
              vw::ArgumentErr() << "A pinhole camera expected.\n");

#if 0
    *pincam = vw::camera::PinholeModel2(camCtr,
                                       camPose.rotation_matrix(),
                                       -100, -100,  // focal lengths. Why it turns out negative?
                                       1000, 1000); // pixel offsets
#else
    *pincam = vw::camera::PinholeModel2(camCtr,
                                       M,
                                       2500, 2500,  // focal lengths.
                                       3500, 3500); // pixel offsets

    std::cout << "--cam is " << *pincam << std::endl;

    std::cout << "control point pixel value: " << pincam->point_to_pixel(controlPt)
              << std::endl;
#endif
  }


#if 1
  std::cout << "---test the logic here with gcp!!!!!!!" << std::endl;
  for (int ipt = 0; ipt < (int)cnet.size(); ipt++){
    if (cnet[ipt].type() == ControlPoint::GroundControlPoint) continue;

    double minimum_angle = 0;
    vw::ba::triangulate_control_point(cnet[ipt],
                                      opt.camera_models,
                                      minimum_angle);

    //std::cout << "5triangulated: " << cnet[ipt].position() << std::endl;

    for ( size_t j = 0; j < cnet[ipt].size(); j++) {
      size_t j_cam_id = cnet[ipt][j].image_id();

      vw::camera::PinholeModel2 * pincam
        = dynamic_cast<vw::camera::PinholeModel2*>(opt.camera_models[j_cam_id].get());
//       Vector2 pix = cnet[ipt][j].position();
//       std::cout << "5id and pix is " << j_cam_id << ' ' << pix << std::endl;
//       std::cout << "5---position is " << cnet[ipt].position() << std::endl;
//       Vector2 pix2 = pincam->point_to_pixel(cnet[ipt].position());
//       std::cout << "5meas and calc: " << pix << ' ' << pix2 << ' '
//                 << norm_2(pix-pix2) << std::endl;

      Vector3 point = cnet[ipt].position();
      Matrix<double,3,4> M = pincam->camera_matrix();
      double denominator = M(2,0)*point(0) + M(2,1)*point(1) +
        M(2,2)*point(2) + M(2,3);
      //std::cout << "5--den " << j_cam_id << ' ' << denominator << std::endl;

      if (1 || denominator < 0) {
        cnet[ipt].set_position(controlPt);
        Vector3 point = cnet[ipt].position();
        Matrix<double,3,4> M = pincam->camera_matrix();
        double denominator = M(2,0)*point(0) + M(2,1)*point(1) +
          M(2,2)*point(2) + M(2,3);
        //std::cout << "6--den " << j_cam_id << ' ' << denominator << std::endl;
      }
    }
  }

#endif

#endif

}

// A ceres cost function. Templated by the BundleAdjust model. We pass
// in the observation, the model, and the current camera and point
// indices. The result is the residual, the difference in the
// observation and the projection of the point into the camera,
// normalized by pixel_sigma.
template<class ModelT>
struct BaReprojectionError {
  BaReprojectionError(Vector2 const& observation, Vector2 const& pixel_sigma,
                      ModelT * const ba_model, size_t icam, size_t ipt):
    m_observation(observation),
    m_pixel_sigma(pixel_sigma),
    m_ba_model(ba_model),
    m_icam(icam), m_ipt(ipt){}

  template <typename T>
  bool operator()(const T* const camera, const T* const point,
                  T* residuals) const {

    try{

      size_t num_cameras = m_ba_model->num_cameras();
      size_t num_points  = m_ba_model->num_points();
      VW_ASSERT(m_icam < num_cameras,
                ArgumentErr() << "Out of bounds in the number of cameras");
      VW_ASSERT(m_ipt < num_points,
                ArgumentErr() << "Out of bounds in the number of points");

      // Copy the input data to structures expected by the BA model
      typename ModelT::camera_vector_t cam_vec;
      int cam_len = cam_vec.size();
      for (int c = 0; c < cam_len; c++)
        cam_vec[c] = (double)camera[c];

      typename ModelT::point_vector_t  point_vec;
      for (size_t p = 0; p < point_vec.size(); p++)
        point_vec[p]  = (double)point[p];

      // Project the current point into the current camera
      Vector2 prediction = (*m_ba_model)(m_ipt, m_icam, cam_vec, point_vec);

      // The error is the difference between the predicted and observed position,
      // normalized by sigma.
      residuals[0] = (prediction[0] - m_observation[0])/m_pixel_sigma[0];
      residuals[1] = (prediction[1] - m_observation[1])/m_pixel_sigma[1];

    } catch (const camera::PixelToRayErr& e) {
      // Failed to project into the camera
      residuals[0] = T(1e+20);
      residuals[1] = T(1e+20);
      return false;
    }

    return true;
  }

  // Factory to hide the construction of the CostFunction object from
  // the client code.
  static ceres::CostFunction* Create(Vector2 const& observation,
                                     Vector2 const& pixel_sigma,
                                     ModelT * const ba_model,
                                     size_t icam, // camera index
                                     size_t ipt // point index
                                     ){
    return (new ceres::NumericDiffCostFunction<BaReprojectionError,
            ceres::CENTRAL, 2, ModelT::camera_params_n, ModelT::point_params_n>
            (new BaReprojectionError(observation, pixel_sigma,
                                     ba_model, icam, ipt)));

  }

  Vector2 m_observation;
  Vector2 m_pixel_sigma;
  ModelT * const m_ba_model;
  size_t m_icam, m_ipt;
};

// A ceres cost function. Here we float a pinhole camera's intrinsic
// and extrinsic parameters. The result is the residual, the
// difference in the observation and the projection of the point into
// the camera, normalized by pixel_sigma.
template<class ModelT>
struct BaPinholeError {
  BaPinholeError(Vector2 const& observation, Vector2 const& pixel_sigma,
                      ModelT * const ba_model, size_t icam, size_t ipt):
    m_observation(observation),
    m_pixel_sigma(pixel_sigma),
    m_ba_model(ba_model),
    m_icam(icam), m_ipt(ipt){}

  template <typename T>
  bool operator()(const T* const camera,
                  const T* const point,
                  const T* const intrinsic,
                  T* residuals) const {

    try{

      size_t num_cameras = m_ba_model->num_cameras();
      size_t num_points  = m_ba_model->num_points();
      VW_ASSERT(m_icam < num_cameras,
                ArgumentErr() << "Out of bounds in the number of cameras");
      VW_ASSERT(m_ipt < num_points,
                ArgumentErr() << "Out of bounds in the number of points");

      // Copy the input data to structures expected by the BA model
      typename ModelT::camera_intr_vector_t cam_intr_vec;
      asp::concat_extrinsics_intrinsics<ModelT>(camera, intrinsic, cam_intr_vec);
      typename ModelT::point_vector_t  point_vec;
      for (size_t p = 0; p < point_vec.size(); p++)
        point_vec[p]  = (double)point[p];

      // Project the current point into the current camera
      Vector2 prediction = (*m_ba_model)(m_ipt, m_icam, cam_intr_vec, point_vec);

      // The error is the difference between the predicted and observed position,
      // normalized by sigma.
      residuals[0] = (prediction[0] - m_observation[0])/m_pixel_sigma[0];
      residuals[1] = (prediction[1] - m_observation[1])/m_pixel_sigma[1];

    } catch (const camera::PixelToRayErr& e) {
      // Failed to project into the camera
      residuals[0] = T(1e+20);
      residuals[1] = T(1e+20);
      return false;
    }

    return true;
  }

  // Factory to hide the construction of the CostFunction object from
  // the client code.
  static ceres::CostFunction* Create(Vector2 const& observation,
                                     Vector2 const& pixel_sigma,
                                     ModelT * const ba_model,
                                     size_t icam, // camera index
                                     size_t ipt // point index
                                     ){
    return (new ceres::NumericDiffCostFunction<BaPinholeError,
            ceres::CENTRAL, 2, ModelT::camera_params_n, ModelT::point_params_n,
            ModelT::intrinsic_params_n>
            (new BaPinholeError(observation, pixel_sigma,
                                     ba_model, icam, ipt)));

  }

  Vector2 m_observation;
  Vector2 m_pixel_sigma;
  ModelT * const m_ba_model;
  size_t m_icam, m_ipt;
};

// A ceres cost function. The residual is the difference between the
// observed 3D point and the current (floating) 3D point, normalized by
// xyz_sigma. Used only for ground control points.
struct XYZError {
  XYZError(Vector3 const& observation, Vector3 const& xyz_sigma):
    m_observation(observation), m_xyz_sigma(xyz_sigma){}

  template <typename T>
  bool operator()(const T* const point, T* residuals) const {
    for (size_t p = 0; p < m_observation.size(); p++)
      residuals[p] = ((double)point[p] - m_observation[p])/m_xyz_sigma[p];

    return true;
  }

  // Factory to hide the construction of the CostFunction object from
  // the client code.
  static ceres::CostFunction* Create(Vector3 const& observation,
                                     Vector3 const& xyz_sigma){
    return (new ceres::NumericDiffCostFunction<XYZError, ceres::CENTRAL, 3, 3>
            (new XYZError(observation, xyz_sigma)));

  }

  Vector3 m_observation;
  Vector3 m_xyz_sigma;
};

// A ceres cost function. The residual is the difference between the
// original camera center and the current (floating) camera center.
// This cost function prevents the cameras from straying too far from
// their starting point.
template<class ModelT>
struct CamError {
  typedef typename ModelT::camera_vector_t CamVecT;

  CamError(CamVecT const& orig_cam, double weight):
    m_orig_cam(orig_cam), m_weight(weight){}

  template <typename T>
  bool operator()(const T* const cam_vec, T* residuals) const {

    // Note that we allow the position to vary more than the orientation.
    for (size_t p = 0; p < 3; p++)
      residuals[p] = 1e-6*m_weight*(cam_vec[p] - m_orig_cam[p]);
    for (size_t p = 3; p < m_orig_cam.size(); p++)
      residuals[p] = m_weight*(cam_vec[p] - m_orig_cam[p]);

    return true;
  }

  // Factory to hide the construction of the CostFunction object from
  // the client code.
  static ceres::CostFunction* Create(CamVecT const& orig_cam, double weight){
    return (new ceres::NumericDiffCostFunction<CamError, ceres::CENTRAL,
            ModelT::camera_params_n, ModelT::camera_params_n>
            (new CamError(orig_cam, weight)));

  }

  CamVecT m_orig_cam;
  double m_weight;
};

ceres::LossFunction* get_loss_function(Options const& opt){
  double th = opt.robust_threshold;
  ceres::LossFunction* loss_function;
  if      ( opt.cost_function == "l2"     )
    loss_function = NULL;
  else if ( opt.cost_function == "huber"  )
    loss_function = new ceres::HuberLoss(th);
  else if ( opt.cost_function == "cauchy" )
    loss_function = new ceres::CauchyLoss(th);
  else if ( opt.cost_function == "l1"     )
    loss_function = new ceres::SoftLOneLoss(th);
  else{
    vw_throw( ArgumentErr() << "Unknown cost function: " << opt.cost_function
              << " used with solver: " << opt.ba_type << ".\n" );
  }
  return loss_function;
}

// Add residual block without floating intrinsics
template<class ModelT>
typename boost::disable_if<boost::is_same<ModelT,BAPinholeModel>,void>::type
add_residual_block(ModelT & ba_model,
                   Vector2 const& observation, Vector2 const& pixel_sigma,
                   size_t icam, size_t ipt,
                   double * camera, double * point, double * intrinsics,
                   ceres::LossFunction* loss_function,
                   ceres::Problem & problem){

  ceres::CostFunction* cost_function =
    BaReprojectionError<ModelT>::Create(observation, pixel_sigma,
                                        &ba_model, icam, ipt);
  problem.AddResidualBlock(cost_function, loss_function, camera, point);

}

// Add residual block floating the intrinsics
template<class ModelT>
typename boost::enable_if<boost::is_same<ModelT,BAPinholeModel>,void>::type
add_residual_block(ModelT & ba_model,
                   Vector2 const& observation, Vector2 const& pixel_sigma,
                   size_t icam, size_t ipt,
                   double * camera, double * point, double * intrinsics,
                   ceres::LossFunction* loss_function,
                   ceres::Problem & problem){

  ceres::CostFunction* cost_function =
    BaPinholeError<ModelT>::Create(observation, pixel_sigma,
                                   &ba_model, icam, ipt);
  problem.AddResidualBlock(cost_function, loss_function, camera,
                           point, intrinsics);
  //problem.SetParameterBlockConstant(&intrinsics[0]);
  //std::cout << "---fixing intrinsics!!!" << std::endl;
  //std::cout << "--not fixing intrinsics!!!" << std::endl;
//   std::cout << "--obs pix is " << observation << ' ' << pixel_sigma << std::endl;
//   std::cout << "icam and ipt is " << icam << ' ' << ipt << std::endl;
//   for (int i = 0; i < 6; i++)
//     std::cout << "camera " << camera[i] << " ";
//   std::cout << std::endl;
//   std::cout << "point " << point[0] << ' ' << point[1] << ' ' << point[2] << std::endl;
//   std::cout << "intrinsics " << intrinsics[0] << ' ' << intrinsics[1] << ' ' << intrinsics[2] << std::endl;
}

// Use Ceres to do bundle adjustment. The camera and point variables
// are stored in arrays.  The projection of point into camera is
// accomplished by interfacing with the bundle adjustment model. In
// the future this class can be bypassed.
template <class ModelT>
void do_ba_ceres(ModelT & ba_model, Options& opt ){

  // Don't use this, use opt.cnet!!!!
  ControlNetwork & cnet = *(ba_model.control_network().get());

  size_t num_camera_params    = ModelT::camera_params_n;
  size_t num_point_params     = ModelT::point_params_n;
  size_t num_intrinsic_params = ModelT::intrinsic_params_n;
  size_t num_cameras          = ba_model.num_cameras();
  size_t num_points           = ba_model.num_points();

  // The camera adjustment and point variables concatenated into
  // vectors. The camera adjustments start as 0. The points come from
  // the network.
  std::vector<double> cameras_vec(num_cameras*num_camera_params, 0.0);
  std::vector<double> intrinsics_vec(num_intrinsic_params, 0.0);

  if (!opt.have_input_cams) {
    update_cnet_and_init_cams(opt, *opt.cnet);
    init_optimization_vars(ba_model, opt, cameras_vec,  intrinsics_vec);
  }

  // Camera extrinsics and intrinsics
  double* cameras = &cameras_vec[0];
  double * intrinsics = NULL;
  if (num_intrinsic_params > 0)
    intrinsics = &intrinsics_vec[0];

  // Points
  std::vector<double> points_vec(num_points*num_point_params, 0.0);
  for (size_t ipt = 0; ipt < num_points; ipt++){
    for (size_t q = 0; q < num_point_params; q++){
      points_vec[ipt*num_point_params + q] = cnet[ipt].position()[q];
    }
  }
  double* points = &points_vec[0];

  // The camera positions and orientations before we float them
  std::vector<double> orig_cameras_vec = cameras_vec;

  ceres::Problem problem;

  CameraRelationNetwork<JFeature> crn;
  crn.read_controlnetwork(cnet);

  // Add the cost function component for difference of pixel observations
  typedef CameraNode<JFeature>::iterator crn_iter;
  std::cout << "---fixing focal length!!!" << std::endl;

  std::cout << "---number of points is " << num_points << std::endl;

  for ( size_t icam = 0; icam < crn.size(); icam++ ) {
    for ( crn_iter fiter = crn[icam].begin(); fiter != crn[icam].end(); fiter++ ){

      // The index of the 3D point
      size_t ipt = (**fiter).m_point_id;

      VW_ASSERT(icam < num_cameras,
                ArgumentErr() << "Out of bounds in the number of cameras");
      VW_ASSERT(ipt < num_points,
                ArgumentErr() << "Out of bounds in the number of points");

      // The observed value for the projection of point with index ipt into
      // the camera with index icam.
      Vector2 observation = (**fiter).m_location;
      Vector2 pixel_sigma = (**fiter).m_scale;

      // Each observation corresponds to a pair of a camera and a point
      // which are identified by indices icam and ipt respectively.
      double * camera = cameras + icam * num_camera_params;
      double * point  = points  + ipt * num_point_params;

      ceres::LossFunction* loss_function = get_loss_function(opt);

      add_residual_block(ba_model, observation, pixel_sigma, icam, ipt,
                         camera, point, intrinsics, loss_function, problem);
    }

  }

  // Add ground control points
  for (size_t ipt = 0; ipt < num_points; ipt++){
    if (cnet[ipt].type() != ControlPoint::GroundControlPoint) continue;

    Vector3 observation = cnet[ipt].position();
    Vector3 xyz_sigma = cnet[ipt].sigma();

    ceres::CostFunction* cost_function =
      XYZError::Create(observation, xyz_sigma);

    ceres::LossFunction* loss_function = get_loss_function(opt);

    double * point  = points  + ipt * num_point_params;
    if (opt.have_input_cams)
      problem.AddResidualBlock(cost_function, loss_function, point);

    Vector3 Q;
    for (size_t i = 0; i < num_point_params; i++) {
      std::cout << point[i] << ' ';
      if (i < 3) Q[i] = point[i];
    }
    std::cout << std::endl;
    std::cout << "lon lat height: "
              << opt.datum.cartesian_to_geodetic(Q) << std::endl;
    std::cout << "-- pixel obs is " << observation << std::endl;

    for (ControlPoint::const_iterator it = cnet[ipt].begin(); it != cnet[ipt].end() ; it++) {
      std::cout << "measure position is " << (*it).position() << std::endl;
      std::cout << "measure id is " << (*it).image_id() << std::endl;
    }
  }

  // Add camera constraints
  if (opt.camera_weight > 0){
    for (size_t icam = 0; icam < num_cameras; icam++){

      typename ModelT::camera_vector_t orig_cam;
      for (size_t q = 0; q < num_camera_params; q++)
        orig_cam[q] = orig_cameras_vec[icam * num_camera_params + q];

      ceres::CostFunction* cost_function =
        CamError<ModelT>::Create(orig_cam, opt.camera_weight);

      ceres::LossFunction* loss_function = get_loss_function(opt);

      double * camera  = cameras  + icam * num_camera_params;
      if (opt.have_input_cams)
        problem.AddResidualBlock(cost_function, loss_function, camera);
    }
  }

  // Solve the problem
  ceres::Solver::Options options;
  options.gradient_tolerance = 1e-16;
  options.function_tolerance = 1e-16;
  options.max_num_iterations = opt.max_iterations;
  options.minimizer_progress_to_stdout = (opt.report_level >= vw::ba::ReportFile);

  if (opt.stereo_session_string == "isis")
    options.num_threads = 1;
  else
    options.num_threads = opt.num_threads;

  options.linear_solver_type = ceres::SPARSE_SCHUR;
  //options.ordering_type = ceres::SCHUR;
  //options.eta = 1e-3; // FLAGS_eta;
  //options->max_solver_time_in_seconds = FLAGS_max_solver_time;
  //options->use_nonmonotonic_steps = FLAGS_nonmonotonic_steps;
  //if (FLAGS_line_search) {
  //  options->minimizer_type = ceres::LINE_SEARCH;
  //}

  ceres::Solver::Summary summary;
  ceres::Solve(options, &problem, &summary);
  vw_out() << summary.FullReport() << "\n";
  if (summary.termination_type == ceres::NO_CONVERGENCE){
    // Print a clarifying message, so the user does not think that
    // the algorithm failed.
    vw_out() << "Found a valid solution, but did not reach the actual minimum. It is suggested to increase the number of iterations."
             << std::endl;
  }

  // Copy the latest version of the optimized variables into
  // ba_model. Need this to be able to save the model to disk later.
  typename ModelT::camera_intr_vector_t concat;
  for (size_t icam = 0; icam < num_cameras; icam++){
    asp::concat_extrinsics_intrinsics<ModelT>(&cameras_vec[icam*num_camera_params],
                                              intrinsics,
                                              concat); // output goes here
    ba_model.set_A_parameters(icam, concat);
  }

}

template <class AdjusterT>
void do_ba(typename AdjusterT::model_type & ba_model,
           typename AdjusterT::cost_type const& cost_function,
           Options const& opt) {

  AdjusterT bundle_adjuster(ba_model, cost_function, false, false);

  if ( opt.lambda > 0 )
    bundle_adjuster.set_lambda( opt.lambda );

  std::string iterCameraFile = opt.out_prefix + "iterCameraParam.txt";
  std::string iterPointsFile = opt.out_prefix + "iterPointsParam.txt";

  //Clearing the monitoring text files to be used for saving camera params
  if (opt.save_iteration){
    fs::remove(iterCameraFile);
    fs::remove(iterPointsFile);

    // Write the starting locations
    vw_out() << "Writing: " << iterCameraFile << std::endl;
    vw_out() << "Writing: " << iterPointsFile << std::endl;
    ba_model.bundlevis_cameras_append(iterCameraFile);
    ba_model.bundlevis_points_append(iterPointsFile);
  }

  BundleAdjustReport<AdjusterT>
    reporter( "Bundle Adjust", ba_model, bundle_adjuster,
              opt.report_level );

  double abs_tol = 1e10, rel_tol=1e10;
  double overall_delta = 2;
  int no_improvement_count = 0;
  while ( true ) {
    // Determine if it is time to quit
    if ( bundle_adjuster.iterations() >= opt.max_iterations ) {
      reporter() << "Triggered 'Max Iterations'\n";
      break;
    } else if ( abs_tol < 0.01 ) {
      reporter() << "Triggered 'Abs Tol " << abs_tol << " < 0.01'\n";
      break;
    } else if ( rel_tol < 1e-6 ) {
      reporter() << "Triggered 'Rel Tol " << rel_tol << " < 1e-10'\n";
      break;
    } else if ( no_improvement_count > 4 ) {
      reporter() << "Triggered break, unable to improve after "
                 << no_improvement_count << " iterations\n";
      break;
    }

    overall_delta = bundle_adjuster.update( abs_tol, rel_tol );
    reporter.loop_tie_in();

    // Writing Current Camera Parameters to file for later reading
    if (opt.save_iteration) {
      ba_model.bundlevis_cameras_append(iterCameraFile);
      ba_model.bundlevis_points_append(iterPointsFile);
    }
    if ( overall_delta == 0 )
      no_improvement_count++;
    else
      no_improvement_count = 0;
  }
  reporter.end_tie_in();

}

void save_cnet_as_csv(Options& opt, std::string const& cnetFile){

  // Save the input control network in the csv file format used by ground
  // control points.

  if (opt.datum.name() == UNSPECIFIED_DATUM)
    vw_throw( ArgumentErr() << "No datum was specified. "
              << "Cannot save control network as csv.\n" );

  vw_out() << "Writing: " << cnetFile << std::endl;
  std::ofstream ofs(cnetFile.c_str());
  ofs.precision(17);

  int count = 0;
  ControlNetwork & cnet = *opt.cnet.get();
  for ( ControlNetwork::const_iterator iter = cnet.begin();
        iter != cnet.end(); ++iter ) {

    // If to dump only gcp
    //if ( (*iter).type() != ControlPoint::GroundControlPoint ) continue;

    count++;

    // lon,lat,height
    Vector3 llr = opt.datum.cartesian_to_geodetic((*iter).position());

    // convert to lat,lon,height
    std::swap(llr[0], llr[1]);

    Vector3 sigma = (*iter).sigma();

    ofs << count  << ' ' << llr[0] << ' ' << llr[1] << ' ' << llr[2] << ' ';
    ofs << sigma[0]  << ' ' << sigma[1] << ' ' << sigma[2] << ' ';

    for ( ControlPoint::const_iterator measure = (*iter).begin();
          measure != (*iter).end(); ++measure ) {

      ofs << opt.image_files[measure->image_id()] << ' '
          << measure->position()[0] << ' ' << measure->position()[1] << ' '
          << measure->sigma()[0]    << ' ' << measure->sigma()[1];

      if ( measure+1 != (*iter).end())
        ofs << ' ';
      else
        ofs << std::endl;
    }
  }

  return;
}

// Use given cost function. Switch based on solver.
template<class ModelType, class CostFunType>
void do_ba_costfun(CostFunType const& cost_fun, Options& opt){

  ModelType ba_model(opt.camera_models, opt.cnet);

  if ( opt.ba_type == "ceres" ) {
    do_ba_ceres<ModelType>(ba_model, opt);
  } else if ( opt.ba_type == "robustsparse" ) {
    do_ba<AdjustRobustSparse< ModelType,CostFunType> >(ba_model, cost_fun, opt);
  } else if ( opt.ba_type == "robustref" ) {
    do_ba<AdjustRobustRef< ModelType,CostFunType> >(ba_model, cost_fun, opt);
  } else if ( opt.ba_type == "sparse" ) {
    do_ba<AdjustSparse< ModelType, CostFunType > >(ba_model, cost_fun, opt);
  }else if ( opt.ba_type == "ref" ) {
    do_ba<AdjustRef< ModelType, CostFunType > >(ba_model, cost_fun, opt);
  }

  // Save the models to disk
  for (size_t icam = 0; icam < ba_model.num_cameras(); icam++){
    std::string adjust_file = asp::bundle_adjust_file_name(opt.out_prefix,
                                                           opt.image_files[icam]);
    vw_out() << "Writing: " << adjust_file << std::endl;
    ba_model.write_adjustment(icam, adjust_file);
  }

}

// Do BA with given model. Switch based on cost function.
template<class ModelType>
void do_ba_with_model(Options& opt){

  if ( opt.cost_function == "cauchy" ) {
    do_ba_costfun<ModelType, CauchyError>
      (CauchyError(opt.robust_threshold), opt);
  }else if ( opt.cost_function == "pseudohuber" ) {
    do_ba_costfun<ModelType, PseudoHuberError>
      (PseudoHuberError(opt.robust_threshold), opt );
  } else if ( opt.cost_function == "huber" ) {
    do_ba_costfun<ModelType, HuberError>
      (HuberError(opt.robust_threshold), opt );
  } else if ( opt.cost_function == "l1" ) {
    do_ba_costfun<ModelType, L1Error>( L1Error(), opt );
  } else if ( opt.cost_function == "l2" ) {
    do_ba_costfun<ModelType, L2Error>( L2Error(), opt );
  }else{
    vw_throw( ArgumentErr() << "Unknown cost function: " << opt.cost_function
              << ". Options are: Cauchy, PseudoHuber, Huber, L1, L2.\n" );
  }

}

void handle_arguments( int argc, char *argv[], Options& opt ) {
  po::options_description general_options("");
  general_options.add_options()
//     ("cnet,c", po::value(&opt.cnet_file),
//      "Load a control network from a file (optional).")
    ("output-prefix,o", po::value(&opt.out_prefix), "Prefix for output filenames.")
    ("bundle-adjuster", po::value(&opt.ba_type)->default_value("Ceres"),
     "Choose a solver from: Ceres, RobustSparse, RobustRef, Sparse, Ref.")
    ("cost-function", po::value(&opt.cost_function)->default_value("Cauchy"),
     "Choose a cost function from: Cauchy, PseudoHuber, Huber, L1, L2.")
    ("robust-threshold", po::value(&opt.robust_threshold)->default_value(0.5),
     "Set the threshold for robust cost functions.")
    ("datum", po::value(&opt.datum_str)->default_value(""), "Use this datum (needed only if ground control points are used). [WGS_1984, D_MOON (radius is assumed to be 1,737,400 meters), D_MARS (radius is assumed to be 3,396,190 meters), etc.]")
    ("semi-major-axis", po::value(&opt.semi_major)->default_value(0),
     "Explicitly set the datum semi-major axis in meters (needed only if ground control points are used).")
    ("semi-minor-axis", po::value(&opt.semi_minor)->default_value(0),
     "Explicitly set the datum semi-minor axis in meters (needed only if ground control points are used).")
    ("session-type,t", po::value(&opt.stereo_session_string)->default_value(""),
     "Select the stereo session type to use for processing. Options: pinhole isis dg rpc. Usually the program can select this automatically by the file extension.")
    ("min-matches", po::value(&opt.min_matches)->default_value(30),
     "Set the minimum  number of matches between images that will be considered.")
    ("max-iterations", po::value(&opt.max_iterations)->default_value(1000),
     "Set the maximum number of iterations.")
    ("overlap-limit", po::value(&opt.overlap_limit)->default_value(3),
     "Limit the number of subsequent images to search for matches to the current image to this value.")
    ("camera-weight", po::value(&opt.camera_weight)->default_value(1.0),
     "The weight to give to the constraint that the camera positions/orientations stay close to the original values (only for the Ceres solver).")
    ("lambda,l", po::value(&opt.lambda)->default_value(-1),
     "Set the initial value of the LM parameter lambda (ignored for the Ceres solver).")
    ("report-level,r",po::value(&opt.report_level)->default_value(10),
     "Use a value >= 20 to get increasingly more verbose output.");
//     ("save-iteration-data,s", "Saves all camera information between iterations to output-prefix-iterCameraParam.txt, it also saves point locations for all iterations in output-prefix-iterPointsParam.txt.");
  general_options.add( asp::BaseOptionsDescription(opt) );

  po::options_description positional("");
  positional.add_options()
    ("input-files", po::value(&opt.image_files));

  po::positional_options_description positional_desc;
  positional_desc.add("input-files", -1);

  std::string usage("<images> <cameras> <optional ground control points> -o <output prefix> [options]");
  bool allow_unregistered = false;
  std::vector<std::string> unregistered;
  po::variables_map vm =
    asp::check_command_line(argc, argv, opt, general_options, general_options,
                            positional, positional_desc, usage,
                             allow_unregistered, unregistered);

  if ( opt.image_files.empty() )
    vw_throw( ArgumentErr() << "Missing input image files.\n"
              << usage << general_options );
  opt.gcp_files = asp::extract_gcps( opt.image_files );
  opt.camera_files = asp::extract_cameras( opt.image_files );

  if ( opt.overlap_limit <= 0 )
    vw_throw( ArgumentErr() << "Must allow search for matches between "
              << "at least each image and its subsequent one.\n"
              << usage << general_options );

  if ( opt.camera_weight < 0.0 )
    vw_throw( ArgumentErr() << "The camera weight must be non-negative.\n"
              << usage << general_options );

  // See if we start with initial cameras or no cameras
  opt.have_input_cams
    = (!opt.camera_files.empty()) || asp::images_are_cubes(opt.image_files);

  if (!opt.gcp_files.empty() || !opt.have_input_cams){
    // Need to read the datum if we have gcps.
    if (opt.datum_str != ""){
      // If the user set the datum, use it.
      opt.datum.set_well_known_datum(opt.datum_str);
    }else if (opt.semi_major > 0 && opt.semi_minor > 0){
      // Otherwise, if the user set the semi-axes, use that.
      opt.datum = cartography::Datum("User Specified Datum",
                                     "User Specified Spheroid",
                                     "Reference Meridian",
                                     opt.semi_major, opt.semi_minor, 0.0);
    }else{
      if (!opt.gcp_files.empty())
        vw_throw( ArgumentErr() << "When ground control points are used, "
                  << "the datum must be specified.\n"
                  << usage << general_options );
      else if (!opt.have_input_cams)
        vw_throw( ArgumentErr() << "When there is no input camera information, "
                  << "the datum must be specified.\n"
                  << usage << general_options );
    }
    vw_out() << "Will use datum: " << opt.datum << std::endl;

  }

  if ( opt.out_prefix.empty() )
    vw_throw( ArgumentErr() << "Missing output prefix.\n"
              << usage << general_options  );

  // Create the output directory
  asp::create_out_dir(opt.out_prefix);

  // Turn on logging to file
  asp::log_to_file(argc, argv, "", opt.out_prefix);

  opt.save_iteration = vm.count("save-iteration-data");
  boost::to_lower( opt.stereo_session_string );
  boost::to_lower( opt.ba_type );
  boost::to_lower( opt.cost_function );
  if ( !( opt.ba_type == "ceres"        ||
          opt.ba_type == "robustsparse" ||
          opt.ba_type == "robustref"    ||
          opt.ba_type == "sparse"       ||
          opt.ba_type == "ref"
          ) )
    vw_throw( ArgumentErr() << "Unknown bundle adjustment version: " << opt.ba_type
              << ". Options are: [Ceres, RobustSparse, RobustRef, Sparse, Ref]\n" );
}

int main(int argc, char* argv[]) {

  Options opt;
  try {
    handle_arguments( argc, argv, opt );

    int num_images = opt.image_files.size();
    // Create the stereo session. Try to auto-guess the session type.
    if (num_images <= 1)
      vw_throw( ArgumentErr() << "Must have at least two image "
                << "files to do bundle adjustment.\n" );

    // If there are no camera files, then the image files have the
    // camera information.
    if (opt.camera_files.empty()){
      for (int i = 0; i < num_images; i++)
        opt.camera_files.push_back("");
    }

    // Sanity check
    if (num_images != (int)opt.camera_files.size())
      vw_throw(ArgumentErr() << "Must have as many cameras as we have images.\n");

    if (!opt.have_input_cams){
      if (opt.stereo_session_string != "" && opt.stereo_session_string != "pinhole"){
        vw_throw( ArgumentErr()
                  << "No input cameras were provided. "
                  << "The only supported stereo session is Pinhole.\n");
      }
      opt.stereo_session_string = "pinhole";
    }

    // Create the stereo session. This will attempt to identify the
    // session type.
    typedef boost::scoped_ptr<asp::StereoSession> SessionPtr;
    SessionPtr session
      (asp::StereoSessionFactory::create(opt.stereo_session_string, opt,
                                         opt.image_files[0], opt.image_files[1],
                                         opt.camera_files[0], opt.camera_files[1],
                                         opt.out_prefix));

    // Read in the camera model and image info for the input images.
    for (int i = 0; i < num_images; i++){
      vw_out(DebugMessage,"asp") << "Loading: "
                                 << opt.image_files[i] << ' '
                                 << opt.camera_files[i] << "\n";
      if (opt.have_input_cams){
        opt.camera_models.push_back(session->camera_model(opt.image_files[i],
                                                          opt.camera_files[i]));
      }
      else{
        // Create new cameras from scratch.
        opt.camera_models.push_back(boost::shared_ptr<vw::camera::CameraModel>
                                    (new PinholeModel2()));
      }
    }

    // Create the match points
    for (int i = 0; i < num_images; i++){
      for (int j = i+1; j <= std::min(num_images-1, i+opt.overlap_limit); j++){
        std::string image1 = opt.image_files[i];
        std::string image2 = opt.image_files[j];
        std::string match_filename
          = ip::match_filename(opt.out_prefix, image1, image2);
        boost::shared_ptr<DiskImageResource>
          rsrc1( DiskImageResource::open(image1) ),
          rsrc2( DiskImageResource::open(image2) );
        float nodata1, nodata2;
        session->get_nodata_values(rsrc1, rsrc2, nodata1, nodata2);
        try{
          // IP matching may not succeed for all pairs
          session->ip_matching(image1, image2, nodata1, nodata2, match_filename,
                               opt.camera_models[i].get(),
                               opt.camera_models[j].get());
        } catch ( const std::exception& e ){
          vw_out(WarningMessage) << e.what() << std::endl;
        }
      }
    }

    opt.cnet.reset( new ControlNetwork("BundleAdjust") );
    if ( opt.cnet_file.empty() ) {
      build_control_network( opt.have_input_cams,
                             (*opt.cnet), opt.camera_models,
                             opt.image_files, opt.min_matches,
                             opt.out_prefix);

      if (opt.have_input_cams)
        add_ground_control_points( (*opt.cnet), opt.image_files,
                                   opt.gcp_files.begin(), opt.gcp_files.end(),
                                   opt.datum);

      //opt.cnet->write_binary(opt.out_prefix + "-control");
      //save_cnet_as_csv(opt, opt.out_prefix + "-cnet.csv");

    } else  {
      vw_out() << "Loading control network from file: "
               << opt.cnet_file << "\n";

      // Deciding which Control Network we have
      std::vector<std::string> tokens;
      boost::split( tokens, opt.cnet_file, boost::is_any_of(".") );
      if ( tokens.back() == "net" ) {
        // An ISIS style control network
        opt.cnet->read_isis( opt.cnet_file );
      } else if ( tokens.back() == "cnet" ) {
        // A VW binary style
        opt.cnet->read_binary( opt.cnet_file );
      } else {
        vw_throw( IOErr() << "Unknown Control Network file extension, \""
                  << tokens.back() << "\"." );
      }
    }

    if (opt.have_input_cams) {
      do_ba_with_model<BundleAdjustmentModel>(opt);
      return 0;
    }


    BAPinholeModel ba_model(opt.camera_models, opt.cnet);
    ControlNetwork & cnet = *opt.cnet.get(); // alias

    do_ba_ceres<BAPinholeModel>(ba_model, opt);

    // Now is the right time to add the gcp
    std::cout << "---cnet size " << (int)cnet.size() << std::endl;
    if (!opt.have_input_cams)
      add_ground_control_points( (*opt.cnet), opt.image_files,
                                 opt.gcp_files.begin(), opt.gcp_files.end(),
                                 opt.datum);

    std::cout << "--cnet size " << (int)cnet.size() << std::endl;

    ba_model.get_pinhole_models(opt.camera_models);
#if 1
    // Use GCP to transform the cameras to the correct coordinate system

    for (size_t icam = 0; icam < opt.camera_models.size(); icam++){
      vw::camera::PinholeModel2 * pincam
        = dynamic_cast<vw::camera::PinholeModel2*>(opt.camera_models[icam].get());
      VW_ASSERT(pincam != NULL,
                vw::ArgumentErr() << "A pinhole camera expected.\n");

      std::cout << "---before rotation, camera is " << *pincam << std::endl;
    }

    int num_gcp = 0;
    for (int ipt = 0; ipt < (int)cnet.size(); ipt++){
      if (cnet[ipt].type() != ControlPoint::GroundControlPoint) continue;
      num_gcp++;
    }
    Eigen::Matrix3Xd in(3, num_gcp), out(3, num_gcp);
    int num_good_gcp = 0;
    for (int ipt = 0; ipt < (int)cnet.size(); ipt++){
      if (cnet[ipt].type() != ControlPoint::GroundControlPoint) continue;

      ControlPoint cp_new = cnet[ipt];
      // Making minimum_angle below big may throw away valid points at this stage
      // really???
      double minimum_angle = 0;
      vw::ba::triangulate_control_point(cp_new,
                                        opt.camera_models,
                                        minimum_angle);
      Vector3 inp = cp_new.position();
      Vector3 outp  = cnet[ipt].position();
      std::cout << "---triangulated: " << cnet[ipt].position() << ' '
                << cp_new.position() << std::endl;
      if (inp == Vector3() || outp == Vector3())
        continue;

#if 1
      for ( size_t j = 0; j < cnet[ipt].size(); j++) {
        size_t j_cam_id = cnet[ipt][j].image_id();

        vw::camera::PinholeModel2 * pincam
          = dynamic_cast<vw::camera::PinholeModel2*>(opt.camera_models[j_cam_id].get());
        Vector2 pix = cnet[ipt][j].position();
        Vector3 point = cp_new.position();
        //           std::cout << "3id and pix is " << j_cam_id << ' ' << pix << std::endl;
        //           std::cout << "3---position is " << point << std::endl;
        //           Vector2 pix2 = pincam->point_to_pixel(point);
        //           std::cout << "3meas and calc: " << pix << ' ' << pix2 << ' '
        //                     << norm_2(pix-pix2) << std::endl;

        Matrix<double,3,4> M = pincam->camera_matrix();
        double denominator = M(2,0)*point(0) + M(2,1)*point(1) +
          M(2,2)*point(2) + M(2,3);
        //std::cout << "3--den " << j_cam_id << ' ' << denominator << std::endl;

      }
#endif

      in.col(num_good_gcp)  << inp[0], inp[1], inp[2];
      out.col(num_good_gcp) << outp[0], outp[1], outp[2];
      std::cout << "--in is " << in.col(num_good_gcp).transpose() << std::endl;
      std::cout << "--ou is " << out.col(num_good_gcp).transpose() << std::endl;
      num_good_gcp++;
    }

    in.conservativeResize(Eigen::NoChange_t(), num_good_gcp);
    out.conservativeResize(Eigen::NoChange_t(), num_good_gcp);
    std::cout << "---must check here that we have at least 3 points!!!"
              << std::endl;

    Eigen::Affine3d T = Find3DAffineTransform(in, out);
    double scale = pow(T.linear().determinant(), 1.0 / 3.0);

    std::cout << "--det is " << T.linear().determinant() << std::endl;
    std::cout << "Transform to world coordinates." << std::endl;
    std::cout << "Rotation:\n" << T.linear() / scale << std::endl;
    std::cout << "Scale:\n" << scale << std::endl;
    std::cout << "Translation:\n" << T.translation().transpose()
              << std::endl;

    for (int i = 0; i < in.cols(); i++) {
      Eigen::Vector3d ip = T*in.col(i);
      Eigen::Vector3d op = out.col(i);

      //std::cout << "vals and err: " << ip.transpose() << ' ' << op.transpose()
      //          << ' ' << (ip-op).norm() << std::endl;
    }

    // Represent the transform as y = scale*A*x + b
    Eigen::MatrixXd eA = T.linear()/scale;
    Eigen::Vector3d eb = T.translation();

    // Store in VW matrices and vectors
    vw::Matrix3x3 A;
    for (size_t r = 0; r < A.rows(); r++) {
      for (size_t c = 0; c < A.cols(); c++) {
        A(r, c) = eA(r, c);
      }
    }
    vw::Vector3 b;
    for (size_t i = 0; i < b.size(); i++)
      b[i] = eb[i];


    // Apply the transform to the cameras
    std::cout << "---do the xyz points" << std::endl;
    int num_cams = opt.camera_models.size();
    for (int icam = 0; icam < num_cams; icam++){
      vw::camera::PinholeModel2 * pincam
        = dynamic_cast<vw::camera::PinholeModel2*>(opt.camera_models[icam].get());
      VW_ASSERT(pincam != NULL,
                vw::ArgumentErr() << "A pinhole camera expected.\n");

      vw::Vector3 position = pincam->camera_center();
      vw::Quat pose = pincam->camera_pose();

      vw::Vector2 fl = pincam->focal_length();
      vw::Vector2 po = pincam->point_offset();

      // New position and rotation
      position = scale*A*position + b;
      vw::Quat aq(A);
      pose = aq*pose;

      *pincam = vw::camera::PinholeModel2(position,
                                          pose.rotation_matrix(),
                                          fl[0], fl[1],  // focal lengths
                                          po[0], po[1]); // pixel offsets

      std::cout << "model: " << *pincam << std::endl;
    }


    ba_model.put_pinhole_models(opt.camera_models);
    ba_model.get_pinhole_models(opt.camera_models);
    ba_model.put_pinhole_models(opt.camera_models);


#if 1
    for (int ipt = 0; ipt < (int)cnet.size(); ipt++){

      std::string v;
      if (cnet[ipt].type() == ControlPoint::GroundControlPoint)
        v = "++++control ";
      else
        v = "----regular ";

      if (cnet[ipt].type() != ControlPoint::GroundControlPoint) {
        double minimum_angle = 0;
        vw::ba::triangulate_control_point(cnet[ipt],
                                          opt.camera_models,
                                          minimum_angle);
      }

      //std::cout << v << "triangulated: " << cnet[ipt].position() << std::endl;

      for ( size_t j = 0; j < cnet[ipt].size(); j++) {
        size_t j_cam_id = cnet[ipt][j].image_id();

        vw::camera::PinholeModel2 * pincam
          = dynamic_cast<vw::camera::PinholeModel2*>(opt.camera_models[j_cam_id].get());
        //Vector2 pix = cnet[ipt][j].position();
        //           std::cout << v << "id and pix is " << j_cam_id << ' ' << pix << std::endl;
        //           std::cout << v << "---position is " << cnet[ipt].position() << std::endl;
        //           Vector2 pix2 = pincam->point_to_pixel(cnet[ipt].position());
        //           std::cout << v << "meas and calc: " << pix << ' ' << pix2 << ' '
        //                     << norm_2(pix-pix2) << std::endl;

        Vector3 point = cnet[ipt].position();
        Matrix<double,3,4> M = pincam->camera_matrix();
        double denominator = M(2,0)*point(0) + M(2,1)*point(1) +
          M(2,2)*point(2) + M(2,3);
        //std::cout << v << "--den " << j_cam_id << ' ' << denominator << std::endl;

      }
    }
#endif



    for (int ipt = 0; ipt < (int)cnet.size(); ipt++){
      if (cnet[ipt].type() != ControlPoint::GroundControlPoint) continue;

      ControlPoint cp_new = cnet[ipt];
      double minimum_angle = 0.0001;
      vw::ba::triangulate_control_point(cp_new,
                                        opt.camera_models,
                                        minimum_angle);
      Vector3 inp  = cnet[ipt].position();
      Vector3 outp = cp_new.position();
      std::cout << "---triangulated: " << cnet[ipt].position() << ' '
                << cp_new.position()  << ' '
                << cnet[ipt].position() - cp_new.position() << ' '
                << norm_2(cnet[ipt].position() - cp_new.position())<< std::endl;

      std::cout << "lon lat height diff: " << opt.datum.cartesian_to_geodetic(cnet[ipt].position()) -  opt.datum.cartesian_to_geodetic( cp_new.position())  << std::endl;

    }

#endif

#if 0
    int num_cams = opt.camera_models.size();
    for (int icam = 0; icam < num_cams; icam++){
      vw::camera::PinholeModel2 * pincam
        = dynamic_cast<vw::camera::PinholeModel2*>(opt.camera_models[icam].get());
      VW_ASSERT(pincam != NULL,
                vw::ArgumentErr() << "A pinhole camera expected.\n");
      std::cout << "model: " << *pincam << std::endl;
    }
#endif


    // Save the camera models to disk
    std::vector<std::string> cam_files;
    for (int icam = 0; icam < (int)opt.camera_models.size(); icam++){
      std::string cam_file
        = asp::bundle_adjust_file_name(opt.out_prefix, opt.image_files[icam]);
      cam_file = fs::path(cam_file).replace_extension("pinhole").string();
      cam_files.push_back(cam_file);
    }
    ba_model.write_camera_models(cam_files);

  } ASP_STANDARD_CATCHES;
}
