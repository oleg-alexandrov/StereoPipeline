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

/// \file JitterSolveRigUtils.h

// Functions invoked in jitter_solve.cc that need a rig.

#ifndef __ASP_CAMERA_JITTER_SOLVE_RIG_UTILS_H__
#define __ASP_CAMERA_JITTER_SOLVE_RIG_UTILS_H__

#include <asp/Camera/CsmModel.h>

#include <string>
#include <map>
#include <vector>

class UsgsAstroLsSensorModel;
class UsgsAstroFrameSensorModel;

namespace rig {
  class RigSet;
}

namespace asp {

typedef std::map<int, std::map<double, int>> TimestampMap;

// Create enum for senor type, which can be frame or linescan
enum RigSensorType {RIG_LINESCAN_SENSOR, RIG_FRAME_SENSOR}; 

// For each camera this will info info that will tie it to the rig
struct RigCamInfo {
  int sensor_id; // The sensor id in the rig set
  RigSensorType sensor_type; // frame or linescan
  
  // The time at the image center for linescan and mid group for frame. Here the
  // group is formed of all camera images acquired within a contigious time with
  // the same frame sensor.
  double mid_group_time;
  
  // For linescan: the starting and ending time for positions/orientations.
  // For frame: both are the current pose time.
  double beg_pose_time, end_pose_time;
  
  // The index of the camera in opt.camera_models
  int cam_index;
  // The index of the reference camera in opt.camera_models. The reference
  // camera is assumed to be linescan.
  int ref_cam_index; 
  
  std::string image_file, camera_file;
  
  // The group name for the all cameras acquired with all sensors on a rig
  // within a contigious time interval. Acquisitions at at a different time
  // intervals and/or a different rig must have different group names.
  std::string rig_group;
   
  // Declare the constructor
  RigCamInfo();
};

// Find the time bounds for a linescan sensor
void linescanTimeBounds(UsgsAstroLsSensorModel const* ls_model, 
                        // Outputs
                        double & beg_time, double & end_time);

// Book-keeping needed to tie each camera to the rig
void populateRigCamInfo(rig::RigSet const& rig,
                        std::vector<std::string> const& image_files,
                        std::vector<std::string> const& camera_files,
                        std::vector<asp::CsmModel*> const& csm_models,
                        std::map<int, int> const& cam2group,
                        bool use_initial_rig_transforms,
                        // Outputs
                        std::vector<RigCamInfo> & rig_cam_info,
                        std::vector<double>     & ref_to_curr_sensor_vec,
                        TimestampMap & timestamp_map);

// Given a reference linescan camera and the transform from it to the current
// camera, find the current camera to world transform as an array.
void linescanToCurrSensorTrans( UsgsAstroLsSensorModel const & ref_ls_cam,
                               double curr_time,
                               double const* ref_to_curr_trans,
                               // Output
                               double * cam2world_arr);

// Given a frame linescan camera and the transform from it to the current
// camera, find the current camera to world transform as an array.
void frameToCurrSensorTrans(std::vector<double>       const& frame_params,
                            asp::RigCamInfo           const& rig_cam_info,
                            std::map<int, int>        const& cam2group,
                            TimestampMap              const& timestamp_map,
                            double                    const* ref_to_curr_trans,
                            // Output
                            double * cam2world_arr);

// Given a reference linescan camera and the transform from it to the current
// linescan camera, update the the current camera poses within the given range.
void updateLinescanWithRig( UsgsAstroLsSensorModel const & ref_ls_cam,
                           double const* ref_to_curr_trans,
                           UsgsAstroLsSensorModel & curr_ls_cam, // update this
                           // Range of quat and position indices to update.
                           // The default is to update all.
                           int beg_quat_index = -1, int end_quat_index = -1,
                           int beg_pos_index = -1, int end_pos_index = -1);

// Update the rig with the optimized transforms
void updateRig(std::vector<double> const& ref_to_curr_sensor_vec,
               rig::RigSet & rig);

// Find the times and indices bracketing a given time
bool timestampBrackets(double time, 
                  std::map<double, int> const& timestamps,
                  // Outputs
                  double & time1, double & time2,
                  int & index1, int & index2);

// Find the timestamps bracketing the camera with the given rig_cam_info.
// This is a wrapper around the above function.
bool timestampBrackets(asp::RigCamInfo     const & rig_cam_info,
                        std::map<int, int> const & cam2group,
                        TimestampMap       const & timestamp_map,
                        // Outputs
                        double & beg_ref_time, double & end_ref_time, 
                        int & beg_ref_index, int & end_ref_index);

// Given two poses, each as x, y, z, qx, qy, qz, qw, linearly interpolate between them.
void interpPose(double time1, double time2, 
                double const* pose1, double const* pose2, double time, 
                // Output
                double* pose);

// Given reference poses at time1 and time2, a time in between, and the transform
// from the rig from the reference to the current sensor, find the current sensor
// pose (camera-to-world) at the given time.
void interpCurrPose(double time1, double time2, double time, 
                    double const* pose1, double const* pose2, 
                    double const* ref_to_curr_trans,
                    // Output
                    double * cam2world_arr);

} // end namespace asp

#endif //__ASP_CAMERA_JITTER_SOLVE_RIG_UTILS_H__
