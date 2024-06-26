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

#ifndef ASP_GOTCHA_CDENSIFY_H
#define ASP_GOTCHA_CDENSIFY_H

#include <asp/Gotcha/CProcBlock.h>
#include <asp/Gotcha/CDensifyParam.h>

#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/calib3d/calib3d.hpp>

#include <fstream>
#include <iostream>
#include <ostream>
#include <stdio.h>
#include <stdlib.h>
#include <termios.h>

namespace gotcha {

class CDensify: public CProcBlock {
public:
  CDensify();
  CDensify(CDensifyParam paramDense, std::vector<CTiePt> const& vecTPs,
           cv::Mat imgL, cv::Mat imgR,
           cv::Mat input_dispX, cv::Mat input_dispY, cv::Mat Mask);

  void setParameters(CDensifyParam paramDense, std::vector<CTiePt> const& vecTPs,
                     cv::Mat imgL, cv::Mat imgR,
                     cv::Mat input_dispX, cv::Mat input_dispY, cv::Mat Mask);
  
    int performDensitification(cv::Mat & output_dispX, cv::Mat & cv_output_dispY);
    int getNumTotTps() const {return m_vectpAdded.size();}
    static bool compareTP(CTiePt tpX, CTiePt tpY){
        return (tpX.m_fSimVal < tpY.m_fSimVal);
    }

    std::string getProcType(){
        return "GOTCHA_FROM_DISP";
    }
private:
    //void loadImages();
    std::vector<CTiePt> getIntToFloatSeed(std::vector<CTiePt>& vecTPSrc); // get integer Seed point pairs from a float seed point pair
    bool doGotcha(const cv::Mat& matImgL, const cv::Mat& matImgR, std::vector<CTiePt>& vectpSeeds,
                  const CGOTCHAParam& paramGotcha, std::vector<CTiePt>& vectpAdded);
    bool doTileGotcha(const cv::Mat& matImgL, const cv::Mat& matImgR, const std::vector<CTiePt>& vectpSeeds,
                      const CGOTCHAParam& paramGotcha, std::vector<CTiePt>& mvectpAdded,
                      const cv::Rect_<float> rectTileL, cv::Mat& matSimMap, std::vector<bool>& pLUT); //IMARS
    void removePtInLUT(std::vector<CTiePt>& vecNeiTp, const std::vector<bool>& pLUT, const int nWidth); //IMARS
    void removeOutsideImage(std::vector<CTiePt>& vecNeiTp, const cv::Rect_<float> rectTileL, const cv::Rect_<float> rectImgR);
    void getNeighbour(const CTiePt tp, std::vector<CTiePt>& vecNeiTp, const int nNeiType, const cv::Mat& matSim);
    void getDisffusedNei(std::vector<CTiePt>& vecNeiTp, const CTiePt tp, const cv::Mat& matSim);
    void breakIntoSubRect(cv::Rect_<float> rectParent, std::vector< cv::Rect_<float> >& vecRes);
    void makeTiles(std::vector< cv::Rect_<float> >& vecRectTiles, int nMin);
    bool isHavingTP(std::vector<CTiePt>& vecNeiTp, CTiePt tp); // used in diffused neighbour
    bool doPGotcha(int nNeiType);
    int getTotPyramidLev(int nszPatch);
    void makeDataProducts();
    bool saveResult(cv::Mat & output_dispX, cv::Mat & cv_output_dispY);
    bool saveProjLog (std::string strFile);
    bool saveLog();
    bool saveResLog(std::string strFile);
    bool loadTPForDensification(std::string strTPFile);

private:
  // inputs
  CDensifyParam m_paramDense;
  cv::Mat m_imgL, m_imgR;
  cv::Mat m_input_dispX, m_input_dispY;
  cv::Mat m_Mask;
  std::vector<CTiePt> m_vecTPs;
  
  // outputs
  std::vector<CTiePt> m_vectpAdded; // final tp result (i.e. seed tps + newly added tps)
  // Disparity maps:
  // the disparity of a point p_i at x,y in the left image is stored at
  // (x, y) in the disparity maps, e.g., (m_matDisMapX and m_matDisMapY).
  // For example, the position of the corresponding point p'_i in the right image can be defined as
  // (x+x_dist, y+y_dist), where 'xdist' and 'ydist' are value at (x,y) in m_matDisMapX and m_matDisMapY, respectively.
  // When disparity values are not known, the vaule of the disparity map is set to zero.
  cv::Mat m_matDisMapX;
  cv::Mat m_matDisMapY;
  cv::Mat m_matDisMapSim; // matching score map: smaller score is better score and -1 indicate unknown
  
  cv::Mat m_matHL;        // used when the rectified images are used
  cv::Mat m_matHR;
  int nNumSeedTPs;
  double procTime;
};

} // end namespace gotcha

#endif // ASP_GOTCHA_CDENSIFY_H
