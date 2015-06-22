

// Copyright (c) 2012, 2013, 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/sfm.hpp"
#include "openMVG/image/image.hpp"

using namespace openMVG;
using namespace openMVG::cameras;
using namespace openMVG::geometry;
using namespace openMVG::sfm;

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <fstream>

int main(int argc, char **argv)
{
  CmdLine cmd;

  std::string sSfM_Data_Filename;
  std::string sOutDir = "";

  cmd.add( make_option('i', sSfM_Data_Filename, "sfmdata") );
  cmd.add( make_option('o', sOutDir, "outdir") );

  try {
      if (argc == 1) throw std::string("Invalid command line parameter.");
      cmd.process(argc, argv);
  } catch(const std::string& s) {
      std::cerr << "Usage: " << argv[0] << '\n'
      << "[-i|--sfmdata filename, the SfM_Data file to convert]\n"
      << "[-o|--outdir path, (will be prefixed by Ori-)]\n"
      << std::endl;

      std::cerr << s << std::endl;
      return EXIT_FAILURE;
  }

  std::cout << " You called : " <<std::endl
            << argv[0] << std::endl
            << "--sfmdata " << sSfM_Data_Filename << std::endl
            << "--outdir " << sOutDir << std::endl;

  // Use Apero/MicMac prefix
  sOutDir = "Ori-" + sOutDir;
  // Create output dir
  if (!stlplus::folder_exists(sOutDir))
    stlplus::folder_create( sOutDir );

  // Read the SfM scene
  SfM_Data sfm_data;
  if (!Load(sfm_data, sSfM_Data_Filename, ESfM_Data(ALL))) {
    std::cerr << std::endl
      << "The input SfM_Data file \""<< sSfM_Data_Filename << "\" cannot be read." << std::endl;
    return EXIT_FAILURE;
  }

  for(Views::const_iterator iter = sfm_data.GetViews().begin();
      iter != sfm_data.GetViews().end(); ++iter)
  {
    const View * view = iter->second.get();
    if (!sfm_data.IsPoseAndIntrinsicDefined(view))
        continue;

    // Valid view, we can ask a pose & intrinsic data
    const Pose3 pose = sfm_data.GetPoseOrDie(view);
    Intrinsics::const_iterator iterIntrinsic = sfm_data.GetIntrinsics().find(view->id_intrinsic);
    const IntrinsicBase * cam = iterIntrinsic->second.get();

    if (!cameras::isPinhole(cam->getType()))
        continue;
    const Pinhole_Intrinsic * pinhole_cam = static_cast<const Pinhole_Intrinsic *>(cam);

    // We compute the 'AltiSol' an the 'Depth' param. Suboptimal for now.
    double sumDepths = 0., sumZ = 0.;
    int anInc = 0;
    for(Landmarks::const_iterator iterL = sfm_data.GetLandmarks().begin();
            iterL != sfm_data.GetLandmarks().end(); ++iterL)
    {
        for(Observations::const_iterator iterObs = iterL->second.obs.begin();
                iterObs != iterL->second.obs.end(); ++iterObs)
        {
            if(iter->second->id_view == view->id_view)
            {
                ++anInc;
                sumZ += iterL->second.X[2];
                sumDepths += pose.depth(iterL->second.X);
            }
        }
    }

    const double meanDepth = sumDepths / (double) anInc;
    const double meanZ = sumZ / (double) anInc;

    // Extrinsic
    const Vec3 c = pose.center();
    const Mat3 R = pose.rotation().transpose();
    // Intrinsic
    const double f = pinhole_cam->focal();
    const Vec2 pp = pinhole_cam->principal_point();
    // Image size in px
    const int w = pinhole_cam->w();
    const int h = pinhole_cam->h();

    // We can now create the .xml file for this View in the output dir
    std::ofstream outfile( stlplus::create_filespec(
                sOutDir, std::string("Orientation-") + stlplus::filename_part(view->s_Img_path), "xml" ).c_str() );
    // See ParamChantierPhotogram.xml in MicMac distrib for full specs. The doc is also useful !
    outfile
        << "<?xml version=\"1.0\" ?>\n"
        << "<ExportAPERO>\n"
        << "    <OrientationConique>\n"
        << "        <OrIntImaM2C>\n"
        << "            <I00>0 0</I00>\n"
        << "            <V10>1 0</V10>\n"
        << "            <V01>0 1</V01>\n"
        << "        </OrIntImaM2C>\n"
        << "        <TypeProj>eProjStenope</TypeProj>\n"
        << "        <ZoneUtileInPixel>true</ZoneUtileInPixel>\n"
        << "        <Interne>\n"
        << "            <KnownConv>eConvApero_DistM2C</KnownConv>\n"
        << "            <PP>" << pp(0) << " " << pp(1) <<"</PP>\n"
        << "            <F>" << f << "</F>\n"
        << "            <SzIm>" << w << " " << h << "</SzIm>\n"
        << "            <CalibDistortion>\n"
        << "                <ModRad>\n"
        << "                    <CDist>" << pp(0) << " " <<  pp(1) << "</CDist>\n"
        << "                    <CoeffDist>" << "</CoeffDist>\n" //TODO: loop through all Ks and nest a <CoefDist> each time
                                                                 // it's the same formula up to the normalization applied in OMVG
                                                                // so raw Ks from OMVG are to big for MM.
        << "                </ModRad>\n"
        << "            </CalibDistortion>\n"
        << "        </Interne>\n"
        << "        <Externe>\n"
        << "            <AltiSol>" << meanZ << "</AltiSol>\n"
        << "            <Profondeur>" << meanDepth << "</Profondeur>\n"
        << "            <KnownConv>eConvApero_DistM2C</KnownConv>\n"
        << "            <Centre>" << c(0) << " " << c(1) << " " << c(2) << "</Centre>\n"
        << "            <IncCentre>1 1 1</IncCentre>\n"
        << "            <ParamRotation>\n"
        << "                <CodageMatr>\n"
        << "                    <L1>" << R(0,0) << " " << R(0,1) << " " << R(0,2) <<  "</L1>\n"
        << "                    <L2>" << R(1,0) << " " << R(1,1) << " " << R(1,2) << "</L2>\n"
        << "                    <L3>" << R(2,0) << " " << R(2,1) << " " << R(2,2) << "</L3>\n"
        << "                </CodageMatr>\n"
        << "            </ParamRotation>\n"
        << "        </Externe>\n"
        << "        <ConvOri>\n"
        << "            <KnownConv>eConvApero_DistM2C</KnownConv>\n"
        << "        </ConvOri>\n"
        << "    </OrientationConique>\n"
        << "</ExportAPERO>\n";
        outfile.close();
  }

  //TODO: export landmarks or a match file as an "Homol" directory (needed to use the C3DC pipeline)
  return EXIT_SUCCESS;
}
