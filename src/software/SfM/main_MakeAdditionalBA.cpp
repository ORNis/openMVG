
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/sfm.hpp"
#include "openMVG/sfm/sfm_data_utils.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

using namespace openMVG;
using namespace openMVG::sfm;


int main(int argc, char **argv)
{
  using namespace std;
  std::cout << "Compute Additional BA(s)" << std::endl;

  CmdLine cmd;


  std::string sSfM_Data_Filename;
  std::string sOutFile = "";
  float weight = 0.0f;

  cmd.add( make_option('i', sSfM_Data_Filename, "input_file") );
  cmd.add( make_option('o', sOutFile, "output_file") );
  cmd.add( make_option('w', weight, "weight") );
  cmd.add( make_switch('s', "split"));

  try {
    if (argc == 1) throw std::string("Invalid command line parameter.");
    cmd.process(argc, argv);
  } catch(const std::string& s) {
    std::cerr << "Usage: " << argv[0] << '\n'
    << "[-i|--input_file] path to a SfM_Data scene\n"
    << "[-o|--output_file] file where the output data will be stored "
    <<    "(i.e. path/sfm_data_structure.bin)\n"
    << "[-w|--weight] weight of the CGPs in the BA (0.0 by default: GCPs are not used)\n "
    << "\n[Optional]\n"
    << "[-s|--split] (switch) split grouped intrinsic before BA (OFF by default)\n"
    << std::endl;

    std::cerr << s << std::endl;
    return EXIT_FAILURE;
  }


  // Load input SfM_Data scene
  SfM_Data sfm_data;
  if (!Load(sfm_data, sSfM_Data_Filename, ESfM_Data(ALL))) {
    std::cerr << std::endl
    << "The input SfM_Data file \""<< sSfM_Data_Filename << "\" cannot be read." << std::endl;
    return EXIT_FAILURE;
  }

  std::cout
    << "Loaded a sfm_data scene with:\n"
    << " #views: " << sfm_data.GetViews().size() << "\n"
    << " #poses: " << sfm_data.GetPoses().size() << "\n"
    << " #intrinsics: " << sfm_data.GetIntrinsics().size() <<  "\n"
    << " #tracks: " << sfm_data.GetLandmarks().size()
    << std::endl;

  if(cmd.used('s')) {
    SplitSharedIntrinsics(sfm_data);
  }

  std::cout << "Bundle adjustment..." << std::endl;

  Bundle_Adjustment_Ceres bundle_adjustment_obj;
  Control_Point_Parameter control_point_opt;

  if(weight > 0.0f) {
    control_point_opt.bUse_control_points = true;
    control_point_opt.weight = weight;
  }

  bundle_adjustment_obj.Adjust
    (
      sfm_data,
      Optimize_Options(
        cameras::Intrinsic_Parameter_Type::ADJUST_ALL,
        Extrinsic_Parameter_Type::ADJUST_ALL,
        Structure_Parameter_Type::ADJUST_ALL,
        control_point_opt
      )
    );

  if (stlplus::extension_part(sOutFile) != "ply") {
    Save(sfm_data,
      stlplus::create_filespec(
        stlplus::folder_part(sOutFile),
        stlplus::basename_part(sOutFile), "ply"),
      ESfM_Data(ALL));
  }

  if (Save(sfm_data, sOutFile, ESfM_Data(ALL)))
    return EXIT_SUCCESS;
  else
    return EXIT_FAILURE;
}
