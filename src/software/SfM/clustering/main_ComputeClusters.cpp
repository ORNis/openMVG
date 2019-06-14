// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2016 Pierre Moulon, Srivathsan Murali, Romuald Perrot.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/image/image_io.hpp"
#include "openMVG/sfm/sfm.hpp"
#include "openMVG/system/timer.hpp"
#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <algorithm>
#include <cstdlib>
#include <iomanip>
#include <iterator>
#include <memory>
#include <string>
#include <vector>

#include <opencv2/ml/ml.hpp>

#include "domSetLibrary/domset.h"
#include "domSetLibrary/types.h"

#ifdef OPENMVG_USE_OPENMP
#include <omp.h>
#endif

using namespace openMVG;
using namespace openMVG::cameras;
using namespace openMVG::geometry;
using namespace openMVG::image;
using namespace openMVG::sfm;


template<unsigned int N>
class BucketImage
{
public:
  BucketImage(int width = 0, int height = 0): height_(height), width_(width){
    size_cell_x_ = width_ / double(NUM_CELL_ONE_DIM);
    size_cell_y_ = height_ / double(NUM_CELL_ONE_DIM);
  }
  
  void insert(const Vec2 &pos, const Vec3 & point)
  {
    size_t idx = std::floor(pos.x() / size_cell_x_);
    size_t idy = std::floor(pos.y() / size_cell_y_);
    points_in_cells[idx][idy].push_back(point);
  }

  std::vector<Vec3> get_mean_points(unsigned int threshold) const
  {
    std::vector<Vec3> result;

    for(int i = 0; i < NUM_CELL_ONE_DIM; ++i)
    {
      for(int j = 0; j < NUM_CELL_ONE_DIM; ++j)
      {
        if(points_in_cells[i][j].size() < threshold)
          continue;
        
        Vec3 acc = Vec3::Zero();
        for(const auto & pt : points_in_cells[i][j])
        {
          acc += pt;
        }
        result.push_back(acc /= points_in_cells[i][j].size());
      }
    }
    return result;
  }

private:
  static constexpr unsigned int NUM_CELL_LEVEL = 4;
  static constexpr unsigned int NUM_CELL_TOTAL = std::pow(NUM_CELL_LEVEL, N);
  static constexpr unsigned int NUM_CELL_ONE_DIM = std::pow(NUM_CELL_LEVEL/2, N);
  std::array<std::array<std::vector<Vec3>,  NUM_CELL_ONE_DIM>, NUM_CELL_ONE_DIM> points_in_cells;
  int width_, height_;
  double size_cell_x_, size_cell_y_;
};

/**
* @brief convert openMVG sfm_data to domset data
* @param sfm_data openMVG dataset
* @param[out] cameras list of intrinsics
* @param[out] view list of views
* @param[out] points list of points
* @param[out] map_view map between views added and the global index
*  We need this map to handle the fact that OpenMVG view Id can be non
*   contiguous:
*   - in case of missing pose or intrinsic for some view
* @retval true if success
* @retval false if failure
*/
bool domsetImporter(
    const SfM_Data &sfm_data,
    std::vector<nomoko::Camera> &cameras,
    std::vector<nomoko::View> &views,
    std::vector<nomoko::Point> &points,
    std::map<openMVG::IndexT, uint32_t> &map_view)
{

  // Convert OpenMVG data to domset library data
  openMVG::system::Timer loadDataTimer;

  // adding views
  for (const auto &view : sfm_data.GetViews())
  {
    if (sfm_data.IsPoseAndIntrinsicDefined(view.second.get()))
    {
      map_view[view.first] = views.size(); // need to make global

      const openMVG::geometry::Pose3 poseMVG(sfm_data.GetPoseOrDie(view.second.get()));
      nomoko::View v;
      v.rot = poseMVG.rotation().cast<float>();
      v.trans = poseMVG.center().transpose().cast<float>();
      views.push_back(v);
    }
  }

  // adding landmarks
  for (const auto &it_landmark : sfm_data.GetLandmarks())
  {
    const Landmark &landmark = it_landmark.second;
    const Observations &obs = landmark.obs;
    // X, color, obsCount
    std::vector<size_t> vIds;
    for (const auto &it_obs : obs)
    {
      vIds.push_back(map_view[it_obs.first]);
    }

    nomoko::Point p;
    p.pos = landmark.X.transpose().cast<float>();
    p.viewList = vIds;
    points.push_back(p);
  }

  std::cout << std::endl
            << "Number of views  = " << views.size() << std::endl
            << "Number of points = " << points.size() << std::endl
            << "Loading data took (s): "
            << loadDataTimer.elapsed() << std::endl;
  return true;
}

/**
* @brief Export a sfm_data file using a subset of the view of a given sfm_data
* @param sfm_data The whole data set
* @param outFilename Output file name
* @param cluster List of view to consider
* @retval true if success
* @retval false if failure
*/
bool exportData(const SfM_Data &sfm_data,
                const std::string &outFilename,
                const std::set<size_t> &cluster)
{
  SfM_Data cl_sfm_data;
  cl_sfm_data.s_root_path = sfm_data.s_root_path;

  // Copy the view (only the requested ones)
  for (const auto view : sfm_data.GetViews())
  {
    const bool inCluster = cluster.find(view.first) != cluster.end();
    if (inCluster && sfm_data.IsPoseAndIntrinsicDefined(view.second.get()))
    {
      cl_sfm_data.poses[view.first] = sfm_data.GetPoseOrDie(view.second.get());
      cl_sfm_data.views[view.first] = view.second;

      const auto intrinsic = sfm_data.GetIntrinsics().at(view.second.get()->id_intrinsic);
      if (cl_sfm_data.intrinsics.count(view.second.get()->id_intrinsic) == 0)
        cl_sfm_data.intrinsics[view.second.get()->id_intrinsic] = intrinsic;
    }
  }

  // Copy observations that have relation with the considered view
  for (const auto &it_landmark : sfm_data.GetLandmarks())
  {
    const Landmark &landmark = it_landmark.second;
    Observations obs;
    for (const auto &observation : landmark.obs)
    {
      const auto &it = cl_sfm_data.views.find(observation.first);
      if (it != cl_sfm_data.views.end())
      {
        obs[observation.first] = observation.second;
      }
    }
    // Landmark observed in less than 2 view are ignored
    if (obs.size() < 2)
      continue;
    cl_sfm_data.structure[it_landmark.first].X = landmark.X;
    cl_sfm_data.structure[it_landmark.first].obs = obs;
  }

  return Save(cl_sfm_data, outFilename, ESfM_Data(ALL));
}

int main(int argc, char **argv)
{
  using namespace std;
  std::cout << "Dominant Set Clustering" << std::endl
            << std::endl;

  CmdLine cmd;

  std::string sSfM_Data_Filename = "";
  std::string sOutDir = "";
  unsigned int clusterSizeLowerBound = 20;
  unsigned int clusterSizeUpperBound = 30;
  unsigned int clusterOverlap = 4;
  float voxelGridSize = 0.01f;

  cmd.add(make_option('i', sSfM_Data_Filename, "input_file"));
  cmd.add(make_option('o', sOutDir, "outdir"));
  cmd.add(make_option('l', clusterSizeLowerBound, "cluster_size_lower_bound"));
  cmd.add(make_option('u', clusterSizeUpperBound, "cluster_size_upper_bound"));
  cmd.add(make_option('a', clusterOverlap, "cluster_overlap"));
  cmd.add(make_option('v', voxelGridSize, "voxel_grid_size"));

  try
  {
    if (argc == 1)
      throw std::string("Invalid command line parameter.");
    cmd.process(argc, argv);
  }
  catch (const std::string &s)
  {
    std::cerr << "Usage: " << argv[0] << "\n"
              << "[-i|--input_file] path to a SfM_Data scene\n"
              << "[-o|--outdir path] path to output directory\n"
              << "[-l|--cluster_size_lower_bound] lower bound to cluster size (20)\n"
              << "[-u|--cluster_size_upper_bound] upper bound to cluster size (30)\n"
              << "[-a|--cluter_overlap] number of common images between adjascent clusters (4) \n" //TODO test 0
              << "[-v|--voxel_grid_size] voxel grid size (0.01)\n"
              //TODO if opencv...
              << std::endl;

    std::cerr << s << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Params: " << argv[0] << std::endl
            << "[Input file]      = " << sSfM_Data_Filename << std::endl
            << "[Outdir path]     = " << sOutDir << std::endl
            << "[Cluster size:" << std::endl
            << "    Lower bound   = " << clusterSizeLowerBound << std::endl
            << "    Upper bound]  = " << clusterSizeUpperBound << std::endl
            << "[Cluster overlap] = " << clusterOverlap << std::endl
            << "[Voxel grid size] = " << voxelGridSize << std::endl;

  if (sSfM_Data_Filename.empty())
  {
    std::cerr << "\nIt is an invalid file input" << std::endl;
    return EXIT_FAILURE;
  }

  // Prepare output folder
  if (!stlplus::folder_exists(sOutDir))
    if (!stlplus::folder_create(sOutDir))
    {
      std::cerr << "Cannot create: " << sOutDir << std::endl;
      return EXIT_FAILURE;
    }

  SfM_Data sfm_data;
  if (!Load(sfm_data, sSfM_Data_Filename, ESfM_Data(ALL)))
  {
    std::cerr << std::endl
              << "The input SfM_Data file \"" << sSfM_Data_Filename << "\" can't be read."
              << std::endl;
    return EXIT_FAILURE;
  }

  // loading data
  std::vector<nomoko::Camera> cameras; // stores the various camera intrinsic parameters
  std::vector<nomoko::View> views;     // stores the poses for each view
  std::vector<nomoko::Point> points;   // 3d point positions

  std::map<openMVG::IndexT, uint32_t> origViewMap; // need to keep track of original views ids
  if (!domsetImporter(sfm_data, cameras, views, points, origViewMap))
  {
    std::cerr << "Error: can't import data" << std::endl;
    return EXIT_FAILURE;
  }

  //---------------------------------------
  // clustering views process
  //---------------------------------------
  openMVG::system::Timer clusteringTimer;

  nomoko::Domset domset(points, views, cameras, voxelGridSize);
  domset.clusterViews(clusterSizeLowerBound, clusterSizeUpperBound, clusterOverlap);

  std::cout << "Clustering view took (s): "
            << clusteringTimer.elapsed() << std::endl;

  // export to ply to visualize
  const std::string viewOut = sOutDir + "/views.ply";
  domset.exportToPLY(viewOut, true);
  const std::string viewOutBorders = sOutDir + "/viewBorders.ply";
  domset.exportToPLYBorders(viewOutBorders);

  // Retrieve the cluster and export them
  std::vector<std::set<size_t>> finalClusters;
  std::vector<std::set<size_t>> finalClustersNoOverlap;
  {
    const std::vector<std::vector<size_t>> clusters = domset.getClusters();
    const std::vector<std::vector<size_t>> clustersNoOverlap = domset.getClustersNoOverlap();

    // Remap the camera index from contiguous to original view Id
    std::map<openMVG::IndexT, uint32_t> origViewMap_reverse;
    for (const auto &it : origViewMap)
    {
      origViewMap_reverse.insert(std::make_pair(it.second, it.first));
    }
    // For every cluster, remap the view Id
    for (const auto &cl : clusters)
    {
      std::set<size_t> newCl;
      for (const auto vId : cl)
        newCl.insert(origViewMap_reverse[vId]);
      finalClusters.emplace_back(newCl);
    }
    int i = 0;

    // For every cluster, remap the view Id //TODO
    for (const auto &cl : clustersNoOverlap)
    {
      std::set<size_t> newCl;
      for (const auto vId : cl)
        newCl.insert(origViewMap_reverse[vId]);

      finalClustersNoOverlap.emplace_back(newCl);
    }
  }

  const size_t numClusters = finalClusters.size();
  std::cout << "Number of clusters = " << numClusters << std::endl;

  //SVM
  std::map<size_t, size_t> map_view_to_cl;

  for (int i = 0; i < finalClustersNoOverlap.size(); ++i)
  {
    for (const auto id : finalClustersNoOverlap[i])
    {
      map_view_to_cl[id] = i;
    }
  }

  std::map<IndexT, std::vector<Vec3>> map_cl_to_points;
  std::map<IndexT, BucketImage<1>> map_views_to_grid;

  for(const auto & v : sfm_data.GetViews())
  {
    map_views_to_grid[v.first] = BucketImage<1>(v.second->ui_width, v.second->ui_height);
  }

  for (const auto &track : sfm_data.GetLandmarks()) // gridded points
  {
    for (const auto &obs : track.second.obs)
    {
      map_views_to_grid[obs.first].insert(obs.second.x, track.second.X);
    }
  }
 
  int total_points = 0;
  for (const auto &view : map_view_to_cl)
  {
    map_cl_to_points[view.second].emplace_back(
        sfm_data.poses.at(sfm_data.views.at(view.first)->id_pose).center());
    ++total_points;
    for(const auto & pts : map_views_to_grid[view.first].get_mean_points(5))
    {
      map_cl_to_points[view.second].emplace_back(pts);
      ++total_points;
    }
  }

  cv::Mat trainingDataMat = cv::Mat::zeros(total_points, 3, CV_32FC1);
  cv::Mat labelsMat = cv::Mat::zeros(total_points, 1, CV_32SC1);

  int idx = 0;
  for (const auto &cl : map_cl_to_points)
  {
    for (const auto &pt : cl.second)
    {
      labelsMat.at<int>(idx) = cl.first;
      trainingDataMat.at<float>(idx, 0) = static_cast<float>(pt.x());
      trainingDataMat.at<float>(idx, 1) = static_cast<float>(pt.y());
      trainingDataMat.at<float>(idx, 2) = static_cast<float>(pt.z());
      ++idx;
    }
  }

  std::cout << "[ Computing SVM parameters ] " << std::endl;
  auto svm = cv::ml::SVM::create();
  svm->setType(cv::ml::SVM::C_SVC);
  svm->setKernel(cv::ml::SVM::RBF);
  svm->setTermCriteria(cv::TermCriteria(cv::TermCriteria::MAX_ITER, 1000, 1e-8));
  auto train_data = cv::ml::TrainData::create(trainingDataMat, cv::ml::ROW_SAMPLE, labelsMat);
  svm->trainAuto(train_data);
  svm->save("svm_param.yml");

#ifdef OPENMVG_USE_OPENMP
#pragma omp parallel for
#endif
  for (int i = 0; i < numClusters; ++i)
  {
    std::stringstream filename;
    filename << sOutDir << "/sfm_data";
    filename << std::setw(4) << std::setfill('0') << i;
    filename << ".bin";

#ifdef OPENMVG_USE_OPENMP
#pragma omp critical
#endif
    {
      std::stringstream ss;
      ss << "Writing cluster to " << filename.str() << std::endl;
      std::cout << ss.str();
    }

    if (!exportData(sfm_data, filename.str(), finalClusters[i]))
    {
      std::stringstream str;
      str << "Could not write cluster : " << filename.str() << std::endl;
      std::cerr << str.str();
    }
  }
  return EXIT_SUCCESS;
}
