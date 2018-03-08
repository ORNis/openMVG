// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Romain JANVIER

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/tracks_selection/sfm_data_tracks_selection.hpp"

#include "openMVG/sfm/pipelines/sfm_features_provider.hpp"
#include "openMVG/system/timer.hpp"


namespace openMVG {
namespace sfm {

SfM_Data_Cui_Tracks_Selection::SfM_Data_Cui_Tracks_Selection(const SfM_Data &sfm_data,
                                                             const graph::indexedGraph &epipolar_graph,
                                                             const std::shared_ptr<Features_Provider> &features_provider)
        :
        SfM_Data_Tracks_Selection_Basis(sfm_data), epipolar_graph_(epipolar_graph),
        features_provider_(features_provider) {}


Landmarks SfM_Data_Cui_Tracks_Selection::select() {

  initAdjascentMap();
  buildTrackStatistic();
  buildViewStatistic();

  std::set<IndexT> selected_tracks;
  std::cout << "MST computation" << std::endl;
  system::Timer timer;
  double curr_percent_selected = 0.0;
  //TODO: this could potentially turn into an infinite loop
  for (int curr_iter_mst = 0;
       curr_iter_mst < num_trees_ || curr_percent_selected <= percent_selected_; ++curr_iter_mst) {
    size_t num_views = view_stats_.size();
    ViewStats curr_layer;
    std::set<IndexT> selected_views;
    curr_layer.push_back(view_stats_[0]); // at the top of the tree we insert the view with the best cost
    selected_views.insert(view_stats_[0].second);

    while (!curr_layer.empty() && selected_views.size() != num_views) {
      ViewStats next_layer;
      std::sort(std::begin(curr_layer), std::end(curr_layer));;
      for (auto &curr_viewStat : curr_layer) {
        const std::set<IndexT> &adj_views = adjascent_map_[curr_viewStat.second];
        // We compute the diff between adjascent views and already selected views
        std::set<IndexT> valid_adj_views;
        std::set_difference(std::begin(adj_views), std::end(adj_views),
                            std::begin(selected_views), std::end(selected_views),
                            std::inserter(valid_adj_views, std::begin(valid_adj_views)));

        // if the diff is empty we go to the next view in the layer
        if (valid_adj_views.empty())
          continue;

        for (auto &curr_track : tracks_invertedList_[curr_viewStat.second]) {
          if (!curr_track->selected) {

            const std::set<IndexT> &curr_vis_set = curr_track->visibility_set;

            std::set<IndexT> valid_observed_views;

            std::set_intersection(std::begin(curr_vis_set), std::end(curr_vis_set),
                                  std::begin(valid_adj_views), std::end(valid_adj_views),
                                  std::inserter(valid_observed_views, std::begin(valid_observed_views)));

            if (valid_observed_views.empty()) // In this case the track cover none of the valid_ngb_views of the curr_view
              continue;

            // diff valid_views and selected_views to get new views covered by the track
            std::vector<IndexT> new_views;
            std::set_difference(std::begin(valid_observed_views), std::end(valid_observed_views),
                                std::begin(selected_views), std::end(selected_views),
                                std::back_inserter(new_views));

            if (new_views.empty()) // there is no additional views covered by the track
              continue;

            // The track is valid, insert the track and the new views
            selected_tracks.insert(curr_track->id);
            curr_track->selected = true;
            for (const auto &nv : new_views) {
              selected_views.insert(nv);
              valid_adj_views.erase(nv);
              auto nv_viewstat = std::find_if(std::begin(view_stats_), std::end(view_stats_),
                                              [&nv](std::pair<double, IndexT> &v_stat) -> bool {
                                                  return v_stat.second == nv;
                                              });
              // TODO: make more checks
              next_layer.push_back(*nv_viewstat);
            }
          }
          // TODO: the break semantic is maybe not the best, but for now it works
          if (valid_adj_views.empty())
            break;

          if (selected_views.size() == num_views)
            break;
        }

        if (selected_views.size() == num_views)
          break;
      }
      std::swap(curr_layer, next_layer);
    }
    curr_percent_selected = selected_tracks.size() / static_cast<double>(sfm_data_.structure.size()) * 100;

  }
  std::cout << "-- End MST computation in " << timer.elapsedMs() << " ms" << "\n";
  std::cout << " --> num selected tracks " << selected_tracks.size() << "\n";
  std::cout << " --> % selected tracks " << curr_percent_selected << std::endl;

  Landmarks result;

  for (const auto track_id : selected_tracks) {
    result[track_id] = sfm_data_.structure.at(track_id);
  }

  return result;
}


void SfM_Data_Cui_Tracks_Selection::initAdjascentMap() {
  using GraphT = lemon::ListGraph;
  // TODO filter orphans => views without pose and intrinsics // connected components // Externalize that?
  // TODO prefer and edge based iteration and a cache with vertex flagged as invalid
  const GraphT &graph = epipolar_graph_.g;
  const GraphT::NodeMap <IndexT> *node_map = epipolar_graph_.node_map_id.get();

  adjascent_map_.clear();
  for (GraphT::NodeIt nIter(graph); nIter != lemon::INVALID; ++nIter) {
    IndexT view_id = (*node_map)[nIter];
    std::set<IndexT> curr_node_set;
    for (GraphT::IncEdgeIt itNgb = lemon::ListGraph::IncEdgeIt(graph, nIter);
         itNgb != lemon::INVALID;
         ++itNgb) {
      GraphT::Node ngb_node = graph.oppositeNode(nIter, itNgb);
      curr_node_set.insert((*node_map)[ngb_node]);
    }
    adjascent_map_[view_id] = curr_node_set;
  }
}

void SfM_Data_Cui_Tracks_Selection::buildTrackStatistic() {
  // clean up previous statistics
  reproj_invertedList_.clear();
  tracks_invertedList_.clear();

  std::cout << "\n"
            << "Tracks statistic computation" << std::endl;

  system::Timer timer;
  uint32_t count_gross_outliers = 0;

  for (auto indexed_landmark = sfm_data_.structure.cbegin();
       indexed_landmark != sfm_data_.structure.cend(); ++indexed_landmark) {

    const auto &curr_landmark = indexed_landmark->second;
    const auto &observations = curr_landmark.obs;

    // Init vectors
    size_t track_size = observations.size();
    std::vector<double> feature_scales;
    feature_scales.reserve(track_size);
    std::vector<double> reprojection_errors;
    reprojection_errors.reserve(track_size);
    std::vector<uint32_t> max_size_image;
    max_size_image.reserve(track_size);
    std::set<IndexT> visibility_set;


    for (const auto &curr_obs : observations) {
      const IndexT &obs_view_id = curr_obs.first;
      const View *obs_view = sfm_data_.GetViews().at(obs_view_id).get();
      const Pose3 &obs_pose = sfm_data_.GetPoseOrDie(obs_view);
      const cameras::IntrinsicBase *obs_intrinsics = sfm_data_.GetIntrinsics().at(obs_view->id_intrinsic).get();
      //const features::PointFeature & obs_feature = features_provider_->getFeatures(obs_view_id).at(curr_obs.second.id_feat);
      const double scale = 1.0;
      // features_provider_->getScales(obs_view_id).at(curr_obs.second.id_feat);
      // TODO: scales handling is implemented in EgoSfM branch but it's not elegant (see the Slack discussion about that)

      max_size_image.push_back(std::max(obs_view->ui_height, obs_view->ui_width));

      feature_scales.push_back(scale);

      double reproj_error = obs_intrinsics->residual(obs_pose(curr_landmark.X), curr_obs.second.x).norm();
      reprojection_errors.push_back(reproj_error);

      reproj_invertedList_[obs_view_id].push_back(reproj_error);

      visibility_set.insert(obs_view_id);
    }

    double cost = reprojection_cost_iterable(reprojection_errors, mu_);


    if (cost < validity_cost_threshold_ * mean_iterable(max_size_image)) {
      std::shared_ptr<TrackStats> curr_track_stat(
              new TrackStats(indexed_landmark->first, track_size, mean_iterable(feature_scales), cost,
                             visibility_set));

      for (const auto &visible_view_id : visibility_set) {
        tracks_invertedList_[visible_view_id].insert(curr_track_stat);
      }

    } else {
      ++count_gross_outliers;
    }
  }

  std::cout << "-- End Tracks statistics computation in " << timer.elapsedMs() << " ms" << "\n";
  std::cout << " --> Found " << count_gross_outliers << " outliers" << std::endl;
}


void SfM_Data_Cui_Tracks_Selection::buildViewStatistic() {
  std::cout << "View statistics computation" << std::endl;
  // Compute viewStats
  for (const auto &reprojIL : reproj_invertedList_) {
    size_t num_adj_views = adjascent_map_[reprojIL.first].size();
    view_stats_.emplace_back(
            -(num_adj_views + std::exp(-reprojection_cost_iterable(reprojIL.second, mu_))),
            reprojIL.first);
  }
  // sort viewstats by cost
  std::sort(std::begin(view_stats_), std::end(view_stats_));
  std::cout << "-- End View statistics computation" << std::endl;
}


SfM_Data_Batched_Tracks_Selection::SfM_Data_Batched_Tracks_Selection(
        const SfM_Data &sfm_data,
        const std::unique_ptr<openMVG::tracks::SharedTrackVisibilityHelper> &shared_track_visibility_helper) :
        SfM_Data_Tracks_Selection_Basis(sfm_data), shared_track_visibility_helper_(shared_track_visibility_helper) {};


void SfM_Data_Batched_Tracks_Selection::buildTrackStatistic() {

  std::map<IndexT, TrackStats> map_track_stat;

  std::cout << "\n"
            << "Tracks statistics computation" << std::endl;
  system::Timer timer;

  // 0- We populate the coverage map with already calibrated cameras
  for (const auto &view_it : sfm_data_.GetViews()) {
    const View *v = view_it.second.get();
    if (sfm_data_.IsPoseAndIntrinsicDefined(v))
      map_view_coverage_[v->id_view] = 0;
  }

  //TODO: openMP
  // 1- We compute the Track stats for inliers tracks
  for (auto indexed_landmark = sfm_data_.structure.cbegin();
       indexed_landmark != sfm_data_.structure.cend(); ++indexed_landmark) {

    const auto &curr_landmark = indexed_landmark->second;
    const auto &observations = curr_landmark.obs;

    size_t track_size = observations.size();
    std::vector<double> reprojection_errors;
    reprojection_errors.reserve(track_size);
    std::set<IndexT> visibility_set;

    for (const auto &curr_obs : observations) {
      const IndexT &obs_view_id = curr_obs.first;
      const View *obs_view = sfm_data_.GetViews().at(obs_view_id).get();
      const Pose3 &obs_pose = sfm_data_.GetPoseOrDie(obs_view);
      const cameras::IntrinsicBase *obs_intrinsics = sfm_data_.GetIntrinsics().at(obs_view->id_intrinsic).get();
      double reproj_error = obs_intrinsics->residual(obs_pose(curr_landmark.X), curr_obs.second.x).norm();
      reprojection_errors.push_back(reproj_error);
      visibility_set.insert(obs_view_id);
    }

    // The reprojection cost is computed in the fashion of the Track selection ( mean(reproj) + mu * sigma(reproj) )
    // instead of the simple mean(reproj) exposed in Batched incremental paper, it's clearly better
    // TODO: this could be a parameter
    map_track_stat.insert({indexed_landmark->first, TrackStats(indexed_landmark->first, track_size,
                                                                reprojection_cost_iterable(reprojection_errors, mu_),
                                                                visibility_set)});
  }

  //2- we compute possible next candidate to resection
  // and we augment each TrackStats

  // Code from Pierre Moulon
  // Collect the views that does not have any 3D pose
  const std::set<IndexT> view_with_no_pose = [&] {
      std::set<IndexT> idx;
      for (const auto &view_it : sfm_data_.GetViews()) {
        const View *v = view_it.second.get();
        const IndexT id_pose = v->id_pose;
        if (sfm_data_.GetPoses().count(id_pose) == 0)
          idx.insert(view_it.first);
      }
      return idx;
  }();

  const IndexT pose_before = sfm_data_.GetPoses().size();

  // Get the track ids of the reconstructed landmarks //TODO can be computed upstream in 1-
  const std::vector<IndexT> reconstructed_trackId = [&] {
      std::vector<IndexT> tracks_ids;
      tracks_ids.reserve(sfm_data_.GetLandmarks().size());
      std::transform(sfm_data_.GetLandmarks().cbegin(), sfm_data_.GetLandmarks().cend(),
                     std::back_inserter(tracks_ids),
                     stl::RetrieveKey());
      std::sort(tracks_ids.begin(), tracks_ids.end());
      return tracks_ids;
  }();

  // List the view that have a sufficient 2D-3D coverage for robust pose estimation
#pragma omp parallel
  for (const auto &view_id : view_with_no_pose) {
#ifdef OPENMVG_USE_OPENMP
#pragma omp single nowait
#endif
    {
      // List the track related to the current view_id
      openMVG::tracks::STLMAPTracks view_tracks;
      shared_track_visibility_helper_->GetTracksInImages({view_id}, view_tracks);
      std::set<IndexT> view_tracks_ids;
      tracks::TracksUtilsMap::GetTracksIdVector(view_tracks, &view_tracks_ids);

      // Get the ids of the already reconstructed tracks
      const std::set<IndexT> track_id_for_resection = [&] {
          std::set<IndexT> track_id;
          std::set_intersection(view_tracks_ids.cbegin(), view_tracks_ids.cend(),
                                reconstructed_trackId.cbegin(), reconstructed_trackId.cend(),
                                std::inserter(track_id, track_id.begin()));
          return track_id;
      }();
      const double track_ratio = track_id_for_resection.size() / static_cast<float>(view_tracks_ids.size() + 1);
      if (track_ratio > track_ratio_ && track_id_for_resection.size() > 12) { // 12 = avoid to weakly supported views
#pragma omp critical
        {
          // the view is valid, the id is added to the coverage map.
          map_view_coverage_[view_id] = 0;
          for (auto track_id : track_id_for_resection) {
            //Beware, original paper increase also the size of the track, but tests shows that it's not a good idea
            map_track_stat.at(track_id).visibility_set.insert(view_id);
          }
        }
      }
    }
  }

  //3- Eventually we sort all TrackStats in a set
  for (auto ts : map_track_stat) {
    track_set_.insert(ts.second);
  }

  std::cout << "-- End Tracks statistics computation in " << timer.elapsedMs() << " ms" << "\n";

}


Landmarks SfM_Data_Batched_Tracks_Selection::select() {
  Landmarks selected_tracks;
  track_set_.clear();
  map_view_coverage_.clear();

  if (sfm_data_.GetLandmarks().empty())
    return selected_tracks;

  buildTrackStatistic();

  std::cout << "Track selection computation" << std::endl;
  system::Timer timer;

  uint32_t num_view_covered = 0;
  const uint64_t num_view_to_cover = map_view_coverage_.size();
  for (const auto &curr_track : track_set_) //TODO: while
  {
    bool selectable = false;
    for (const auto &view_id : curr_track.visibility_set) {
      if (map_view_coverage_[view_id] < coverage_) {
        selectable = true;
        if (++map_view_coverage_[view_id] >= coverage_)
          ++num_view_covered;
      }
    }
    if (selectable)
      selected_tracks.insert({curr_track.id, sfm_data_.structure.at(curr_track.id)});
    if (num_view_covered == num_view_to_cover)
      break;
  }
  std::cout << "-- End track selection in " << timer.elapsedMs() << " ms" << "\n";
  std::cout << " --> num selected tracks " << selected_tracks.size() << std::endl;

  return selected_tracks;
}

} // namespace sfm
} // namespace openMVG
