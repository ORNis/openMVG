// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Romain JANVIER

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_DATA_TRACKS_SELECTION_HPP
#define OPENMVG_SFM_DATA_TRACKS_SELECTION_HPP

#include "openMVG/graph/graph.hpp"

#include "openMVG/sfm/pipelines/sfm_features_provider.hpp"
#include "openMVG/sfm/pipelines/sfm_matches_provider.hpp"

#include "openMVG/sfm/sfm_data.hpp"

#include "openMVG/sfm/tracks_selection/sfm_data_track_statistics_cui.hpp"


namespace openMVG { namespace sfm { struct SfM_Data; }}

namespace openMVG {
namespace sfm {

///- Base class for implementation of tracks selection methods
class SfM_Data_Tracks_Selection_Basis {

public:
    explicit SfM_Data_Tracks_Selection_Basis(const SfM_Data &sfm_data) : sfm_data_(sfm_data) {};

    virtual Landmarks select() = 0;

protected:
    const SfM_Data &sfm_data_;
};

///- Track Selection method exposed in:
/// Cui H., Shen S., Hu. Z., Tracks Selection for Robust, Efficient and Scalable Large-Scale Structure from Motion.
/// Pattern Recognition, 2017.
class SfM_Data_Cui_Tracks_Selection : SfM_Data_Tracks_Selection_Basis {

public:
    SfM_Data_Cui_Tracks_Selection(const SfM_Data &sfm_data,
                                  const graph::indexedGraph &epipolar_graph,
                                  const std::shared_ptr<Features_Provider> & features_provider);

    Landmarks select() override;


    void setValidityCostThreshold(const double threshold)
    {
        validity_cost_threshold_ = threshold;
    }

    void setNumMSTRuns(const uint32_t runs)
    {
        num_trees_ = runs;
    }

    void setSigmaMultiplier(const double mu)
    {
        mu_ = mu;
    }

    void setPercentSelected(const double percent)
    {
        percent_selected_ = percent;
    }

private:

    void initAdjascentMap();

    void buildTrackStatistic();

    void buildViewStatistic();

    const graph::indexedGraph &epipolar_graph_;
    std::map<IndexT, std::set<IndexT >> adjascent_map_;
    const std::shared_ptr<Features_Provider> features_provider_;

    ReprojectionInvertedList reproj_invertedList_;
    TracksInvertedList tracks_invertedList_;
    ViewStats view_stats_;

    /// Parameters
    // sigma multiplier
    double mu_ = 3.0;

    // number of mst run
    uint32_t num_trees_ = 30;

    // alternative stop criterion (not in the paper)
    double percent_selected_ = 0.0;

    // a valid track must have a cost < max_image_size * validity_cost_threshold_
    double validity_cost_threshold_ = 0.3;
};



///- Track Selection method exposed in:
/// Cui H., Shen S., Gao X., Hu. Z.,
/// Batched Incremental Structure-from-Motion. International Conference on 3D Vision (3DV) 2017.
class SfM_Data_Batched_Tracks_Selection : SfM_Data_Tracks_Selection_Basis {


public:
    SfM_Data_Batched_Tracks_Selection(const SfM_Data &sfm_data);

    Landmarks select() override;

    void setCoverage(const uint32_t coverage)
    {
        coverage_ = coverage;
    }

private:

    void buildTrackStatistic();

    /// parameter
    // number of selected tracks per view
    uint32_t coverage_ = 100;

    // sorted set of tracks
    std::set<TrackStats, SimpleTrackStatComparator> track_set_;
};

} // namespace sfm
} // namespace openMVG

#endif //OPENMVG_SFM_DATA_TRACKS_SELECTION_HPP
