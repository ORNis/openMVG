// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Romain JANVIER

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/tracks_selection/sfm_data_track_statistics_cui.hpp"
namespace openMVG {
namespace sfm {

TrackStats::TrackStats(IndexT id_,
                       size_t size_,
                       double scale_,
                       double reprojection_cost_,
                       const std::set<IndexT> &visibility_set_) :
        id(id_), size(size_), reprojection_cost(reprojection_cost_), visibility_set(visibility_set_) {
  // It's not likely that scale == rhs.scale (because it is an high precision float value) so we truncate it
  scale = std::round(scale_ * 10.0) / 10.0;
};

TrackStats::TrackStats(IndexT id_,
                       size_t size_,
                       double reprojection_cost_,
                       const std::set<IndexT> &visibility_set_) :
        id(id_), size(size_), reprojection_cost(reprojection_cost_), visibility_set(visibility_set_) {};

std::ostream &operator<<(std::ostream &Stream, const TrackStats &ts) {
  Stream << ts.id << " size: " << ts.size << " scale: " << ts.scale << " cost: " << ts.reprojection_cost;
  return Stream;
}

} // namespace sfm
} // namespace openMVG
