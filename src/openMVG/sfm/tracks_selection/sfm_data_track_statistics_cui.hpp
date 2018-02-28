// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Romain JANVIER

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_DATA_TRACK_STATISTICS_CUI_HPP
#define OPENMVG_SFM_DATA_TRACK_STATISTICS_CUI_HPP

#include "openMVG/types.hpp"

#include <memory>

namespace openMVG {
namespace sfm {

struct TrackStats {

    TrackStats(const IndexT &id_,
               const size_t size_,
               const double scale_,
               const double cost_,
               const std::set<IndexT> &visibility_set_);

    friend std::ostream &operator<<(std::ostream &Stream, const TrackStats &ts);

    // Is it an already selected track?
    bool selected = false;

    std::set<IndexT> visibility_set;

    // Id of the track (id in the Landmarks map)
    IndexT id = 0;

    // Number of obs : it is basically visibility_set.size() but we cache it
    size_t size = 0;

    // Scale of the feature
    double scale = 0.0;

    // cost
    double reprojection_cost = 0.0;
};


struct TrackStatPointerComparator {
    // /!\ sign are inverted because we do not want to reverse order
    bool operator()(const std::shared_ptr<TrackStats> &lhs, const std::shared_ptr<TrackStats> &rhs) const {
      if (rhs->size < lhs->size)
        return true;
      else if (rhs->size > lhs->size)
        return false;

      if (lhs->scale > rhs->scale)
        return true;
      else if (lhs->scale < rhs->scale)
        return false;

      if (lhs->reprojection_cost < rhs->reprojection_cost)
        return true;
      else if (lhs->reprojection_cost > rhs->reprojection_cost)
        return false;
    }
};

// view_id, set of TrackStat
using TracksInvertedList = std::map<IndexT, std::set<std::shared_ptr<TrackStats>, TrackStatPointerComparator> >;

// view_id, reprojection_error
using ReprojectionInvertedList = std::map<IndexT, std::vector<double> >; //TODO: per view reproj error

// view_id & cost
using ViewStats = std::vector<std::pair<double, IndexT> >;

// helpers
// TODO: should be placed somewhere else?
template<typename Container>
double mean_iterable(const Container &data) {
  const double sum = std::accumulate(std::begin(data), std::end(data), 0.0);
  return sum / data.size();
}

// compute reprojection cost of a track as mean(data) + mu * sigma(data)
template<typename Container>
double reprojection_cost_iterable(const Container &data, const double mu) {
  const double mean = mean_iterable(data);
  const double sq_sum = std::inner_product(std::begin(data), std::end(data), std::begin(data), 0.0);
  const double sigma = std::sqrt(sq_sum / data.size() - mean * mean);

  return mean + mu * sigma;
}


} // namespace sfm
} // namespace openMVG

#endif //OPENMVG_SFM_DATA_TRACK_STATISTICS_CUI_HPP
