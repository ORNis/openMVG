// Copyright (c) 2016 Pierre MOULON & Romain JANVIER

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef OPENMVG_FEATURES_SCALING_DECORATOR_IMAGE_DESCRIBER_HPP
#define OPENMVG_FEATURES_SCALING_DECORATOR_IMAGE_DESCRIBER_HPP

#include "openMVG/features/image_describer.hpp"
#include "openMVG/features/regions_factory.hpp"
#include "openMVG/image/image.hpp"
#include "openMVG/numeric/numeric.h"

namespace openMVG {
namespace features {

class Scaling_Decorator_Image_describer : public Homography_Decorator_Image_describer
{

public:
  Scaling_Decorator_Image_describer(std::unique_ptr<Image_describer> id, double scale): Homography_Decorator_Image_describer(std::move(id)) {
    Mat3 H;
    H << scale, 0, 0,
        0, scale, 0,
        0, 0, 1;
    this->H_ = H;
  };

};

} // namespace features
} // namespace openMVG
#endif // OPENMVG_FEATURES_SCALING_DECORATOR_IMAGE_DESCRIBER_HPP