// Copyright (c) 2016 Pierre MOULON & Romain JANVIER

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef OPENMVG_FEATURES_HOMOGRAPHY_DECORATOR_IMAGE_DESCRIBER_HPP
#define OPENMVG_FEATURES_HOMOFRAPHY_IMAGE_DESCRIBER_HPP

#include "openMVG/features/image_describer.hpp"
#include "openMVG/features/regions_factory.hpp"
#include "openMVG/numeric/numeric.h"
#include <cereal/cereal.hpp>


namespace openMVG {
namespace features {

class Homography_Decorator_Image_describer : public Image_describer
{

public:

  Homography_Decorator_Image_describer() {}; //FIXME ?

  Homography_Decorator_Image_describer(std::unique_ptr<Image_describer> id): img_describer_(std::move(id)) {
    H_.setIdentity();
  };

  Homography_Decorator_Image_describer(std::unique_ptr<Image_describer> id, const Mat3 & H): img_describer_(std::move(id)), H_(H) {};

  bool Set_configuration_preset(EDESCRIBER_PRESET preset) {
    return img_describer_->Set_configuration_preset(preset);
  };

  bool Describe(const image::Image<unsigned char> & image,
                std::unique_ptr<Regions> &regions,
                const image::Image<unsigned char> * mask = nullptr ) {

    image::Image<unsigned char> transformedImage;
    image::Image<unsigned char> transformedMask;
    bool bOk;

    Vec2 topLeft(0, 0);
    Vec2 topRight(image.Width()-1, 0);
    Vec2 bottomRight(image.Width()-1, image.Height()-1);
    Vec2 bottomLeft(0, image.Height()-1);

    // Check for mirror transformation
    // Miror transformations are non sense here
    if(!(image::ApplyH_AndCheckOrientation(H_, topLeft.x(), topLeft.y()) &&
       image::ApplyH_AndCheckOrientation(H_,topRight.x(), topRight.y()) &&
       image::ApplyH_AndCheckOrientation(H_, bottomRight.x(), bottomRight.y()) &&
       image::ApplyH_AndCheckOrientation(H_, bottomLeft.x(), bottomLeft.y()))) return false;

    // Computing dimensions of the transformed image
    int width = static_cast<int>(floor(max(max(max(topLeft.x(), topRight.x()), bottomLeft.x()), bottomRight.x()))) +1;
    int height = static_cast<int>(floor(max(max(max(topLeft.y(), topRight.y()), bottomLeft.y()), bottomRight.y()))) +1;

    // Warping
    transformedImage.resize(width, height);
    image::Warp(image, H_.inverse(), transformedImage);

    if(mask)
    {
      transformedMask.resize(width, height);
      image::Warp(*mask, H_.inverse(), transformedMask);
      bOk = img_describer_->Describe(transformedImage, regions, mask);
    } else {
      bOk = img_describer_->Describe(transformedImage, regions, nullptr);
    }

    if(!bOk)
      return false;

    regions->ApplyTransformToRegionPositions(H_.inverse(), image.Width(), image.Height());
    return true;
  }

  void Allocate(std::unique_ptr<Regions> &regions) const {
    img_describer_->Allocate(regions);
  };

protected:
  std::unique_ptr<Image_describer> img_describer_;
  Mat3 H_;
};

} // namespace features
} // namespace openMVG
#endif // OPENMVG_FEATURES_HOMOGRAPHY_DECORATOR_IMAGE_DESCRIBER_HPP