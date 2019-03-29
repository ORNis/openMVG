// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2016 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MULTIVIEW_TRIANGULATION_HPP
#define OPENMVG_MULTIVIEW_TRIANGULATION_HPP

#include "openMVG/numeric/eigen_alias_definition.hpp"

namespace openMVG
{

/**
* @brief Linear DLT triangulation
* @param P1 First camera projection matrix
* @param P2 Second camera projection matrix
* @param x1 bearing vector of the landmark observation in the first camera
* @param x2 bearing vector of the landmark observation in the second camera
* @param[out] X_homogeneous Homogeneous triangulated point
* @see HZ 12.2 pag.312
* @ref Multiple View Geometry - Richard Hartley, Andrew Zisserman - second edition
*/
void TriangulateDLT
(
  const Mat34 &P1,
  const Vec3 &x1,
  const Mat34 &P2,
  const Vec3 &x2,
  Vec4 *X_homogeneous
);

/**
* @brief Linear DLT triangulation
* @param P1 First camera projection matrix
* @param P2 Second camera projection matrix
* @param x1 bearing vector of the landmark observation in the first camera
* @param x2 bearing vector of the landmark observation in the second camera
* @param[out] X_euclidean Euclidean triangulated point
* @see HZ 12.2 pag.312
* @ref Multiple View Geometry - Richard Hartley, Andrew Zisserman - second edition
*/
void TriangulateDLT
( const Mat34 &P1,
  const Vec3 &x1,
  const Mat34 &P2,
  const Vec3 &x2,
  Vec3 *X_euclidean
);

/**
* @brief Optimal L1 Angular triangulation
* @brief Minimize the L1 norm of angular errors
* @param R0 First Camera rotation matrix
* @param t0 First Camera translation vector
* @param x0 bearing vector of the landmark observation in the first camera
* @param R1 Second Camera rotation matrix
* @param t1 Second Camera translation vector
* @param x1 bearing vector of the landmark observation in the second camera
* @param[out] X_euclidean Euclidean triangulated point
* @ref S.H. Lee, J. Civera - Closed-Form Optimal Triangulation Based on Angular Errors - https://arxiv.org/pdf/1903.09115.pdf
*/
void TriangulateL1Angular
(
  const Mat3 &R0,
  const Vec3 &t0,
  const Vec3 &x0,
  const Mat3 &R1,
  const Vec3 &t1,
  const Vec3 &x1,
  Vec3 *X_euclidean
);


/**
* @brief Optimal LInfinity Angular triangulation
* @brief Minimize the LInfinity norm of angular errors
* @param R0 First Camera rotation matrix
* @param t0 First Camera translation vector
* @param x0 bearing vector of the landmark observation in the first camera
* @param R1 Second Camera rotation matrix
* @param t1 Second Camera translation vector
* @param x1 bearing vector of the landmark observation in the second camera
* @param[out] X_euclidean Euclidean triangulated point
* @ref S.H. Lee, J. Civera - Closed-Form Optimal Triangulation Based on Angular Errors - https://arxiv.org/pdf/1903.09115.pdf
*/
void TriangulateLInfinityAngular
(
  const Mat3 &R0,
  const Vec3 &t0,
  const Vec3 &x0,
  const Mat3 &R1,
  const Vec3 &t1,
  const Vec3 &x1,
  Vec3 *X_euclidean
);


} // namespace openMVG

#endif  // OPENMVG_MULTIVIEW_TRIANGULATION_HPP
