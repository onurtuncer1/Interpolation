/* ---------------------------------------------------------------------------*
  Copyright 2022 MILTEKSAN

  Use of this software is restricted to MILTEKSAN

  Written by Melina Aero, Istanbul, Turkey
  Contact onur.tuncer@melina-aero.com
------------------------------------------------------------------------------*/

#ifndef POSE_H
#define POSE_H

#include <glm/glm.hpp>

struct Pose {

    glm::dvec3 Point;
    glm::quaternion Orientation;
}

#endif // POSE_H