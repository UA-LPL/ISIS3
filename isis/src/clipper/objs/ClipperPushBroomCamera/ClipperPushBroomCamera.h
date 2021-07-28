#ifndef ClipperPushBroomCamera_h
#define ClipperPushBroomCamera_h

/** This is free and unencumbered software released into the public domain.
The authors of ISIS do not claim copyright on the contents of this file.
For more details about the LICENSE terms and the AUTHORS, you will
find files of those names at the top level of this repository. **/

/* SPDX-License-Identifier: CC0-1.0 */

#include "LineScanCamera.h"

#include <QString>

#include "VariableLineScanCameraDetectorMap.h"

namespace Isis {
  /**
   * This is the camera model for the Europa Clipper Push Broom Camera
   *
   */
  class ClipperPushBroomCamera : public LineScanCamera {
    public:
      ClipperPushBroomCamera(Cube &cube);

      ~ClipperPushBroomCamera();

      virtual int CkFrameId() const;
      virtual int CkReferenceId() const;
      virtual int SpkReferenceId() const;
  };
};

#endif
