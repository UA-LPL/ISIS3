# Installation of the PSMRTS system using FetchContent
include(FetchContent)

# This downloads PSMRTS into ISIS3/build/_deps
set(PSMRTS_BUILD_APPS      ON)
set(PSMRTS_BUILD_CAPI_APPS ON)
set(PSMRTS_BUILD_SHARED    ON)
FetchContent_Declare(
  psmrts
  GIT_REPOSITORY https://github.com/UA-LPL/psmrts.git
  GIT_TAG        b85e10f56a0e3617be3f145cff93862706049491
  OVERRIDE_FIND_PACKAGE
)

# Configure the PSMRTS system
FetchContent_MakeAvailable(psmrts)
