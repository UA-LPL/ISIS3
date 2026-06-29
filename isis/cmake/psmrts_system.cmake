# Installation of the PSMRTS system using FetchContent
include(FetchContent)

# This downloads PSMRTS into ISIS3/build/_deps
set(PSMRTS_BUILD_APPS      ON)
set(PSMRTS_BUILD_CAPI_APPS ON)
set(PSMRTS_BUILD_SHARED    ON)
FetchContent_Declare(
  psmrts
  GIT_REPOSITORY https://github.com/UA-LPL/psmrts.git
  # GIT_TAG        3712a77078ceac14b5f657a5941be77ffa84a00c
  GIT_TAG        feature/psmrts-isis-tracer
  OVERRIDE_FIND_PACKAGE
)

# Configure the PSMRTS system
FetchContent_MakeAvailable(psmrts)
