# CMake module for find_package(BulletFloat64) double precision package
# Finds include directory and all applicable libraries
#
#[=[
This Bullet find module provides a target definition for the double precision
version of the Bullet library. The target is provided as an imported library
interface. It is derived from CMAKE varibles provided in the Bullet package
configuration. The Bullet package is assumed to be provided in the Conda 
environment but could be used by other compatible Bullet configurations that
provide a double precision version the Bullet package.

Note the only thing that differentiates the Conda package manager version is
the explicit find_package() call that is specific to Conda. In the Conda
Bullet package, the double precision package configuration is provided in 
a separate CMAKE file, namely BulletConfig-float64.cmake. Other package
managers use different techniques. However any package manager that provides
consistent definitions of the Bullet CMAKE variables can use the 
add_bullet_double_target() macro to create the Bullet::Bullet_double target.

This module provides the following elements:

add_bullet_double_target()   # Macro to generate the bullet Target
Bullet::Bullet_double        # INTERFACE IMPORTED library target
Bullet_double_FOUND          # Module found variable

To add this target to ISIS, just append the Bullet::Bullet_double target
to ALLLIBS variable in ISIS3/isis/CMakeLists.txt file.
]=]

macro(add_bullet_double_target)

  if (NOT TARGET Bullet::Bullet_double)
    set(Bullet_double_FOUND FALSE)
    if (BULLET_FOUND OR Bullet_FOUND )
      if(NOT ${BULLET_DEFINITIONS} MATCHES ".*-DBT_USE_DOUBLE_PRECISION.*")
        message(FATAL_ERROR "Bullet does not appear to be built with double "
                  "precision, current definitions: ${BULLET_DEFINITIONS}")
      endif()
      message(STATUS "Bullet Compile Definitions: ${BULLET_DEFINITIONS}")

      # This configuration ensures the Bullet variable definitions are
      # also conformant.
      add_library(Bullet::Bullet_double INTERFACE IMPORTED)
      set_target_properties(Bullet::Bullet_double
        PROPERTIES
          INTERFACE_COMPILE_DEFINITIONS "${BULLET_DEFINITIONS}"
          INTERFACE_INCLUDE_DIRECTORIES "${BULLET_ROOT_DIR}/${BULLET_INCLUDE_DIRS}"
          INTERFACE_LINK_DIRECTORIES "${BULLET_ROOT_DIR}/${BULLET_LIBRARY_DIR}"
          INTERFACE_LINK_LIBRARIES "${BULLET_LIBRARIES}"
      )
      set(Bullet_double_FOUND TRUE)
      message(STATUS "Bullet Target Created: Bullet::Bullet_double")
    endif()
  endif()

endmacro()

# Call the Conda specific Bullet double precision package
find_package(Bullet REQUIRED CONFIGS BulletConfig-float64.cmake)

# Create/confirm the Conda Bullet double precision interface target.
# Note this target needs to be added to the end of the ALLLIBS variable
# in ISIS3/isis/CMakeLists.txt file.
add_bullet_double_target()
