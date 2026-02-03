cmake_minimum_required(VERSION 3.20)
project(LePIC_3D LANGUAGES Fortran)

# Put .mod files in build/mod
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

# Default build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

# MPI (works with mpiifx / mpifort wrappers too)
find_package(MPI REQUIRED COMPONENTS Fortran)

# --- Fortran flags (keep simple; add -qopenmp if you actually use OMP) ---
set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -O3")

# Collect module sources
set(PIC_SOURCES
  modules/core/mod_constants.f90
  modules/core/mod_utils.f90

  modules/config/mod_config.f90
  modules/config/mod_readConditions.f90

  modules/geometry/mod_domain.f90
  modules/geometry/mod_generateBoundary.f90
  modules/geometry/mod_boundary.f90

  modules/core/mod_state.f90
)

# Build modules as an object library (nice for big codes)
add_library(pic_modules_obj OBJECT ${PIC_SOURCES})

target_include_directories(pic_modules_obj
  PUBLIC
    ${CMAKE_Fortran_MODULE_DIRECTORY}
)

target_link_libraries(pic_modules_obj PUBLIC MPI::MPI_Fortran)

# Main executable
add_executable(run_min
  modules/main.f90
  $<TARGET_OBJECTS:pic_modules_obj>
)

target_include_directories(run_min
  PUBLIC
    ${CMAKE_Fortran_MODULE_DIRECTORY}
)

target_link_libraries(run_min PUBLIC MPI::MPI_Fortran)
