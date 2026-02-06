cmake_minimum_required(VERSION 3.20)
project(LePIC_3D LANGUAGES Fortran)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

find_package(MPI REQUIRED COMPONENTS Fortran)

set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -O3")

set(PIC_SOURCES
  modules/core/mod_constants.f90
  modules/core/mod_utils.f90
  modules/core/mod_part_info.f90

  modules/config/mod_config.f90
  modules/config/mod_readConditions.f90

  modules/particles/mod_chemistryState.f90
  modules/particles/mod_reactions.f90
  modules/particles/mod_reactionsDB.f90

  modules/fields/mod_poissonDecomposition.f90
  modules/fields/mod_charge_weights.f90

  modules/geometry/mod_domain.f90
  modules/fields/mod_fields.f90
  modules/geometry/mod_generateBoundary.f90
  modules/geometry/mod_boundary.f90

  modules/core/mod_debug_checks.f90

  modules/core/mod_state.f90
)

add_library(pic_modules_obj OBJECT ${PIC_SOURCES})

target_include_directories(pic_modules_obj
  PUBLIC
    ${CMAKE_Fortran_MODULE_DIRECTORY}
)

target_link_libraries(pic_modules_obj PUBLIC MPI::MPI_Fortran)

add_executable(run_min
  modules/main.f90
  $<TARGET_OBJECTS:pic_modules_obj>
)

target_include_directories(run_min
  PUBLIC
    ${CMAKE_Fortran_MODULE_DIRECTORY}
)

target_link_libraries(run_min PUBLIC MPI::MPI_Fortran)
