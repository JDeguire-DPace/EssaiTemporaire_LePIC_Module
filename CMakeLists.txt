cmake_minimum_required(VERSION 3.23)
project(LePIC_3D LANGUAGES Fortran)

find_package(MPI REQUIRED COMPONENTS Fortran)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM")
  add_compile_options(-O3)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  add_compile_options(-O3)
endif()

set(MODULES_DIR ${CMAKE_SOURCE_DIR}/modules)

# legacy includes like include 'mg.h'
include_directories(${MODULES_DIR}/legacy)

# ----------------------------
# Core modules (now also includes PoissonDecomp because core USEs it)
# ----------------------------
add_library(pic_core_obj OBJECT
  ${MODULES_DIR}/config/mod_config.f90
  ${MODULES_DIR}/config/mod_readConditions.f90

  ${MODULES_DIR}/core/mod_utils.f90
  ${MODULES_DIR}/core/mod_constants.f90
  ${MODULES_DIR}/core/mod_part_info.f90
  ${MODULES_DIR}/core/mod_state.f90
  ${MODULES_DIR}/core/mod_debug_checks.f90

  # >>> moved here so mod_poisson_decomp.mod exists for core users <<<
  ${MODULES_DIR}/fields/mod_PoissonDecomposition.f90
  ${MODULES_DIR}/fields/mod_rho_stub.f90

  ${MODULES_DIR}/geometry/mod_domain.f90
  ${MODULES_DIR}/geometry/mod_generateBoundary.f90
  ${MODULES_DIR}/geometry/mod_boundary.f90

  ${MODULES_DIR}/fields/mod_fields.f90
  ${MODULES_DIR}/fields/mod_charge_weights.f90

  ${MODULES_DIR}/particles/mod_chemistryState.f90
  ${MODULES_DIR}/particles/mod_particles.f90
  ${MODULES_DIR}/particles/mod_reactions.f90
  ${MODULES_DIR}/particles/mod_reactionsDB.f90

  ${MODULES_DIR}/simulation/mod_simulation.f90
)
target_link_libraries(pic_core_obj PUBLIC MPI::MPI_Fortran)
target_include_directories(pic_core_obj PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})

# ----------------------------
# Poisson stub driver (calls pdesolver_)
# ----------------------------
add_library(poisson_stub_driver_obj OBJECT
  ${MODULES_DIR}/fields/mod_poisson_stub_driver.f90
)
target_link_libraries(poisson_stub_driver_obj PUBLIC MPI::MPI_Fortran)
target_include_directories(poisson_stub_driver_obj PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})

add_dependencies(poisson_stub_driver_obj pic_core_obj)
target_include_directories(poisson_stub_driver_obj PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})
# ----------------------------
# Legacy multigrid (provides pdesolver_)
# ----------------------------
add_library(legacy_mg STATIC
  ${MODULES_DIR}/legacy/pdesolver.f90
  ${MODULES_DIR}/legacy/mg.f90
  ${MODULES_DIR}/legacy/sors.f90
)
target_link_libraries(legacy_mg PUBLIC MPI::MPI_Fortran)
target_include_directories(legacy_mg PUBLIC
  ${MODULES_DIR}/legacy
  ${CMAKE_Fortran_MODULE_DIRECTORY}
)
add_dependencies(legacy_mg pic_core_obj)

# ----------------------------
# Executables
# ----------------------------
add_executable(run_min
  ${MODULES_DIR}/main/main.f90
)
target_link_libraries(run_min PRIVATE
  pic_core_obj
  MPI::MPI_Fortran
)

add_executable(run_poisson_quadPotential
  ${MODULES_DIR}/main/main_poisson_stub.f90
)
target_link_libraries(run_poisson_quadPotential PRIVATE
  pic_core_obj
  poisson_stub_driver_obj
  legacy_mg
  MPI::MPI_Fortran
)
