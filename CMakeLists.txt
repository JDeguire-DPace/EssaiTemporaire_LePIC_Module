cmake_minimum_required(VERSION 3.23)
project(LePIC_3D LANGUAGES Fortran)

find_package(MPI REQUIRED COMPONENTS Fortran)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

# ----------------------------
# Compiler + OpenMP flags
# ----------------------------
if(CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM|Intel")
  add_compile_options(-O3 -qopenmp)
  set(OMP_LINK_FLAG "-qopenmp")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  add_compile_options(-O3 -fopenmp)
  set(OMP_LINK_FLAG "-fopenmp")
else()
  message(WARNING "Unknown Fortran compiler: ${CMAKE_Fortran_COMPILER_ID}. OpenMP may not link correctly.")
  set(OMP_LINK_FLAG "")
endif()

set(MODULES_DIR ${CMAKE_SOURCE_DIR}/modules)

include_directories(${MODULES_DIR}/legacy)

# ============================================================
# 0) Core dependency library
# ============================================================
add_library(core_deps STATIC
  ${MODULES_DIR}/core/mod_constants.f90
  ${MODULES_DIR}/core/mod_utils.f90
  ${MODULES_DIR}/core/mod_RNG.f90
)
target_link_libraries(core_deps PUBLIC MPI::MPI_Fortran)
target_include_directories(core_deps PUBLIC
  ${CMAKE_Fortran_MODULE_DIRECTORY}
  ${MODULES_DIR}/legacy
)
target_link_options(core_deps PRIVATE ${OMP_LINK_FLAG})

# ============================================================
# 1) Config base first
# ============================================================
add_library(config_base STATIC
  ${MODULES_DIR}/config/mod_config.f90
)
target_link_libraries(config_base PUBLIC MPI::MPI_Fortran core_deps)
target_include_directories(config_base PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})
target_link_options(config_base PRIVATE ${OMP_LINK_FLAG})

# ============================================================
# 2) Geometry base
# ============================================================
add_library(geometry_base STATIC
  ${MODULES_DIR}/geometry/mod_domain.f90
)
target_link_libraries(geometry_base PUBLIC MPI::MPI_Fortran core_deps config_base)
target_include_directories(geometry_base PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})
target_link_options(geometry_base PRIVATE ${OMP_LINK_FLAG})
add_dependencies(geometry_base config_base)

# ============================================================
# 3) Fields base
# ============================================================
add_library(fields_base STATIC
  ${MODULES_DIR}/fields/mod_fields.f90
  ${MODULES_DIR}/fields/mod_density.f90
)
target_link_libraries(fields_base PUBLIC MPI::MPI_Fortran core_deps config_base geometry_base)
target_include_directories(fields_base PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})
target_link_options(fields_base PRIVATE ${OMP_LINK_FLAG})
add_dependencies(fields_base geometry_base)

# ============================================================
# 4) Core object library
# ============================================================
add_library(pic_core_obj OBJECT
  ${MODULES_DIR}/config/mod_readConditions.f90

  ${MODULES_DIR}/core/mod_intro.f90
  ${MODULES_DIR}/core/mod_part_info.f90
  ${MODULES_DIR}/core/mod_debug_checks.f90
  ${MODULES_DIR}/core/mod_state.f90

  ${MODULES_DIR}/fields/mod_PoissonDecomposition.f90
  ${MODULES_DIR}/fields/mod_rho_stub.f90
  ${MODULES_DIR}/fields/mod_charge_weights.f90
  ${MODULES_DIR}/fields/mod_magneticField.f90

  ${MODULES_DIR}/geometry/mod_generateBoundary.f90
  ${MODULES_DIR}/geometry/mod_boundary.f90

  ${MODULES_DIR}/particles/mod_chemistryState.f90
  ${MODULES_DIR}/particles/mod_particles.f90
  ${MODULES_DIR}/particles/mod_particle_loader.f90
  ${MODULES_DIR}/particles/mod_reactions.f90
  ${MODULES_DIR}/particles/mod_reactionsDB.f90

  ${MODULES_DIR}/simulation/mod_simParams.f90
  ${MODULES_DIR}/simulation/mod_simulation.f90

  ${MODULES_DIR}/io/mod_output_2d.f90
)
target_link_libraries(pic_core_obj PUBLIC
  MPI::MPI_Fortran
  core_deps
  config_base
  geometry_base
  fields_base
)
target_include_directories(pic_core_obj PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})
add_dependencies(pic_core_obj core_deps config_base geometry_base fields_base)

# ----------------------------
# Poisson stub driver object
# ----------------------------
add_library(poisson_stub_driver_obj OBJECT
  ${MODULES_DIR}/fields/mod_poisson_stub_driver.f90
)
target_link_libraries(poisson_stub_driver_obj PUBLIC
  MPI::MPI_Fortran
  core_deps
  config_base
  geometry_base
  fields_base
)
target_include_directories(poisson_stub_driver_obj PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})
add_dependencies(poisson_stub_driver_obj pic_core_obj)

# ----------------------------
# Legacy multigrid
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
target_link_options(legacy_mg PRIVATE ${OMP_LINK_FLAG})

# ----------------------------
# Executables
# ----------------------------
add_executable(run_min
  ${MODULES_DIR}/main/main.f90
)
target_link_libraries(run_min PRIVATE
  pic_core_obj
  fields_base
  geometry_base
  config_base
  core_deps
  MPI::MPI_Fortran
)
target_link_options(run_min PRIVATE ${OMP_LINK_FLAG})

add_executable(run_poisson_quadPotential
  ${MODULES_DIR}/main/main_poisson_stub.f90
)
target_link_libraries(run_poisson_quadPotential PRIVATE
  pic_core_obj
  fields_base
  geometry_base
  config_base
  core_deps
  poisson_stub_driver_obj
  legacy_mg
  MPI::MPI_Fortran
)
target_link_options(run_poisson_quadPotential PRIVATE ${OMP_LINK_FLAG})
