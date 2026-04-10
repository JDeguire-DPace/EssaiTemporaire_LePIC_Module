cmake_minimum_required(VERSION 3.23)
project(LePIC_3D LANGUAGES Fortran)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

# ------------------------------------------------------------
# Detect whether the selected Fortran compiler is an MPI wrapper
# ------------------------------------------------------------
get_filename_component(FORTRAN_COMPILER_NAME "${CMAKE_Fortran_COMPILER}" NAME)

set(USING_MPI_WRAPPER FALSE)
if(FORTRAN_COMPILER_NAME MATCHES "^(mpiifx|mpiifort|mpifort|mpif90|mpifc)$")
  set(USING_MPI_WRAPPER TRUE)
endif()

message(STATUS "Fortran compiler          : ${CMAKE_Fortran_COMPILER}")
message(STATUS "Fortran compiler ID       : ${CMAKE_Fortran_COMPILER_ID}")
message(STATUS "Fortran compiler basename : ${FORTRAN_COMPILER_NAME}")
message(STATUS "Using MPI wrapper         : ${USING_MPI_WRAPPER}")

# ------------------------------------------------------------
# MPI handling
# ------------------------------------------------------------
if(USING_MPI_WRAPPER)
  message(STATUS "MPI will be provided by compiler wrapper: ${FORTRAN_COMPILER_NAME}")
else()
  find_package(MPI REQUIRED COMPONENTS Fortran)
  message(STATUS "MPI_Fortran_FOUND       : ${MPI_Fortran_FOUND}")
  message(STATUS "MPI_Fortran_COMPILER    : ${MPI_Fortran_COMPILER}")
  message(STATUS "MPI_Fortran_LIBRARIES   : ${MPI_Fortran_LIBRARIES}")
  message(STATUS "MPI_Fortran_INCLUDE_DIRS: ${MPI_Fortran_INCLUDE_DIRS}")
endif()

# ------------------------------------------------------------
# Compiler / OpenMP flags
# ------------------------------------------------------------
set(OMP_COMPILE_FLAG "")
set(OMP_LINK_FLAG "")

if(CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM|Intel")
  set(OMP_COMPILE_FLAG "-qopenmp")
  set(OMP_LINK_FLAG    "-qopenmp")
  add_compile_options(-O3 ${OMP_COMPILE_FLAG})
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(OMP_COMPILE_FLAG "-fopenmp")
  set(OMP_LINK_FLAG    "-fopenmp")
  add_compile_options(-O3 ${OMP_COMPILE_FLAG})
else()
  message(WARNING
    "Unknown Fortran compiler: ${CMAKE_Fortran_COMPILER_ID}. "
    "OpenMP flags may need adjustment."
  )
  add_compile_options(-O3)
endif()

# ------------------------------------------------------------
# Common interface target
# ------------------------------------------------------------
add_library(build_config INTERFACE)

target_include_directories(build_config INTERFACE
  ${CMAKE_Fortran_MODULE_DIRECTORY}
)

if(OMP_LINK_FLAG)
  target_link_options(build_config INTERFACE ${OMP_LINK_FLAG})
endif()

if(NOT USING_MPI_WRAPPER)
  target_link_libraries(build_config INTERFACE MPI::MPI_Fortran)
  target_include_directories(build_config INTERFACE ${MPI_Fortran_INCLUDE_DIRS})
endif()

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------
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
target_link_libraries(core_deps PUBLIC build_config)
target_include_directories(core_deps PUBLIC
  ${CMAKE_Fortran_MODULE_DIRECTORY}
  ${MODULES_DIR}/legacy
)

# ============================================================
# 1) Config base
# ============================================================
add_library(config_base STATIC
  ${MODULES_DIR}/config/mod_config.f90
)
target_link_libraries(config_base PUBLIC build_config core_deps)
target_include_directories(config_base PUBLIC
  ${CMAKE_Fortran_MODULE_DIRECTORY}
)

# ============================================================
# 2) Geometry base
# ============================================================
add_library(geometry_base STATIC
  ${MODULES_DIR}/geometry/mod_domain.f90
)
target_link_libraries(geometry_base PUBLIC
  build_config
  core_deps
  config_base
)
target_include_directories(geometry_base PUBLIC
  ${CMAKE_Fortran_MODULE_DIRECTORY}
)
add_dependencies(geometry_base config_base)

# ============================================================
# 3) Fields base
# ============================================================
add_library(fields_base STATIC
  ${MODULES_DIR}/fields/mod_fields.f90
  ${MODULES_DIR}/fields/mod_density.f90
)
target_link_libraries(fields_base PUBLIC
  build_config
  core_deps
  config_base
  geometry_base
)
target_include_directories(fields_base PUBLIC
  ${CMAKE_Fortran_MODULE_DIRECTORY}
)
add_dependencies(fields_base geometry_base)

# ============================================================
# 4) Core object library
# ============================================================
add_library(pic_core_obj OBJECT
  # -------------------------
  # Config / input
  # -------------------------
  ${MODULES_DIR}/config/mod_readConditions.f90

  # -------------------------
  # Core
  # -------------------------
  ${MODULES_DIR}/core/mod_intro.f90
  ${MODULES_DIR}/core/mod_part_info.f90
  ${MODULES_DIR}/core/mod_debug_checks.f90

  # -------------------------
  # Fields
  # -------------------------
  ${MODULES_DIR}/fields/mod_PoissonDecomposition.f90
  ${MODULES_DIR}/fields/mod_PoissonSolver_legacy.f90
  ${MODULES_DIR}/fields/mod_charge_weights.f90
  ${MODULES_DIR}/fields/mod_chargeDeposition.f90
  ${MODULES_DIR}/fields/mod_magneticField.f90
  ${MODULES_DIR}/fields/mod_electricField.f90

  # -------------------------
  # Geometry / boundaries
  # -------------------------
  ${MODULES_DIR}/geometry/mod_generateBoundary.f90
  ${MODULES_DIR}/geometry/mod_boundary.f90

  # -------------------------
  # Particles / chemistry
  # -------------------------
  ${MODULES_DIR}/particles/mod_chemistryState.f90
  ${MODULES_DIR}/particles/mod_particles.f90
  ${MODULES_DIR}/particles/mod_particle_sorting.f90
  ${MODULES_DIR}/particles/mod_particleMover.f90
  ${MODULES_DIR}/particles/mod_particleBC.f90
  ${MODULES_DIR}/particles/mod_particle_loader.f90
  ${MODULES_DIR}/particles/mod_collisions.f90
  ${MODULES_DIR}/particles/mod_reactions.f90
  ${MODULES_DIR}/particles/mod_reactionsDB.f90
  ${MODULES_DIR}/particles/mod_heating.f90
  ${MODULES_DIR}/particles/mod_planeMoments.f90

  # -------------------------
  # Simulation
  # -------------------------
  ${MODULES_DIR}/simulation/mod_simParams.f90

  # -------------------------
  # IO
  # -------------------------
  ${MODULES_DIR}/io/mod_output_2d.f90

  # -------------------------
  # High-level orchestration
  # -------------------------
  ${MODULES_DIR}/core/mod_state.f90
  ${MODULES_DIR}/simulation/mod_simulation.f90
)
target_link_libraries(pic_core_obj PUBLIC
  build_config
  core_deps
  config_base
  geometry_base
  fields_base
)
target_include_directories(pic_core_obj PUBLIC
  ${CMAKE_Fortran_MODULE_DIRECTORY}
)
add_dependencies(pic_core_obj
  core_deps
  config_base
  geometry_base
  fields_base
)

# ============================================================
# 5) Legacy multigrid / Poisson core
# ============================================================
add_library(legacy_mg STATIC
  ${MODULES_DIR}/legacy/pdesolver.f90
  ${MODULES_DIR}/legacy/mg.f90
  ${MODULES_DIR}/legacy/sors.f90
)
target_link_libraries(legacy_mg PUBLIC
  build_config
)
target_include_directories(legacy_mg PUBLIC
  ${MODULES_DIR}/legacy
  ${CMAKE_Fortran_MODULE_DIRECTORY}
)
add_dependencies(legacy_mg pic_core_obj)

# ============================================================
# 6) Main executable
# ============================================================
add_executable(run_min
  ${MODULES_DIR}/main/main.f90
)

target_link_libraries(run_min PRIVATE
  build_config
  pic_core_obj
  fields_base
  geometry_base
  config_base
  core_deps
  legacy_mg
)

target_include_directories(run_min PRIVATE
  ${CMAKE_Fortran_MODULE_DIRECTORY}
  ${MODULES_DIR}/legacy
)

set_target_properties(run_min PROPERTIES
  OUTPUT_NAME "run_min"
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
)

message(STATUS "Executable will be: ${CMAKE_BINARY_DIR}/run_min")