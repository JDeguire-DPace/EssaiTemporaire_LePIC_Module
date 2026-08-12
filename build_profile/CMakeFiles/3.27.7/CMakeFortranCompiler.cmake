set(CMAKE_Fortran_COMPILER "/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v4/Compiler/intel2025/intelmpi/2021.16.1/mpi/2021.16/bin/mpiifx")
set(CMAKE_Fortran_COMPILER_ARG1 "")
set(CMAKE_Fortran_COMPILER_ID "IntelLLVM")
set(CMAKE_Fortran_COMPILER_VERSION "2025.2.0")
set(CMAKE_Fortran_COMPILER_WRAPPER "")
set(CMAKE_Fortran_PLATFORM_ID "Linux")
set(CMAKE_Fortran_SIMULATE_ID "")
set(CMAKE_Fortran_COMPILER_FRONTEND_VARIANT "GNU")
set(CMAKE_Fortran_SIMULATE_VERSION "")




set(CMAKE_AR "/cvmfs/soft.computecanada.ca/gentoo/2023/x86-64-v3/usr/bin/ar")
set(CMAKE_Fortran_COMPILER_AR "CMAKE_Fortran_COMPILER_AR-NOTFOUND")
set(CMAKE_RANLIB "/cvmfs/soft.computecanada.ca/gentoo/2023/x86-64-v3/usr/bin/ranlib")
set(CMAKE_TAPI "CMAKE_TAPI-NOTFOUND")
set(CMAKE_Fortran_COMPILER_RANLIB "CMAKE_Fortran_COMPILER_RANLIB-NOTFOUND")
set(CMAKE_COMPILER_IS_GNUG77 )
set(CMAKE_Fortran_COMPILER_LOADED 1)
set(CMAKE_Fortran_COMPILER_WORKS TRUE)
set(CMAKE_Fortran_ABI_COMPILED TRUE)

set(CMAKE_Fortran_COMPILER_ENV_VAR "FC")

set(CMAKE_Fortran_COMPILER_SUPPORTS_F90 1)

set(CMAKE_Fortran_COMPILER_ID_RUN 1)
set(CMAKE_Fortran_SOURCE_FILE_EXTENSIONS f;F;fpp;FPP;f77;F77;f90;F90;for;For;FOR;f95;F95;f03;F03;f08;F08)
set(CMAKE_Fortran_IGNORE_EXTENSIONS h;H;o;O;obj;OBJ;def;DEF;rc;RC)
set(CMAKE_Fortran_LINKER_PREFERENCE 20)
set(CMAKE_Fortran_LINKER_DEPFILE_SUPPORTED )
if(UNIX)
  set(CMAKE_Fortran_OUTPUT_EXTENSION .o)
else()
  set(CMAKE_Fortran_OUTPUT_EXTENSION .obj)
endif()

# Save compiler ABI information.
set(CMAKE_Fortran_SIZEOF_DATA_PTR "8")
set(CMAKE_Fortran_COMPILER_ABI "ELF")
set(CMAKE_Fortran_LIBRARY_ARCHITECTURE "x86_64-unknown-linux-gnu")

if(CMAKE_Fortran_SIZEOF_DATA_PTR AND NOT CMAKE_SIZEOF_VOID_P)
  set(CMAKE_SIZEOF_VOID_P "${CMAKE_Fortran_SIZEOF_DATA_PTR}")
endif()

if(CMAKE_Fortran_COMPILER_ABI)
  set(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_Fortran_COMPILER_ABI}")
endif()

if(CMAKE_Fortran_LIBRARY_ARCHITECTURE)
  set(CMAKE_LIBRARY_ARCHITECTURE "x86_64-unknown-linux-gnu")
endif()





set(CMAKE_Fortran_IMPLICIT_INCLUDE_DIRECTORIES "/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v4/Compiler/intel2025/intelmpi/2021.16.1/mpi/2021.16/include/mpi;/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v4/Compiler/intel2025/intelmpi/2021.16.1/mpi/2021.16/include;/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v4/Compiler/gcccore/ucx/1.19.0/include;/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v4/Compiler/intel2025/aocl-lapack/5.1/include;/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v4/Compiler/intel2025/aocl-blas/5.1/include/blis;/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v4/Compiler/intel2025/aocl-blas/5.1/include;/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v4/Compiler/intel2025/flexiblas/3.4.5/include/flexiblas;/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/Core/imkl/2025.2.0/mkl/2025.2/include/fftw;/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/Core/imkl/2025.2.0/mkl/2025.2/include;/cvmfs/restricted.computecanada.ca/easybuild/software/2023/x86-64-v3/Core/intel/2025.2.0/tbb/2022.2/include;/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v4/Compiler/gcccore/python/3.11.5/include;/cvmfs/restricted.computecanada.ca/easybuild/software/2023/x86-64-v3/Core/intel/2025.2.0/compiler/2025.2/opt/compiler/include/intel64;/cvmfs/restricted.computecanada.ca/easybuild/software/2023/x86-64-v3/Core/intel/2025.2.0/compiler/2025.2/opt/compiler/include;/cvmfs/restricted.computecanada.ca/easybuild/software/2023/x86-64-v3/Core/intel/2025.2.0/compiler/2025.2/include;/cvmfs/restricted.computecanada.ca/easybuild/software/2023/x86-64-v3/Core/intel/2025.2.0/compiler/2025.2/lib/clang/21/include;/cvmfs/soft.computecanada.ca/gentoo/2023/x86-64-v3/usr/include")
set(CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES "mpifort;mpi;dl;rt;pthread;ifport;ifcoremt;imf;svml;m;ipgo;irc;pthread;svml;dl;c;gcc;gcc_s;irc_s;dl;c")
set(CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES "/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v4/Compiler/intel2025/intelmpi/2021.16.1/mpi/2021.16/lib;/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v4/Compiler/intel2025/intelmpi/2021.16.1/mpi/2021.16/opt/mpi/libfabric/lib;/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v4/Compiler/gcccore/ucx/1.19.0/lib;/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v4/Compiler/intel2025/aocl-lapack/5.1/lib;/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v4/Compiler/intel2025/aocl-blas/5.1/lib;/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v4/Compiler/intel2025/flexiblas/3.4.5/lib64;/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/Core/imkl/2025.2.0/mkl/2025.2/lib;/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/Core/imkl/2025.2.0/compiler/2025.2/lib;/cvmfs/restricted.computecanada.ca/easybuild/software/2023/x86-64-v3/Core/intel/2025.2.0/tbb/2022.2/lib;/cvmfs/restricted.computecanada.ca/easybuild/software/2023/x86-64-v3/Core/intel/2025.2.0/compiler/2025.2/lib;/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v4/Compiler/gcccore/python/3.11.5/lib;/cvmfs/restricted.computecanada.ca/easybuild/software/2023/x86-64-v3/Core/intel/2025.2.0/compiler/2025.2/lib/clang/21/lib/x86_64-unknown-linux-gnu;/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v4/Compiler/gcccore/ucx/1.19.0/lib64;/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v4/Compiler/intel2025/aocl-lapack/5.1/lib64;/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v4/Compiler/intel2025/aocl-blas/5.1/lib64;/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v4/Compiler/gcccore/python/3.11.5/lib64;/cvmfs/soft.computecanada.ca/gentoo/2023/x86-64-v3/usr/lib/gcc/x86_64-pc-linux-gnu/14;/cvmfs/soft.computecanada.ca/gentoo/2023/x86-64-v3/usr/lib64;/cvmfs/soft.computecanada.ca/gentoo/2023/x86-64-v3/lib64;/cvmfs/soft.computecanada.ca/gentoo/2023/x86-64-v3/usr/x86_64-pc-linux-gnu/lib;/cvmfs/soft.computecanada.ca/gentoo/2023/x86-64-v3/usr/lib;/cvmfs/soft.computecanada.ca/gentoo/2023/x86-64-v3/lib")
set(CMAKE_Fortran_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")
