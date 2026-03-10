// ! ==========================================================
// ! MINIMAL particle_info.h shim for legacy load_part.f90
// ! Modern build: globals are provided by mod_legacy_particle_globals.
// ! So this header MUST NOT declare np_cell, vt0, Nm, Pabs, etc.
// ! ==========================================================

//       integer, parameter :: npart = 32

//       real(kind=8)     :: charge(npart), mass(npart), Ti(npart)
//       character(len=6) :: pname(npart)

//       data charge / npart * 0.0d0 /
//       data mass   / npart * 1.0d0 /
//       data Ti     / npart * 0.0d0 /
//       data pname  / npart * '      ' /