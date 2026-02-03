file(REMOVE_RECURSE
  "libpic_modules.a"
  "libpic_modules.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang Fortran)
  include(CMakeFiles/pic_modules.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
