#
# Minimal compiler setup for warning flags
#
if(CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC" OR CMAKE_Fortran_COMPILER_ID MATCHES "NVIDIA")

  add_compile_options("-Minform=warn;-Minfo")

elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")

    add_compile_options("-Wall;-Wextra;-ffixed-line-length-none;-ffree-line-length-none")

else()

  add_compile_options("-Wall;-Wextra")

endif()
