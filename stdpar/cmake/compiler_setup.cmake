#
# Minimal compiler setup for warning flags
#
if(CMAKE_CXX_COMPILER_ID MATCHES "NVHPC" OR CMAKE_CXX_COMPILER_ID MATCHES "PGI")

  add_compile_options("-Minform=warn;-Minfo")

else()

  add_compile_options("-Wall;-Wextra")

endif()
