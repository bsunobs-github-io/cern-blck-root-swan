set(ROOT_ARCHITECTURE linux)
set(ROOT_PLATFORM linux)

if(CMAKE_SYSTEM_PROCESSOR MATCHES x86_64)
  message(STATUS "Found a 64bit system")
  if(CMAKE_COMPILER_IS_GNUCXX)
    message(STATUS "Found GNU compiler collection")
    set(ROOT_ARCHITECTURE linuxx8664gcc)
  elseif(CMAKE_CXX_COMPILER_ID STREQUAL Clang)
    message(STATUS "Found CLANG compiler")
    set(ROOT_ARCHITECTURE linuxx8664gcc)
  elseif(CMAKE_CXX_COMPILER_ID STREQUAL Intel)
    set(ROOT_ARCHITECTURE linuxx8664icc)
  else()
    message(FATAL_ERROR "There is no Setup for this compiler up to now. Don't know what to do. Stop cmake at this point.")
  endif()
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES aarch64)
  message(STATUS "Found a 64bit ARM system")
  if(CMAKE_COMPILER_IS_GNUCXX)
    message(STATUS "Found GNU compiler collection")
    set(ROOT_ARCHITECTURE linuxarm64)
  elseif(CMAKE_CXX_COMPILER_ID STREQUAL Clang)
    message(STATUS "Found CLANG compiler")
    set(ROOT_ARCHITECTURE linuxarm64)
  else()
    message(FATAL_ERROR "There is no Setup for this compiler up to now. Don't know what to do. Stop cmake at this point.")
  endif()
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES ppc64)
  message(STATUS "Found a 64bit PPC system (ppc64/ppc64le)")
  if(CMAKE_COMPILER_IS_GNUCXX)
    message(STATUS "Found GNU compiler collection")
    set(ROOT_ARCHITECTURE linuxppc64gcc)
  elseif(CMAKE_CXX_COMPILER_ID STREQUAL Clang)
    message(STATUS "Found CLANG compiler")
    set(ROOT_ARCHITECTURE linuxppc64gcc)
  else()
    message(FATAL_ERROR "There is no Setup for this compiler up to now. Don't know what to do. Stop cmake at this point.")
  endif()
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES arm)
  message(STATUS "Found a 32bit ARM system")
  if(CMAKE_COMPILER_IS_GNUCXX)
    message(STATUS "Found GNU compiler collection")
    set(ROOT_ARCHITECTURE linuxarm)
  elseif(CMAKE_CXX_COMPILER_ID STREQUAL Clang)
    message(STATUS "Found CLANG compiler")
    set(ROOT_ARCHITECTURE linuxarm)
  else()
    message(FATAL_ERROR "There is no Setup for this compiler up to now. Don't know what to do. Stop cmake at this point.")
  endif()
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES s390x)
  message(STATUS "Found a 64bit system")
  if(CMAKE_COMPILER_IS_GNUCXX)
    message(STATUS "Found GNU compiler collection")
    set(ROOT_ARCHITECTURE linuxs390xgcc)
  elseif(CMAKE_CXX_COMPILER_ID STREQUAL Clang)
    message(STATUS "Found CLANG compiler")
    set(ROOT_ARCHITECTURE linuxs390xgcc)
  else()
    message(FATAL_ERROR "There is no Setup for this compiler up to now. Don't know what to do. Stop cmake at this point.")
  endif()
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES s390)
  message(STATUS "Found a 31bit system")
  if(CMAKE_COMPILER_IS_GNUCXX)
    message(STATUS "Found GNU compiler collection")
    set(ROOT_ARCHITECTURE linuxs390gcc)
  elseif(CMAKE_CXX_COMPILER_ID STREQUAL Clang)
    message(STATUS "Found CLANG compiler")
    set(ROOT_ARCHITECTURE linuxs390gcc)
  else()
    message(FATAL_ERROR "There is no Setup for this compiler up to now. Don't know what to do. Stop cmake at this point.")
  endif()
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES i686)
  message(STATUS "Found a 32bit system")
  set(FP_MATH_FLAGS "-msse -mfpmath=sse")
  if(CMAKE_COMPILER_IS_GNUCXX)
    message(STATUS "Found GNU compiler collection")
    set(ROOT_ARCHITECTURE linux)
  elseif(CMAKE_CXX_COMPILER_ID STREQUAL Clang)
    message(STATUS "Found CLANG compiler")
    set(ROOT_ARCHITECTURE linux)
  elseif(CMAKE_CXX_COMPILER_ID STREQUAL Intel)
    set(ROOT_ARCHITECTURE linuxicc)
  else()
    message(FATAL_ERROR "There is no Setup for this compiler up to now. Don't know what to do. Stop cmake at this point.")
  endif()
else()
  message(FATAL_ERROR "Unknown processor: ${CMAKE_SYSTEM_PROCESSOR}")
endif()

# JIT must be able to resolve symbols from all ROOT binaries.
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -rdynamic")

if(CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pipe ${FP_MATH_FLAGS} -Wshadow -Wall -W -Woverloaded-virtual -fsigned-char")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pipe -Wall -W")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=legacy")

  set(CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS}")
  set(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS}")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--no-undefined -Wl,--hash-style=\"both\"")

  # Select flags.
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g -DNDEBUG")
  set(CMAKE_CXX_FLAGS_RELEASE        "-O3 -DNDEBUG")
  set(CMAKE_CXX_FLAGS_OPTIMIZED      "-Ofast -DNDEBUG")
  set(CMAKE_CXX_FLAGS_DEBUG          "-g")
  set(CMAKE_CXX_FLAGS_DEBUGFULL      "-g3")
  set(CMAKE_CXX_FLAGS_PROFILE        "-g3 -ftest-coverage -fprofile-arcs")
  set(CMAKE_C_FLAGS_RELWITHDEBINFO   "-O3 -g -DNDEBUG")
  set(CMAKE_C_FLAGS_RELEASE          "-O3 -DNDEBUG")
  set(CMAKE_C_FLAGS_OPTIMIZED        "-Ofast -DNDEBUG")
  set(CMAKE_C_FLAGS_DEBUG            "-g")
  set(CMAKE_C_FLAGS_DEBUGFULL        "-g3 -fno-inline")
  set(CMAKE_C_FLAGS_PROFILE          "-g3 -fno-inline -ftest-coverage -fprofile-arcs")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL Clang)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pipe ${FP_MATH_FLAGS} -Wall -W -Woverloaded-virtual -fsigned-char")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pipe -Wall -W")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=legacy")

  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 8)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wshadow")
  endif()

  set(CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS}")
  set(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS}")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--no-undefined")

  # Select flags.
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g -DNDEBUG")
  set(CMAKE_CXX_FLAGS_RELEASE        "-O2 -DNDEBUG")
  set(CMAKE_CXX_FLAGS_OPTIMIZED      "-O3 -ffast-math -DNDEBUG")
  set(CMAKE_CXX_FLAGS_DEBUG          "-g")
  set(CMAKE_CXX_FLAGS_DEBUGFULL      "-g3")
  set(CMAKE_CXX_FLAGS_PROFILE        "-g3 -ftest-coverage -fprofile-arcs")
  set(CMAKE_C_FLAGS_RELWITHDEBINFO   "-O2 -g -DNDEBUG")
  set(CMAKE_C_FLAGS_RELEASE          "-O2 -DNDEBUG")
  set(CMAKE_C_FLAGS_OPTIMIZED        "-O3 -ffast-math -DNDEBUG")
  set(CMAKE_C_FLAGS_DEBUG            "-g")
  set(CMAKE_C_FLAGS_DEBUGFULL        "-g3")
  set(CMAKE_C_FLAGS_PROFILE          "-g3 -ftest-coverage -fprofile-arcs")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL Intel)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -wd1476")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -restrict")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")

  # Check icc compiler version and set compile flags according to the
  execute_process(COMMAND ${CMAKE_CXX_COMPILER} -v
                  ERROR_VARIABLE _icc_version_info ERROR_STRIP_TRAILING_WHITESPACE)

  string(REGEX REPLACE "(^V|^icc[ ]v|^icpc[ ]v)ersion[ ]([0-9]+)\\.[0-9]+.*" "\\2" ICC_MAJOR "${_icc_version_info}")
  string(REGEX REPLACE "(^V|^icc[ ]v|^icpc[ ]v)ersion[ ][0-9]+\\.([0-9]+).*" "\\2" ICC_MINOR "${_icc_version_info}")

  message(STATUS "Found ICC major version ${ICC_MAJOR}")
  message(STATUS "Found ICC minor version ${ICC_MINOR}")

  if(ICC_MAJOR GREATER 9 OR ICC_MAJOR EQUAL 9)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -wd1572")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -wd1572")
  endif()

  if(ICC_MAJOR GREATER 11 OR ICC_MAJOR EQUAL 11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -wd279")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -wd279")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--no-undefined")
  endif()

  if(ICC_MAJOR GREATER 14 OR ICC_MAJOR EQUAL 14)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -wd873 -wd2536")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -wd873 -wd2536")
  endif()

  if(ICC_MAJOR GREATER 15 OR ICC_MAJOR EQUAL 15)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -wd597 -wd1098 -wd1292 -wd1478 -wd3373")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -wd597 -wd1098 -wd1292 -wd1478 -wd3373")
  endif()

  set(CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS}")
  set(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS}")

  # Select flags.
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -fp-model precise -g -DNDEBUG")
  set(CMAKE_CXX_FLAGS_RELEASE        "-O2 -fp-model precise -DNDEBUG")
  set(CMAKE_CXX_FLAGS_OPTIMIZED      "-O3 -DNDEBUG")
  set(CMAKE_CXX_FLAGS_DEBUG          "-O0 -g")
  set(CMAKE_C_FLAGS_RELWITHDEBINFO   "-O2 -fp-model precise -g -DNDEBUG")
  set(CMAKE_C_FLAGS_RELEASE          "-O2 -fp-model precise -DNDEBUG")
  set(CMAKE_C_FLAGS_OPTIMIZED        "-O3 -DNDEBUG")
  set(CMAKE_C_FLAGS_DEBUG            "-O0 -g")
endif()
