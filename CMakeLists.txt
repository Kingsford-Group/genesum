cmake_minimum_required (VERSION 2.8)

project(genesum)
set(CPACK_PACKAGE_VERSION "0.5.0")
set(CPACK_PACKAGE_VERSION_MAJOR "0")
set(CPACK_PACKAGE_VERSION_MINOR "5")
set(CPACK_PACKAGE_VERSION_PATCH "0")
set(CPACK_GENERATOR "TGZ")
set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_PACKAGE_VENDOR "Carnegie Mellon University")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "genesum - gene-level summary of transcript-level expression estimates")
set(CPACK_PACKAGE_NAME
    "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
set(CPACK_SOURCE_PACKAGE_FILE_NAME
    "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}-Source")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

## Set the standard required compile flags
set (CMAKE_CXX_FLAGS "-g -funroll-loops -fPIC -O3 -D_FILE_OFFSET_BITS=64 -DHAVE_ANSI_TERM -DHAVE_SSTREAM -Wall -std=c++11 -Wno-unused-local-typedefs")

## OSX is strange (some might say, stupid in this regard). Deal with it's quirkines here.
if (APPLE)
    # To allow ourselves to build a dynamic library, we have to tell the compiler
    # that, yes, the symbols will be around at runtime.
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -undefined dynamic_lookup")
    set (GENESUM_LINKER_FLAGS "-all_load")
    # In order to "think different", we also have to use non-standard suffixes
    # for our shared libraries
    set(SHARED_LIB_EXTENSION "dylib")
    #set(HAVE_LOGGER FALSE)
else()
    # We're in sane linux world
    set (SHARED_LIB_EXTENSION "so")
    set (GENESUM_LINKER_FLAGS "")
endif()

##
# Compiler-specific C++11 activation.
# http://stackoverflow.com/questions/10984442/how-to-detect-c11-support-of-a-compiler-with-cmake
##
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
    execute_process(
        COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
    if (NOT (GCC_VERSION VERSION_GREATER 4.7 OR GCC_VERSION VERSION_EQUAL 4.7))
        message(FATAL_ERROR "${PROJECT_NAME} requires g++ 4.7 or greater.")
    endif ()

    if (NOT APPLE OR
            GCC_VERSION VERSION_GREATER 4.8.2 OR
            GCC_VERSION VERSION_EQUAL 4.8.2)
        #Perhaps the logger now works in both worlds
        set (CMAKE_CXX_FLAGS "-DHAVE_LOGGER ${CMAKE_CXX_FLAGS}")
        set (HAVE_LOGGER TRUE)
    else()
        message(WARNING "Logging requires gcc >= 4.8.2 on OSX; it will be disabled")
    endif()


elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    # On OSX, if we compile with clang, we want to use libc++ instead of libstdc++
    if (APPLE)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
    endif()
else ()
    message(FATAL_ERROR "Your C++ compiler does not support C++11.")
endif ()

##
# Record this top-level path
##
set (GAT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})

include(ExternalProject)
message("Build system will compile libgff")
message("==================================================================")
ExternalProject_Add(libgff
    DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
    DOWNLOAD_COMMAND curl -L https://github.com/Kingsford-Group/libgff/archive/v1.0.tar.gz -o libgff.tar.gz &&
                     tar -xzf libgff.tar.gz &&
                     rm -fr libgff &&
                     mv -f libgff-1.0 libgff
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/libgff
    UPDATE_COMMAND sh -c "mkdir -p <SOURCE_DIR>/build"
    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
    BINARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/libgff/build
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_SOURCE_DIR}/external/install
    )

###
#
# Done building external dependencies.
#
###

set (CPACK_SOURCE_IGNORE_FILES
    "/build/"
    "/bin/"
    "/lib/"
    "/external/"
    "/src/obsolete/"
    "/include/obsolete/"
    ".git/")

message("CPACK_SOURCE_IGNORE_FILES = ${CPACK_SOURCE_IGNORE_FILES}")

# Recurse into Sailfish source directory
add_subdirectory ( src )

# build a CPack driven installer package
include (CPack)
