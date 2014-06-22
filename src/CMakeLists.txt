set ( GENESUM_SRCS
    GeneSum.cpp
    )

include_directories(
    ${GAT_SOURCE_DIR}/include
    ${GAT_SOURCE_DIR}/external
    ${GAT_SOURCE_DIR}/external/install/include
    )

link_directories(
    ${GAT_SOURCE_DIR}/lib
    ${GAT_SOURCE_DIR}/external/install/lib
    )

# Build the genesum executable
add_executable(genesum ${GENESUM_SRCS})

# Link the executable
target_link_libraries(genesum
    pthread
    m
    gff
    )

IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    SET(CMAKE_INSTALL_PREFIX
        "${GAT_SOURCE_DIR}" CACHE PATH "Default install prefix" FORCE
        )
ENDIF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

set(INSTALL_LIB_DIR lib )
set(INSTALL_BIN_DIR bin )
set(INSTALL_INCLUDE_DIR include )

install(DIRECTORY
    ${GAT_SOURCE_DIR}/external/install/lib/
    DESTINATION ${INSTALL_LIB_DIR}
    )

install(TARGETS genesum
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    )

#set(POST_INSTALL_SCRIPT ${GAT_SOURCE_DIR}/cmake/PostInstall.cmake)

install(
    CODE
    "
    execute_process(COMMAND \"${CMAKE_COMMAND}\"
        -DCMAKE_SYSTEM_NAME=${CMAKE_SYSTEM_NAME}
        -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
        -P \"${POST_INSTALL_SCRIPT}\")
    "
    )

include(InstallRequiredSystemLibraries)
add_test( NAME simple_test COMMAND ${CMAKE_COMMAND} -DTOPLEVEL_DIR=${GAT_SOURCE_DIR} -P ${GAT_SOURCE_DIR}/cmake/SimpleTest.cmake )



