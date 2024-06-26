project(dechat)

cmake_minimum_required(VERSION 3.20)
set(CMAKE_CXX_STANDARD 17)
#...............................................................................................................#

SET (GATB_CORE_HOME  ${PROJECT_SOURCE_DIR}/gatb-core)

add_compile_options(-mfma)

set(HIFIASM_PATH ${PROJECT_SOURCE_DIR}/aligner-correct)
message(${PROJECT_SOURCE_DIR}/aligner-correct)

aux_source_directory(${HIFIASM_PATH} HIFIASM_SRC_LIST)
# include_directories(${GraphAlignerTest_include})
include_directories(${HIFIASM_PATH})

find_package(Threads REQUIRED)



# find_package(Protobuf REQUIRED)

find_package(ZLIB REQUIRED)



# include_directories(${Protobuf_INCLUDE_DIRS})

include_directories(${ZLIB_INCLUDE_DIRS})
#...............................................................................................................#
# Define cmake modules directory
FOREACH (path "${CMAKE_CURRENT_SOURCE_DIR}/cmake" "${GATB_CORE_HOME}/gatb-core/cmake")
    IF (EXISTS "${path}")
        SET (CMAKE_MODULE_PATH  "${CMAKE_MODULE_PATH}" "${path}")
    ENDIF()
ENDFOREACH(path)

#...............................................................................................................#


# We don't want to install some GATB-CORE artifacts
SET (GATB_CORE_EXCLUDE_TOOLS     1) # no need to compile dbgh5, etc..
SET (GATB_CORE_EXCLUDE_TESTS     1)
SET (GATB_CORE_EXCLUDE_EXAMPLES  1)


include (GatbCore)
# We add the compilation options for the library
# g++11
add_definitions (${gatb-core-flags})
# We add the gatb-core include directory
include_directories (${gatb-core-includes})
# We add the path for extra libraries
link_directories (${gatb-core-extra-libraries-path})

find_package(Boost COMPONENTS program_options REQUIRED)
# include_directories(${Boost_INCLUDE_DIRS})

#...............................................................................................................#
#Debug
# set(CMAKE_BUILD_TYPE Debug)


#...............................................................................................................#

set (program "dechat")


set (PROGRAM_SOURCE_DIR ${PROJECT_SOURCE_DIR}/src)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)


include_directories(${PROGRAM_SOURCE_DIR})


file (GLOB_RECURSE  ProjectFiles  ${PROGRAM_SOURCE_DIR}/*.cpp)
# aux_source_directory(./src SRC_LIST)

add_executable(${program} ${ProjectFiles} ${HIFIASM_SRC_LIST})


# add_definitions(-DDEBUG2)


target_link_libraries(${program} ${gatb-core-libraries})
# target_link_libraries(demo PRIVATE  -lpthread -lm)
target_link_libraries(${program} 
    boost_program_options 
    # ${PROTOBUF_LIBRARIES}
    ${ZLIB_LIBRARIES}
)


target_link_libraries(${program} ZLIB::ZLIB)

target_link_libraries(${program} Threads::Threads)
# target_link_libraries(${program} ${GraphAligner_PATH}/MEMfinder/lib/memfinder.a)
target_link_libraries(${program} -lpthread -lm)
#...............................................................................................................#
#
# We set the version number
SET (CPACK_PACKAGE_VERSION              "")  

# We have to tell what is the server name
SET (CPACK_GFORGE_PROJECT_NAME          "gatb-tools")

# We set the kind of archive
SET (CPACK_GENERATOR                    "TGZ")
SET (CPACK_SOURCE_GENERATOR             "TGZ")

# We ignore unwanted files for the source archive
SET (CPACK_SOURCE_IGNORE_FILES
    "^${PROJECT_SOURCE_DIR}/.git/"     ;
    "^${PROJECT_SOURCE_DIR}/.gitmodules" ;
    "^${PROJECT_SOURCE_DIR}/.gitignore" ;
    "^${PROJECT_SOURCE_DIR}/build/"  ;
    "^${GATB_CORE_HOME}/.cproject" ;
    "^${GATB_CORE_HOME}/.git/"     ;
    "^${GATB_CORE_HOME}/.project"  ;
    "^${GATB_CORE_HOME}/.gitignore"
)

# We copy the project binary to the 'bin' directory
INSTALL (TARGETS ${PROJECT_NAME} DESTINATION bin)
INSTALL (DIRECTORY "${PROJECT_SOURCE_DIR}/example/" DESTINATION example)
INSTALL (FILES VERSION LICENSE README.md DESTINATION bin/..)

# cmake_system_name for mac is Darwin, i want to change that
if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    set(PRETTY_SYSTEM_NAME "Linux")
elseif (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(PRETTY_SYSTEM_NAME "Mac")
endif()

set (CPACK_PACKAGE_FILE_NAME  ${PROJECT_NAME}-binaries-${VERSION}-${PRETTY_SYSTEM_NAME})

include (CPack)