cmake_minimum_required(VERSION 3.20)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "-std=c++14 -mfma")
project(test)

aux_source_directory(${PROJECT_SOURCE_DIR} SRC_LIST)

add_executable(demo ${SRC_LIST})


target_include_directories(demo PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(demo PRIVATE  -lpthread -lm)


find_package(ZLIB REQUIRED)


if(ZLIB_FOUND)
    target_link_libraries(demo PRIVATE ZLIB::ZLIB)
endif()
