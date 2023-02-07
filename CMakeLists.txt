cmake_minimum_required(VERSION 3.0.2)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_EXTENSIONS YES)

project(celllist)

set(objects src/celllist.cpp)
set(headers include/celllist.h)
set_source_files_properties(${objects} PROPERTIES LANGUAGE CXX)
set_source_files_properties(${headers} PROPERTIES HEADER_FILE_ONLY ON)

add_library(celllist ${objects} ${headers})
target_include_directories(celllist
    PUBLIC 
    $<INSTALL_INTERFACE:include>    
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    )
set_target_properties(celllist PROPERTIES PUBLIC_HEADER "${headers}")
target_compile_options(celllist PRIVATE -Wall -Wextra -Wshadow -Wconversion -Wno-language-extension-token)

include_directories(include/)
add_executable(test src/test.cpp)

