cmake_minimum_required(VERSION 3.0.2)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_EXTENSIONS YES)

project(celllist VERSION 0.1.0 LANGUAGES CXX)

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

# install target
include(GNUInstallDirs)
install(TARGETS celllist
    EXPORT celllist-export
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION include
    )
install(FILES cmake/celllistConfig.cmake DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/celllist)

# export celllist
install(EXPORT celllist-export
    FILE
    celllistTargets.cmake
    NAMESPACE
    celllist::
    DESTINATION
    ${CMAKE_INSTALL_LIBDIR}/cmake/celllist
    )

include_directories(include/)
add_executable(test src/test.cpp)
target_link_libraries(test celllist)


