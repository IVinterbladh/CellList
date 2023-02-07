get_filename_component(celllist_CMAKE_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
include(CMakeFindDependencyMacro)

if(NOT TARGET celllist::celllist)
    include("${celllist_CMAKE_DIR}/celllistTargets.cmake")
endif()
