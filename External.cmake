##################################################################################
# INSTALLED LIBRARIES
##################################################################################
#find_package(PNG)
#if(PNG_FOUND)
#	include_directories(${PNG_INCLUDE_DIR})
#	list(APPEND svg_cpp_plot_libs ${PNG_LIBRARY})
#	list(APPEND svg_cpp_plot_defs "USE_PNG")
#endif(PNG_FOUND)

######################################################################
# EXTERNAL LIBRARIES (GitHub etc...
######################################################################
if (NOT EXTERNAL_INSTALL_LOCATION)
	set(EXTERNAL_INSTALL_LOCATION ${CMAKE_CURRENT_SOURCE_DIR}/external)
endif()
if (NOT IS_DIRECTORY ${EXTERNAL_INSTALL_LOCATION})
	file(MAKE_DIRECTORY ${EXTERNAL_INSTALL_LOCATION})
endif()

include(ExternalProject)
# External include directory
include_directories(${EXTERNAL_INSTALL_LOCATION})
add_custom_target(update)

#ExternalProject_Add(svg-cpp-plot
#  GIT_REPOSITORY https://github.com/adolfomunoz/svg-cpp-plot.git
#  SOURCE_DIR ${EXTERNAL_INSTALL_LOCATION}/svg-cpp-plot
#  UPDATE_DISCONNECTED 1
#  STEP_TARGETS update
#  BUILD_COMMAND ""
#  CONFIGURE_COMMAND ""
#  INSTALL_COMMAND ""
#)
#add_dependencies(update svg-cpp-plot-update)






