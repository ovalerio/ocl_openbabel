CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
project(part1)

#need to include FindOPENCL.cmake to
SET( CMAKE_MODULE_PATH ${part1_SOURCE_DIR}/../cmake/modules)

message("CMake module path: ${CMAKE_MODULE_PATH}\n")

#important to use the right capitalization for the module to be found
FIND_PACKAGE (OpenCL)

INCLUDE_DIRECTORIES( 
#    ${part1_SOURCE_DIR}/../include/opencl10
    ${part1_SOURCE_DIR}/../include/CL
    ${part1_SOURCE_DIR}/../include/viennacl
)

message("include directories: \n")
get_directory_property(includes INCLUDE_DIRECTORIES)
message("${includes}\n")

#provide the source path so we can open our .cl file at runtime
ADD_DEFINITIONS(-DCL_SOURCE_DIR="${part1_SOURCE_DIR}")
#build in debug mode
ADD_DEFINITIONS(-g )

#the source files we want to compile into the library
set (LIBCXXFILES cll.cpp part1.cpp util.cpp)
#set a CMake variable to name the library
SET(cllib part1)
ADD_LIBRARY(${cllib} ${LIBCXXFILES})


#create the executable
SET(EXEC part1.exe)
ADD_EXECUTABLE(${EXEC} main.cpp)

TARGET_LINK_LIBRARIES (${EXEC}
   ${cllib}
   ${OPENCL_LIBRARIES}
)

  install(TARGETS ${EXEC}
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
  )


#create the executable for obenergy applications
set(applications
  obabelx
  obenergyx
  obconformerx
  obminimizex
  viennacl-info
  opencl
  vector
#  eigen-with-viennacl # not working with Eigen 3  needs to be ported
)

foreach(application ${applications})
  add_executable(${application} ${application}.cpp)
  target_link_libraries(${application} ${BABEL_LIBRARY} ${OPENCL_LIBRARIES})
  install(TARGETS ${application}
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
  )
endforeach(application)
  
  