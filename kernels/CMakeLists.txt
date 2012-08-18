CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
project(kernels)

#
# export file: copy it to the build tree on every build invocation and add rule for installation
#
function    (cm_export_file FILE DEST)
  if    (NOT TARGET export-files)
    add_custom_target(export-files ALL COMMENT "Exporting files into build tree")
  endif (NOT TARGET export-files)
  get_filename_component(FILENAME "${FILE}" NAME)
  add_custom_command(TARGET export-files COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/${FILE}" "${CMAKE_CURRENT_BINARY_DIR}/${DEST}/${FILENAME}")
  install(FILES "${FILE}" DESTINATION "${DEST}")
endfunction (cm_export_file)


message("Moving the kernels *.cl files to build directory")

set(kernels
  simple_kernel.cl
  electrostatic_forces_kernel.cl
)

foreach(kernel ${kernels})
    cm_export_file(${kernel} "bin/kernels/")
endforeach(kernel)
  
#cm_export_file("electrostatic_forces_kernel.cl" "kernels/")
