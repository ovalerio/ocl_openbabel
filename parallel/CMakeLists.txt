set(applications
#  obconformersearch_default
#  obconversion_readstring
#  obforcefield_energy
  obenergy_mmff94
  obminimize_mmff94
  obconformer_mmff94
  results_mmff94
  validate_mmff94
  obminimize_gpu
)

foreach(application ${applications})
  add_executable(${application} ${application}.cpp)
  target_link_libraries(${application} ${BABEL_LIBRARY})
  install(TARGETS ${application}
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
  )
endforeach(application)
