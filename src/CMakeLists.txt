macro(add_project_executable exe source)
  add_executable(${exe} ${source})
  target_link_dune_default_libraries(${exe})
  dune_target_enable_all_packages(${exe})
  target_compile_definitions(${exe} PRIVATE DUNE_GRID_PATH=\"${CMAKE_CURRENT_LIST_DIR}/grids\")

  if (dune-common_VERSION VERSION_GREATER_EQUAL 2.11)
    dune_default_include_directories(${exe} PRIVATE)
  endif()
endmacro()


# add_project_executable(boxsurface boxsurface.cc)

# Experimental order of convergence for curve geometric flows
foreach(FLOW 1 2)
  foreach(SURFACE 1 2)
    set(exe "eoc_curve_flow${FLOW}_surf${SURFACE}")
    add_project_executable(${exe} eoc_curve.cc)
    target_compile_definitions(${exe} PRIVATE FLOW=${FLOW} SURFACE=${SURFACE})
  endforeach()
endforeach()

foreach(FLOW 1 2)
  foreach(SURFACE 1 2)
    set(exe "curve_flow${FLOW}_surf${SURFACE}")
    add_project_executable(${exe} curve.cc)
    target_compile_definitions(${exe} PRIVATE FLOW=${FLOW} SURFACE=${SURFACE})
  endforeach()
endforeach()

foreach(FLOW 2)
  foreach(SURFACE 1 2)
    set(exe "curve_sp_flow${FLOW}_surf${SURFACE}")
    add_project_executable(${exe} curve_sp.cc)
    target_compile_definitions(${exe} PRIVATE FLOW=${FLOW} SURFACE=${SURFACE})
  endforeach()
endforeach()



# Experimental order of convergence for surfaces geometric flows
foreach(FLOW 1 2)
  foreach(SURFACE 1 2)
    set(exe "eoc_surface_flow${FLOW}_surf${SURFACE}")
    add_project_executable(${exe} eoc_surface.cc)
    target_compile_definitions(${exe} PRIVATE FLOW=${FLOW} SURFACE=${SURFACE})
  endforeach()
endforeach()

foreach(FLOW 1 2)
  foreach(SURFACE 1 2)
    set(exe "surface_flow${FLOW}_surf${SURFACE}")
    add_project_executable(${exe} surface.cc)
    target_compile_definitions(${exe} PRIVATE FLOW=${FLOW} SURFACE=${SURFACE})
  endforeach()
endforeach()

foreach(FLOW 2)
  foreach(SURFACE 1 2)
    set(exe "surface_sp_flow${FLOW}_surf${SURFACE}")
    add_project_executable(${exe} surface_sp.cc)
    target_compile_definitions(${exe} PRIVATE FLOW=${FLOW} SURFACE=${SURFACE})
  endforeach()
endforeach()
