add_subdirectory( boundary_conditions )
# Propagate lib_srcs list to parent
add_lib_sources(
    soln_type.f90
    grid_type.f90
    bc_type.f90
    point_type.f90
)

