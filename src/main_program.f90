program main_program
  
  use set_constants, only : set_derived_constants
  use fluid_constants, only : set_fluid_constants
  use set_inputs, only : set_derived_inputs
  use file_handling, only : grid_in, grid_out
  use geometry, only : setup_geometry, teardown_geometry
  !use init_problem, only : initialize
  !use namelist, only : read_namelist
  use grid_type, only : grid_t
  use soln_type, only : soln_t
  implicit none
  
  !character(len=100) :: header_str1
  !character(len=100) :: niter
  !integer :: j, pnorm
  type( grid_t )      :: grid
  type( soln_t )      :: soln
  
  call set_derived_constants
  call set_fluid_constants
  !call read_namelist('input.nml')
  call set_derived_inputs
  call setup_geometry(grid,soln)
  
  
  call teardown_geometry(grid,soln)
  close(50)
end program main_program
