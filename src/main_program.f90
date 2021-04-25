program main_program
  
  use set_constants, only : set_derived_constants
  use fluid_constants, only : set_fluid_constants
  use set_inputs, only : set_derived_inputs
  !use set_inputs, only : max_iter, tol, soln_save, res_save
  !use set_inputs, only : leftV, rightV, leftU, rightU, flux_scheme
  !use set_inputs, only : limiter_freeze, res_out, cons
  !use variable_conversion
  !use time_integration
  !use boundary_conds, only : enforce_bndry
  !use limiter_calc, only : select_limiter
  !use flux_calc, only : select_flux, flux_fun
  use other_subroutines
  !use geometry, only : setup_geometry, teardown_geometry
  !use init_problem, only : initialize
  use namelist, only : read_namelist
  use grid_type
  use soln_type
  implicit none
  
  character(len=100) :: header_str1
  character(len=100) :: niter
  integer :: j, pnorm
  type( grid_t )      :: grid
  type( soln_t )      :: soln
  
  !open(50,file='temp.txt',status='unknown')
  
  call set_derived_constants
  call set_fluid_constants
  !call read_namelist('input.nml')
  call set_derived_inputs
  call allocate_grid( grid )
  call deallocate_grid( grid )
  !call setup_geometry(grid,soln)
  
  
  !call teardown_geometry(grid,soln)
  close(50)
end program main_program
