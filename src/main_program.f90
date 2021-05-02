program main_program
  
  use set_constants, only : zero, set_derived_constants
  use fluid_constants, only : set_fluid_constants
  use set_inputs, only : set_derived_inputs, geometry_file
  use set_inputs, only : i_low,i_high,j_low,j_high, n_ghost
  use file_handling, only : grid_in, grid_out
  use geometry, only : setup_geometry, teardown_geometry
  use other_subroutines, only : output_file_headers, output_soln
  use time_integration
  use limiter_calc, only : select_limiter, calc_consecutive_variations, limiter_fun
  !use init_problem, only : initialize
  !use namelist, only : read_namelist
  use grid_type, only : grid_t
  use soln_type, only : soln_t, calc_mms
  use flux_calc, only : select_flux, calc_flux_2D
  implicit none
  
  !character(len=100) :: header_str1
  !character(len=100) :: niter
  integer :: j
  type( grid_t )      :: grid
  type( soln_t )      :: soln  
  real(prec), dimension(20,4) :: r_plus, r_minus, psi_plus, psi_minus
  r_plus = zero
  r_minus =  zero
  psi_plus = zero
  psi_minus =  zero
  call set_derived_constants
  call set_fluid_constants
  call select_flux()
  call select_limiter()
  !call read_namelist('input.nml')
  call set_derived_inputs
  call setup_geometry(grid,soln)
  call calc_mms(grid,soln)
  soln%V(i_low:i_high,j_low:j_high,:) = soln%Vmms
  write(*,*)
  write(*,*) lbound(soln%V,1), ubound(soln%V,1)
  write(*,*) lbound(soln%V,2), ubound(soln%V,2)
  write(*,*)
  do j = j_low,j_high
    write(*,*)'V',j,soln%V(i_low,j,:)
  end do
  write(*,*)
  call calc_flux_2D(soln%V,grid%n_xi,grid%n_eta,soln%F)
  !write(*,*)'Vmms', soln%Vmms(i_low,j_low,:)
  !write(*,*)'V', soln%V(i_low,j_low,:)
  !call output_file_headers
  !call output_soln(grid,soln,1)
  
  call grid_out(geometry_file,grid)
  call teardown_geometry(grid,soln)
  close(50)
end program main_program
