program main_program
  
  use set_constants, only : set_derived_constants
  use fluid_constants, only : set_fluid_constants
  use set_inputs, only : set_derived_inputs, geometry_file
  use set_inputs, only : i_low,i_high,j_low,j_high
  use file_handling, only : grid_in, grid_out
  use geometry, only : setup_geometry, teardown_geometry
  use other_subroutines, only : output_file_headers, output_soln
  !use init_problem, only : initialize
  !use namelist, only : read_namelist
  use grid_type, only : grid_t
  use soln_type, only : soln_t, calc_mms
  use flux_calc
  implicit none
  
  !character(len=100) :: header_str1
  !character(len=100) :: niter
  !integer :: j, pnorm
  type( grid_t )      :: grid
  type( soln_t )      :: soln  

  call set_derived_constants
  call set_fluid_constants
  call select_flux()
  !call read_namelist('input.nml')
  call set_derived_inputs
  call setup_geometry(grid,soln)
  call calc_mms(grid,soln)
  call output_file_headers
  call output_soln(grid,soln,1)
  !! flux at i-1/2 & i+1/2
!  call calc_flux_2D(soln%U(i_low-1:i_high,j_low:j_high+1,:),&
!                    soln%U(i_low:i_high+1,j_low:j_high+1,:),&
!                    grid%n_xi_x(i_low:i_high+1,j_low:j_high+1),&
!                    grid%n_xi_y(i_low:i_high+1,j_low:j_high+1),&
!                    soln%F(i_low:i_high+1,j_low:j_high+1,:))
!  
!  grid%n_xi_x(i_low:i_high,j_low:j_high) = &
!         soln%F(i_low+1:i_high+1,j_low:j_high,1)*&
!         grid%A_xi(i_low+1:i_high+1,j_low:j_high)&
!        -soln%F(i_low:i_high,j_low:j_high,1)*&
!         grid%A_xi(i_low:i_high,j_low:j_high)
!  
!  call calc_flux_2D(soln%U(i_low:i_high+1,j_low-1:j_high,:),&
!                    soln%U(i_low:i_high+1,j_low:j_high+1,:),&
!                    grid%n_eta_x(i_low:i_high+1,j_low:j_high+1),&
!                    grid%n_eta_y(i_low:i_high+1,j_low:j_high+1),&
!                    soln%F(i_low:i_high+1,j_low:j_high+1,:))
!  
!  grid%n_xi_y(i_low:i_high,j_low:j_high) = &
!         soln%F(i_low:i_high,j_low+1:j_high+1,1)*&
!         grid%A_eta(i_low:i_high,j_low+1:j_high+1)&
!        -soln%F(i_low:i_high,j_low:j_high,1)*&
!         grid%A_eta(i_low:i_high,j_low:j_high)
  
  call grid_out(geometry_file,grid)
  call teardown_geometry(grid,soln)
  close(50)
end program main_program
