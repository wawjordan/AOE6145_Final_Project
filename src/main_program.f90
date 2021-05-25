program main_program
  
  use set_constants, only : zero, set_derived_constants
  use fluid_constants, only : set_fluid_constants
  use set_inputs, only : set_derived_inputs, geometry_file
  use set_inputs, only : i_low,i_high,j_low,j_high, n_ghost, cons
  use file_handling, only : grid_in, grid_out
  use geometry, only : setup_geometry, teardown_geometry
  use other_subroutines
  use time_integration
  use limiter_calc, only : select_limiter, calc_consecutive_variations, limiter_fun
  !use init_problem, only : initialize
  !use namelist, only : read_namelist
  use bc_type, only : subsonic_mms_bc
  use grid_type, only : grid_t
  use soln_type, only : soln_t, calc_mms
  use flux_calc, only : select_flux, calc_flux_2D
  use quadrature
  implicit none
  
  !character(len=100) :: header_str1
  !character(len=100) :: niter
  integer :: i, j, k
  type( grid_t )      :: grid
  type( soln_t )      :: soln
  type( subsonic_mms_bc ) :: bndry1
  integer, dimension(:), allocatable :: ind1, ind2
  
  call set_derived_constants
  call set_fluid_constants
  call select_flux()
  call select_limiter()
  !call read_namelist('input.nml')
  call set_derived_inputs
  call setup_geometry(grid,soln)
  
  call calc_mms(grid,soln)
  call output_exact_soln(grid,soln)
  call output_file_headers
  soln%V = soln%Vmms
  !soln%V(i_low:i_high,j_low:j_high,:) = soln%Vmms(i_low:i_high,j_low:j_high,:)
  soln%S = soln%Smms
  
  !call Limit(soln%V,soln%psi_plus,soln%psi_minus)
  
  call output_soln(grid,soln,0)
  do i = 1,140
  
  call prim2cons(soln%U,soln%V)
  !call Limit(soln%V,soln%psi_plus,soln%psi_minus)
  !soln%S = soln%Smms
  call calc_flux_2D(soln,grid,soln%F)
  !call calc_sources(soln,grid)
  call calc_time_step(grid%A_xi,grid%A_eta,grid%n_xi_avg, &
                      grid%n_eta_avg,grid%V,soln%V,soln%dt)
  call explicit_RK(grid,soln%S,soln%dt,soln%F,soln%U,soln%R,1)
  call update_states(soln)
  call calc_de(soln,soln%DE,soln%DEnorm,0,cons)
  !if (mod(i,10)==0) then
  call output_soln(grid,soln,i)
  !end if
  end do
  call grid_out(geometry_file,grid)
  call teardown_geometry(grid,soln)
  close(50)
  write(*,*) 'Program End'
end program main_program
