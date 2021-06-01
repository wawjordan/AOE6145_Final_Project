program main_program
  
  use set_constants, only : zero, set_derived_constants
  use fluid_constants, only : set_fluid_constants
  use set_inputs, only : set_derived_inputs, geometry_file
  use set_inputs, only : i_low,i_high,j_low,j_high, n_ghost, cons
  use set_inputs, only : ig_low,ig_high,jg_low,jg_high
  use file_handling, only : grid_in, grid_out
  use geometry, only : setup_geometry, teardown_geometry
  use variable_conversion, only : prim2cons, cons2prim
  use other_subroutines, only : output_exact_soln, output_file_headers,&
                                output_soln, calc_de
  use time_integration, only : calc_time_step, explicit_RK
  use limiter_calc, only : select_limiter, calc_consecutive_variations, limiter_fun
  !use init_problem, only : initialize
  !use namelist, only : read_namelist
  use bc_type, only : dirichlet_mms_bc_t
  use grid_type, only : grid_t
  use soln_type, only : soln_t, calc_mms
  use flux_calc, only : select_flux, calc_flux_2D
  implicit none
  
  !character(len=100) :: header_str1
  !character(len=100) :: niter
  integer :: i, j, k
  type( grid_t )      :: grid
  type( soln_t )      :: soln
  !type( dirichlet_mms_bc ) :: inlet_left, inlet_bottom
  !type( outflow_bc ) 
  !integer, dimension(:), allocatable :: ind1, ind2
  
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
  soln%U = soln%Umms
  !soln%V(i_low:i_high,j_low:j_high,:) = soln%Vmms(i_low:i_high,j_low:j_high,:)

  soln%S = soln%Smms
  
  !call Limit(soln%V,soln%psi_plus,soln%psi_minus)
  
  !call update_states(soln)
  call output_soln(grid,soln,0)
  do i = 1,250


  call prim2cons(soln%U,soln%V)
  !soln%U(ig_low:i_low-1,jg_low:jg_high,:) = soln%Umms(ig_low:i_low-1,jg_low:jg_high,:)
  !soln%U(i_high+1:ig_high,jg_low:jg_high,:) = soln%Umms(i_high+1:ig_high,jg_low:jg_high,:)
  !soln%U(ig_low:ig_high,jg_low:j_low-1,:) = soln%Umms(ig_low:ig_high,jg_low:j_low-1,:)
  !soln%U(ig_low:ig_high,j_high+1:jg_high,:) = soln%Umms(ig_low:ig_high,j_high+1:jg_high,:)
  write(*,*) i
  !call Limit(soln%V,soln%psi_plus,soln%psi_minus)
  !soln%S = soln%Smms
  call calc_flux_2D(soln,grid,soln%F)
  !call calc_sources(soln,grid)
  call calc_time_step(grid%A_xi,grid%A_eta,grid%n_xi_avg, &
                      grid%n_eta_avg,grid%V,soln%V,soln%dt)
  call explicit_RK(grid,soln%S,soln%dt,soln%F,soln%U,soln%R,1)
  call cons2prim(soln%U,soln%V)
  soln%V(ig_low:i_low-1,jg_low:jg_high,:) = soln%Vmms(ig_low:i_low-1,jg_low:jg_high,:)
  !soln%V(i_high+1:ig_high,jg_low:jg_high,:) = soln%Vmms(i_high+1:ig_high,jg_low:jg_high,:)
  soln%V(ig_low:ig_high,jg_low:j_low-1,:) = soln%Vmms(ig_low:ig_high,jg_low:j_low-1,:)
  !soln%V(ig_low:ig_high,j_high+1:jg_high,:) = soln%Vmms(ig_low:ig_high,j_high+1:jg_high,:)
  !call update_states(soln)
  !call calc_de(soln,soln%DE,soln%DEnorm,0,cons)
  if (mod(i,10)==0) then
  call output_soln(grid,soln,i)
  end if
  end do
  call grid_out(geometry_file,grid)
  call teardown_geometry(grid,soln)
  close(50)
  write(*,*) 'Program End'
end program main_program
