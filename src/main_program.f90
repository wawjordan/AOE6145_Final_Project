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
  allocate(ind1(n_ghost*(i_high-i_low+1)))
  allocate(ind2(n_ghost))


do j = jg_low,j_low-1
ind2(j-jg_low+1) = j
do i = i_low,i_high
k = (j-jg_low)*(i_high-i_low+1)+i
ind1(k) = i
end do
end do
  
  call bndry1%init_mms_bc(ind1,ind2,grid) 
  call bndry1%enforce(soln)
  
  call calc_mms(grid,soln)
  
  !soln%V(i_low:i_high,j_low:j_high,:) = soln%Vmms(i_low:i_high,j_low:j_high,:)
  soln%V = soln%Vmms
! stop 
  do i = 1,20
  call prim2cons(soln%U,soln%V)
  call calc_flux_2D(soln,grid,soln%F)
  call calc_time_step(grid%A_xi,grid%A_eta,grid%n_xi_avg, &
                      grid%n_eta_avg,grid%V,soln%V,soln%dt)
  call bndry1%enforce(soln)
  call explicit_RK(grid,soln%S,soln%dt,soln%F,soln%U,soln%R,4)
  call update_states(soln)
  end do
  !write(*,*)'Vmms', soln%Vmms(i_low,j_low,:)
  !write(*,*) 'i', (i,i=ig_low,ig_high)
  !do j = jg_low, jg_high
  !write(*,*)'j', j
  !end do
  !write(*,*)
  !do j = jg_low, jg_high
  !write(*,*)'V', soln%V(:,j,1)
  !end do
  call output_file_headers
  call output_soln(grid,soln,1)
  
  deallocate(ind1,ind2)
  !call grid_out(geometry_file,grid)
  call teardown_geometry(grid,soln)
  close(50)
end program main_program
