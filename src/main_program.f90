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
  implicit none
  
  !character(len=100) :: header_str1
  !character(len=100) :: niter
  integer :: i, j
  type( grid_t )      :: grid
  type( soln_t )      :: soln
  type( subsonic_mms_bc ) :: bndry1, bndry2, bndry3, bndry4
  !real(prec), dimension(20,4) :: r_plus, r_minus, psi_plus, psi_minus
  integer, dimension(:,:), allocatable :: ind1, ind2, ind3, ind4
  
  !r_plus = zero
  !r_minus =  zero
 ! psi_plus = zero
 ! psi_minus =  zero
  call set_derived_constants
  call set_fluid_constants
  call select_flux()
  call select_limiter()
  !call read_namelist('input.nml')
  call set_derived_inputs
  call setup_geometry(grid,soln)
  
  allocate(ind1(n_ghost*(i_high-i_low),2))
  allocate(ind2(n_ghost*(i_high-i_low+1),2))
  allocate(ind3(n_ghost*(j_high-j_low+1),2))
  allocate(ind4(n_ghost*(j_high-j_low+1),2))
write(*,*) ( ( i, i=i_low,i_high ), j=j_low-1,jg_low,-1)
write(*,*) ( ( j, i=i_low,i_high ), j=j_low-1,jg_low,-1)
write(*,*) lbound(ind1,1), ubound(ind1,1)
write(*,*) lbound(ind1,2), ubound(ind1,2)
do j = j_low-1,jg_low,-1
do i = i_low,i_high
ind1(i,1) = i
ind1(i,2) = j
write(*,*) ind1(i,:)
end do
end do
!ind1(:,1) = [( ( i, i=i_low,i_high ), j=j_low-1,jg_low,-1)]
!ind1(:,2) = [( ( j, i=i_low,i_high ), j=j_low-1,jg_low,-1)]
!ind1(:,1) = [( ( i, i=i_low,i_high ), j=j_low-1,j_low-1)]
!ind1(:,2) = [( ( j, i=i_low,i_high ), j=j_low-1,j_low-1)]
!  
!  ind2(:,1) = [( ( i, i=i_low,i_high ), j=j_high,jg_high)]
!  ind2(:,2) = [( ( j, i=i_low,i_high ), j=j_high,jg_high)]
!  
!  ind3(:,1) = [( ( i, i=i_low,ig_low,-1 ), j=j_low,j_high)]
!  ind3(:,2) = [( ( j, i=i_low,ig_low,-1 ), j=j_low,j_high)]
!  
!  ind4(:,1) = [( ( i, i=i_high,ig_high ), j=j_low,j_high)]
!  ind4(:,2) = [( ( j, i=i_high,ig_high ), j=j_low,j_high)]
  write(*,*) 'i_low= ',i_low
  write(*,*) 'i_high= ',i_high
  write(*,*) 'j_low= ',j_low
  write(*,*) 'j_high= ',j_high
  write(*,*)
!  do i = 1, size(ind1,1)
!  write(*,*)'ind1', ind1(i,1),ind1(i,2)
!  end do
  write(*,*)
!  do i = 1, size(ind2,1)
!  write(*,*)'ind2', ind2(i,1),ind2(i,2)
!  end do
!  write(*,*)
!  do i = 1, size(ind3,1)
!  write(*,*)'ind3', ind3(i,1),ind3(i,2)
!  end do
!  write(*,*)
!  do i = 1, size(ind4,1)
!  write(*,*)'ind4', ind4(i,1),ind4(i,2)
!  end do
  
  call bndry1%init_mms_bc(ind1,grid) 
  call bndry2%init_mms_bc(ind2,grid)
  call bndry3%init_mms_bc(ind3,grid)
  call bndry4%init_mms_bc(ind4,grid)
  
  call bndry1%enforce(soln)
  call bndry2%enforce(soln)
  call bndry3%enforce(soln)
  call bndry4%enforce(soln)
  
  call calc_mms(grid,soln)
  
  soln%V(i_low:i_high,j_low:j_high,:) = soln%Vmms
  
  !call calc_flux_2D(soln%V,grid%n_xi,grid%n_eta,soln%F)
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
  
  deallocate(ind1,ind2,ind3,ind4)
  call grid_out(geometry_file,grid)
  call teardown_geometry(grid,soln)
  close(50)
end program main_program
