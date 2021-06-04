program main_program
  
  use set_precision, only : prec  
  use grid_type, only : grid_t
  use soln_type, only : soln_t, calc_mms
  use point_type!, only : point_t, array_to_point
  use set_constants, only : zero, one, set_derived_constants
  implicit none
  
  type(point_t) :: point1
  
  point1 = (/one, one/)
  
  
   
!  use fluid_constants, only : set_fluid_constants
!  use set_inputs, only : set_derived_inputs, geometry_file
!  use set_inputs, only : i_low,i_high,j_low,j_high, n_ghost, cons, neq
!  use set_inputs, only : ig_low,ig_high,jg_low,jg_high
!  use set_inputs, only : res_save, res_out, soln_save, tol, Lmms
!  use file_handling, only : grid_in, grid_out
!  use geometry, only : setup_geometry, teardown_geometry
!  use variable_conversion, only : prim2cons, cons2prim, update_states, &
!                                  limit_primitives
!  use other_subroutines, only : output_exact_soln, output_file_headers,&
!                                output_soln, output_res, calc_de, Limit, MUSCL_extrap
!  use time_integration, only : calc_time_step, explicit_RK, residual_norms
!  use limiter_calc, only : select_limiter, calc_consecutive_variations, limiter_fun
!  use mms_functions, only : rho_mms, uvel_mms, vvel_mms, press_mms
!  !use init_problem, only : initialize
!  !use namelist, only : read_namelist
!  use bc_type, only : dirichlet_mms_bc_t
!  use grid_type, only : grid_t
!  use soln_type, only : soln_t, calc_mms
!  use flux_calc, only : select_flux, calc_flux_2D, flux_fun
!  implicit none
!  
!  integer :: i, j, k
!  type( grid_t )      :: grid
!  type( soln_t )      :: soln
!  real(prec), dimension(4) :: left, right
!  real(prec) :: nx, ny
!  integer :: i1,j1,k1
!  integer, dimension(4) :: ind
!  real(prec), dimension(4,4) :: Vtmp, psiPtmp, psiMtmp
!
!
!  
!  call set_derived_constants
!  call set_fluid_constants
!  call select_flux()
!  call select_limiter()
!  call set_derived_inputs
!  call setup_geometry(grid,soln)
!  
!  call calc_mms(grid,soln)
!  call output_exact_soln(grid,soln)
!  call output_file_headers
!  soln%V = soln%Vmms
!  
!  call prim2cons(soln%Umms,soln%Vmms)
!  soln%U = soln%Umms
!  soln%S = soln%Smms
!  
!  call output_soln(grid,soln,0)
!  do k = 1,100000
!!==============================================================================
!  soln%V(i_low:i_high,j_low-1,1) = &
!   2*rho_mms(Lmms,grid%x(i_low:i_high,j_low), grid%y(i_low:i_high,j_low) ) &
!   - soln%V(i_low:i_high,j_low,1)
!  soln%V(i_low:i_high,j_low-1,2) = &
!   2*uvel_mms(Lmms,grid%x(i_low:i_high,j_low), grid%y(i_low:i_high,j_low) ) &
!   - soln%V(i_low:i_high,j_low,2)
!  soln%V(i_low:i_high,j_low-1,3) = &
!   2*vvel_mms(Lmms,grid%x(i_low:i_high,j_low), grid%y(i_low:i_high,j_low) ) &
!   - soln%V(i_low:i_high,j_low,3)
!  soln%V(i_low:i_high,j_low-1,4) = &
!   2*press_mms(Lmms,grid%x(i_low:i_high,j_low), grid%y(i_low:i_high,j_low) ) &
!   - soln%V(i_low:i_high,j_low,4)
!   
!  !soln%V(i_low:i_high,j_low-2,:) = soln%V(i_low:i_high,j_low-1,:)
!  !soln%V(i_low:i_high,j_low-2,:) = 2*soln%V(i_low:i_high,j_low-1,:) - &
!  !                                   soln%V(i_low:i_high,j_low,:)
!  
!  soln%V(i_low-1,j_low:j_high,1) = &
!   2*rho_mms(Lmms,grid%x(i_low,j_low:j_high), grid%y(i_low,j_low:j_high) ) &
!   - soln%V(i_low-1,j_low:j_high,1)
!  soln%V(i_low-1,j_low:j_high,2) = &
!   2*uvel_mms(Lmms,grid%x(i_low,j_low:j_high), grid%y(i_low,j_low:j_high) ) &
!   - soln%V(i_low-1,j_low:j_high,2)
!  soln%V(i_low-1,j_low:j_high,3) = &
!   2*vvel_mms(Lmms,grid%x(i_low,j_low:j_high), grid%y(i_low,j_low:j_high) ) &
!   - soln%V(i_low-1,j_low:j_high,3)
!  soln%V(i_low-1,j_low:j_high,4) = &
!   2*press_mms(Lmms,grid%x(i_low,j_low:j_high), grid%y(i_low,j_low:j_high) ) &
!   - soln%V(i_low-1,j_low:j_high,4)
!   
!  !soln%V(i_low-2,j_low:j_high,:) = soln%V(i_low-1,j_low:j_high,:)
!  !soln%V(i_low-2,j_low:j_high,:) = 2*soln%V(i_low-1,j_low:j_high,:) - &
!  !                                 soln%V(i_low,j_low:j_high,:)
!!==============================================================================
!  
!  !soln%V(i_low:i_high,jg_high-1,:) = 2*soln%V(i_low:i_high,j_high,:) - &
!  !                                       soln%V(i_low:i_high,j_high-1,:)
!  !soln%V(ig_high-1,j_low:j_high,:) = 2*soln%V(i_high,j_low:j_high,:) - &
!  !                                       soln%V(i_high-1,j_low:j_high,:)
!  !
!  !soln%V(i_low:i_high,jg_high,:) = 2*soln%V(i_low:i_high,jg_high-1,:) - &
!  !                                       soln%V(i_low:i_high,j_high,:)
!  !soln%V(ig_high,j_low:j_high,:) = 2*soln%V(ig_high-1,j_low:j_high,:) - &
!  !                                       soln%V(i_high,j_low:j_high,:)
! 
! 
!  !soln%V(i_high+1,j_high+1,:) = 2*soln%V(i_high,j_high,:) - &
!  !                                       soln%V(i_high-1,j_high-1,:)
!  !soln%V(i_high+2,j_high+2,:) = 2*soln%V(i_high+1,j_high+1,:) - &
!  !                                       soln%V(i_high,j_high,:)
!  !
!  !soln%V(i_high+1,j_low-1,:) = 2*soln%V(i_high,j_low,:) - &
!  !                                       soln%V(i_high-1,j_low+1,:)
!  !soln%V(i_high+2,j_low-2,:) = 2*soln%V(i_high+1,j_low-1,:) - &
!  !                                       soln%V(i_high,j_low,:)
!  call limit_primitives(soln%V)
!  call prim2cons(soln%U,soln%V)
!  call update_states(soln)
!   
!  call Limit(soln%V,soln%psi_plus,soln%psi_minus)
!  call calc_flux_2D(soln,grid,soln%F)
!  
!  i1 = i_high+1
!  do j1 = j_low,j_high
!    nx = grid%n_xi(i1,j1,1)
!    ny = grid%n_xi(i1,j1,2)
!    !ind = (/ ( k1,k1=i1-2,i1+1 ) /)
!    !Vtmp(1,:) = soln%V(ind(1),j1,:)
!    !Vtmp(2,:) = soln%Vmms(ind(2),j1,:)
!    !Vtmp(3,:) = soln%Vmms(ind(3),j1,:)
!    !Vtmp(4,:) = soln%Vmms(ind(4),j1,:)
!    !psiPtmp = soln%psi_plus(ind,j1,:,1)
!    !psiMtmp = soln%psi_minus(ind,j1,:,1)
!    !call MUSCL_extrap(Vtmp, psiPtmp, psiMtmp, left, right)
!    left = soln%V(i1-1,j1,:)
!    right = soln%Vmms(i1,j1,:)
!    call flux_fun(left,right,nx,ny,soln%F(i1,j1,:,1))
!  end do
!  j1 = j_high+1
!  do i1 = i_low,i_high
!    nx = grid%n_eta(i1,j1,1)
!    ny = grid%n_eta(i1,j1,2)
!    !ind = (/ ( k1,k1=j1-2,j1+1 ) /)
!    !Vtmp(1,:) = soln%V(i1,ind(1),:)
!    !Vtmp(2,:) = soln%V(i1,ind(2),:)
!    !Vtmp(3,:) = soln%V(i1,ind(3),:)
!    !Vtmp(4,:) = soln%Vmms(i1,ind(4),:)
!    !psiPtmp = soln%psi_plus(i1,ind,:,2)
!    !psiMtmp = soln%psi_minus(i1,ind,:,2)
!    !call MUSCL_extrap(Vtmp, psiPtmp, psiMtmp, left, right)
!    left = soln%V(i1,j1-1,:)
!    right = soln%Vmms(i1,j1,:)
!    call flux_fun(left,right,nx,ny,soln%F(i1,j1,:,2))
!  end do
!  
!  !call calc_sources(soln,grid)
!  call calc_time_step(grid%A_xi,grid%A_eta,grid%n_xi_avg, &
!                      grid%n_eta_avg,grid%V,soln%V,soln%dt)
!  
!  call explicit_RK(grid,soln)
!  call cons2prim(soln%U,soln%V)
!  if (k==1) then 
!  call residual_norms(soln%R,soln%Rnorm,2,(/one,one,one,one/))
!  soln%rinit = soln%Rnorm
!  else
!  call residual_norms(soln%R,soln%Rnorm,2,soln%rinit)
!  end if
!  
!  if (mod(k,res_save)==0) then
!  call calc_DE(soln,soln%DE,soln%DEnorm,0,cons)
!  call output_res(soln,k)
!  end if
!  
!  if (mod(k,res_out)==0) then
!  write(*,*) k, soln%Rnorm
!  end if
!  
!  if (mod(k,soln_save)==0) then
!  call calc_DE(soln,soln%DE,soln%DEnorm,0,cons)
!  call output_soln(grid,soln,k)
!  end if
!  if (all(soln%Rnorm<tol)) then
!    exit
!  end if
!
!  end do
!  call calc_DE(soln,soln%DE,soln%DEnorm,0,cons)
!  call output_res(soln,k)
!  call output_soln(grid,soln,k)
!  
!  call grid_out(geometry_file,grid)
!  call teardown_geometry(grid,soln)
!  close(50)
!  write(*,*) 'Program End'
end program main_program
