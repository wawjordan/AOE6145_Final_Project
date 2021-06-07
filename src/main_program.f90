program main_program
  
  use set_precision, only : prec  
  use set_constants, only : zero, one, set_derived_constants
  use fluid_constants, only : set_fluid_constants
  use set_inputs, only : set_derived_inputs, geometry_file, limiter_freeze
  use set_inputs, only : i_low,i_high,j_low,j_high, n_ghost, cons, neq
  use set_inputs, only : ig_low,ig_high,jg_low,jg_high
  use set_inputs, only : res_save, res_out, soln_save, tol, Lmms, CFL
  use file_handling, only : grid_out, output_file_headers, &
                            output_exact_soln, output_soln, output_res
  use geometry, only : setup_geometry, teardown_geometry
  use variable_conversion, only : prim2cons, cons2prim, update_states, &
                                  limit_primitives
  use other_subroutines, only : calc_de, MUSCL_extrap
  use time_integration, only : calc_time_step, explicit_RK, residual_norms
  use limiter_calc, only : select_limiter, limiter_fun, calculate_limiters
  use mms_functions, only : rho_mms, uvel_mms, vvel_mms, press_mms
  use init_problem, only : initialize_MMS
  !use namelist, only : read_namelist
  use grid_type, only : grid_t
  use soln_type, only : soln_t, calc_mms
  use flux_calc, only : select_flux, calc_flux_2D, flux_fun
  implicit none
  
  integer :: i, j, k
  type( grid_t )      :: grid
  type( soln_t )      :: soln
  real(prec), dimension(4) :: left, right
  real(prec) :: nx, ny
  real(prec), dimension(4) :: Rnorm2
  integer :: i1,j1,k1
  logical :: freeze_BC = .false.
!  integer, dimension(4) :: ind
!  real(prec), dimension(4,4) :: Vtmp, psiPtmp, psiMtmp


  
  call set_derived_constants
  call set_fluid_constants
  call select_flux()
  call select_limiter()
  call set_derived_inputs
  call setup_geometry(grid,soln)
  call initialize_MMS(grid,soln)
  call output_exact_soln(grid,soln)
  call output_file_headers  
  call output_soln(grid,soln,0)
!  !call calc_sources(soln,grid)
  
 
 
  do k = 1,100000
!!==============================================================================
  !  soln%V(:,i_low:i_high,jg_high-1) = 2*soln%V(:,i_low:i_high,j_high  ) &
  !                                     - soln%V(:,i_low:i_high,j_high-1)
  !  soln%V(:,i_low:i_high,jg_high)   = 2*soln%V(:,i_low:i_high,j_high-1) &
  !                                     - soln%V(:,i_low:i_high,j_high-2)
  !  soln%V(:,ig_high-1,j_low:j_high) = 2*soln%V(:,i_high  ,j_low:j_high) &
  !                                     - soln%V(:,i_high-1,j_low:j_high)
  !  soln%V(:,ig_high,j_low:j_high)   = 2*soln%V(:,i_high-1,j_low:j_high) &
  !                                     - soln%V(:,i_high-2,j_low:j_high)
     
    call limit_primitives(soln%V)
    call prim2cons(soln%U,soln%V)
    call update_states(soln)
    !if (limiter_freeze .eqv. .false.) then
    !  call calculate_limiters(soln)
    !end if
    call calc_flux_2D(grid,soln)
  
  !i1 = i_high
  !do j1 = j_low,j_high
  !  nx = grid%n_xi(i1,j1,1)
  !  ny = grid%n_xi(i1,j1,2)
  !  !ind = (/ ( k1,k1=i1-2,i1+1 ) /)
  !  !Vtmp(1,:) = soln%V(ind(1),j1,:)
  !  !Vtmp(2,:) = soln%Vmms(ind(2),j1,:)
  !  !Vtmp(3,:) = soln%Vmms(ind(3),j1,:)
  !  !Vtmp(4,:) = soln%Vmms(ind(4),j1,:)
  !  !psiPtmp = soln%psi_plus(ind,j1,:,1)
  !  !psiMtmp = soln%psi_minus(ind,j1,:,1)
  !  !call MUSCL_extrap(Vtmp, psiPtmp, psiMtmp, left, right)
  !  left = soln%V(:,i1,j1)
  !  right = soln%Vmms(:,i1+1,j1)
  !  call flux_fun(left,right,nx,ny,soln%Fxi(:,i1,j1))
  !end do
  !j1 = j_high
  !do i1 = i_low,i_high
  !  nx = grid%n_eta(i1,j1,1)
  !  ny = grid%n_eta(i1,j1,2)
  !  !ind = (/ ( k1,k1=j1-2,j1+1 ) /)
  !  !Vtmp(1,:) = soln%V(i1,ind(1),:)
  !  !Vtmp(2,:) = soln%V(i1,ind(2),:)
  !  !Vtmp(3,:) = soln%V(i1,ind(3),:)
  !  !Vtmp(4,:) = soln%Vmms(i1,ind(4),:)
  !  !psiPtmp = soln%psi_plus(i1,ind,:,2)
  !  !psiMtmp = soln%psi_minus(i1,ind,:,2)
  !  !call MUSCL_extrap(Vtmp, psiPtmp, psiMtmp, left, right)
  !  left = soln%V(:,i1,j1)
  !  right = soln%Vmms(:,i1,j1+1)
  !  call flux_fun(left,right,nx,ny,soln%Feta(:,i1,j1))
  !end do

    call calc_time_step(grid,soln)
    call explicit_RK(grid,soln)
    call cons2prim(soln%U,soln%V)
    if (k==1) then
      call residual_norms(soln%R,soln%Rnorm,2,soln%rinit)
      soln%rinit = soln%Rnorm
      Rnorm2 = soln%Rnorm
    else
      Rnorm2 = soln%Rnorm
      call residual_norms(soln%R,soln%Rnorm,2,soln%rinit)
    end if
    
    if (mod(k,res_save)==0) then
      call calc_DE(soln,soln%DE,soln%DEnorm,cons)
      call output_res(soln,k)
    end if
    
    if (mod(k,res_out)==0) then
      write(*,*) k, soln%Rnorm
    end if
    
    if (mod(k,soln_save)==0) then
      call calc_DE(soln,soln%DE,soln%DEnorm,cons)
      call output_soln(grid,soln,k)
    end if
    
    if (all(soln%Rnorm<1e-2)) then
      limiter_freeze = .true.
    end if
   ! if (all(soln%Rnorm<1e-4).and.(all(soln%Rnorm<Rnorm2))) then
   !   if (CFL<0.5_prec) then
   !     CFL = CFL+1.0e-3_prec
   !   end if
   ! end if
    
    if (all(soln%Rnorm<tol)) then
      exit
    end if
    
  end do
  call calc_DE(soln,soln%DE,soln%DEnorm,cons)
  call output_res(soln,k)
  call output_soln(grid,soln,k)
  
  call grid_out(geometry_file,grid)
  call teardown_geometry(grid,soln)
  close(50)
  write(*,*) 'Program End'
end program main_program
  
