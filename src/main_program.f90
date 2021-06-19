program main_program
  
  use set_precision, only : prec  
  use set_constants, only : pi, zero, one, two, set_derived_constants
  use fluid_constants, only : set_fluid_constants, R_gas
  use set_inputs, only : set_derived_inputs, geometry_file, limiter_freeze, epsM
  use set_inputs, only : imax, i_low,i_high,j_low,j_high, n_ghost, cons, neq
  use set_inputs, only : jmax, ig_low,ig_high,jg_low,jg_high, max_iter
  use set_inputs, only : C_grid, index1, index2
  use set_inputs, only : res_save, res_out, soln_save, tol, Lmms, CFL, isMMS
  use set_inputs, only : alpha, rho_inf, u_inf, p_inf, u0, v0, num_BCs, bounds
  use file_handling, only : grid_out, output_file_headers, output_flux, &
                            output_exact_soln, output_soln, output_res, &
                            output_airfoil_forces
  use geometry, only : setup_geometry, teardown_geometry
  use variable_conversion, only : prim2cons, cons2prim, update_states, &
                                  limit_primitives
  use other_subroutines, only : calc_de, MUSCL_extrap, calc_sources, &
                                surface_MUSCL_extrap, airfoil_forces
  use time_integration, only : calc_time_step, explicit_RK, residual_norms, &
                               calc_residual, SER, RDM
  use limiter_calc, only : select_limiter, limiter_fun, calculate_limiters
  use mms_functions, only : rho_mms, uvel_mms, vvel_mms, press_mms
  use init_problem, only : initialize_MMS, initialize_const, initialize_lin_xi
  use namelist, only : read_namelist
  use grid_type, only : grid_t
  use soln_type, only : soln_t, calc_mms
  use bc_type, only : bc_t, reflect_vec
  use slip_wall_bc_type, only : eta_wall_t
  use dir_bc_type, only : dir_bc_t
  use wake_cut_bc_type, only : xi_cut_t
  use flux_calc, only : select_flux, calc_flux_2D, flux_fun, exact_flux
  implicit none
  
  integer :: i, j, k
  type( grid_t )      :: grid
  type( soln_t )      :: soln
  type(bc_t), dimension(:), allocatable :: bnds
!  type(bc_t) :: Lbnd, Rbnd, Bbnd, Tbnd
 ! type(bc_t)   :: inlet, outlet
 ! type(eta_wall_t) :: top_wall, bottom_wall
!  type(bc_t) :: rear1, rear2, airfoil, farfield, wake
 ! type(dir_bc_t) :: farfield!, rear1, rear2
!  type(eta_wall_t) :: airfoil
!  type(bc_t) :: airfoil
 ! real(prec), allocatable, dimension(:,:) :: Cp
  real(prec) :: nx, ny, press, Cl, Cd
  real(prec), dimension(4) :: Rnorm2, left, right, flux
  integer :: i1,j1,k1, counter
  
  call set_derived_constants
  call set_fluid_constants
  call read_namelist('input.nml')
  call set_derived_inputs
  call setup_geometry(grid,soln)
  call select_flux()
  call select_limiter()
!  call initialize_MMS(grid,soln)
  call initialize_const(grid,soln,&
       (/rho_inf,u0,v0,p_inf/) )
!  call initialize_lin_xi(grid,soln,          &
!       (/one,10.0_prec,zero,100000.0_prec/), &
!       (/one,20.0_prec,zero,100000.0_prec/)  )
  call update_states(soln)
  if (isMMS) then
    call output_exact_soln(grid,soln)
  end if
  call output_file_headers
  call output_soln(grid,soln,0)
! For Cartesian Mesh
!  call Lbnd%set_bc(grid,1,(/1,0/),ig_low,ig_low+1,j_low,j_high)
!  call Rbnd%set_bc(grid,4,(/-1,0/),ig_high-1,ig_high,j_low,j_high)
!  call Bbnd%set_bc(grid,1,(/0,1/),i_low,i_high,jg_low,jg_low+1)
!  call Tbnd%set_bc(grid,4,(/0,-1/),i_low,i_high,jg_high-1,jg_high)

!  call Lbnd%set_bc(grid,1,(/1,0/),ig_low,ig_low+1,j_low,j_high)
!  call Rbnd%set_bc(grid,1,(/-1,0/),ig_high-1,ig_high,j_low,j_high)
!  call Bbnd%set_bc(grid,1,(/0,1/),i_low,i_high,jg_low,jg_low+1)
!  call Tbnd%set_bc(grid,1,(/0,-1/),i_low,i_high,jg_high-1,jg_high)

! For Curvilinear Mesh...
!  call Lbnd%set_bc(grid,1,(/-1,0/),ig_high-1,ig_high,j_low,j_high)
!  call Rbnd%set_bc(grid,4,(/1,0/),ig_low,ig_low+1,j_low,j_high)
!  call Bbnd%set_bc(grid,1,(/0,-1/),i_low,i_high,jg_high-1,jg_high)
!  call Tbnd%set_bc(grid,4,(/0,1/),i_low,i_high,jg_low,jg_low+1)


! For Inlet Mesh
!  call inlet%set_bc(grid,2,(/0,1/),i_low,20,jg_low,jg_low+1)
!  call outlet%set_bc(grid,4,(/0,-1/),ig_high-1,ig_high,j_low,j_high)
!  call top_wall%set_bc(grid,5,21,i_high,j_high,j_high)
!  call bottom_wall%set_bc(grid,5,i_low,i_high,j_low,j_low)
! For Airfoil Mesh

!  call farfield%set_bc(grid,2,i_low,i_high,jg_high-1,jg_high)
!  call rear1%set_bc(grid,3,(/2,1/),ig_high-1,ig_high,j_low,j_high)
!  call rear2%set_bc(grid,3,(/2,-1/),ig_low,ig_low+1,j_low,j_high)
!  call farfield%set_bc(grid,2,i_low,i_high,jg_high-1,jg_high,4)

!  call rear1%set_bc(grid,3,ig_high-1,ig_high,j_low,j_high,2)
!  call rear2%set_bc(grid,3,ig_low,ig_low+1,j_low,j_high,1)
!  call airfoil%set_bc(grid,5,index2+1,imax-index2-1,jg_low,j_low-1,3)
!  call wake%set_bc(grid,6,ig_low,index2,jg_low,j_low-1,3)
  allocate(bnds(num_BCs))
  do i = 1,num_BCs
    call bnds(i)%set_bc(grid,bounds(1,i), bounds(2,i), bounds(3,i), &
                             bounds(4,i), bounds(5,i), bounds(6,i) )
  end do
!  call bnds(1)%set_bc(grid,2,i_low,i_high,jg_high-1,jg_high,4)
!  call bnds(2)%set_bc(grid,3,ig_high-1,ig_high,j_low,j_high,2)
!  call bnds(3)%set_bc(grid,3,ig_low,ig_low+1,j_low,j_high,1)
!  call bnds(4)%set_bc(grid,5,index2+1,imax-index2-1,jg_low,j_low-1,3)
!  call bnds(5)%set_bc(grid,6,ig_low,index2,jg_low,j_low-1,3)
  
!  write(*,*) 2, i_low, i_high, jg_high-1, jg_high,4
!  write(*,*) 3, ig_high-1, ig_high, j_low, j_high,2
!  write(*,*) 3, ig_low, ig_low+1, j_low, j_high,1
!  write(*,*) 5, index2+1, imax-index2-1, jg_low, j_low-1,3
!  write(*,*) 6, ig_low, index2, jg_low, j_low-1,3
!  write(*,*) imax, jmax
!  stop
!  call wake2%set_bc(grid,10,ig_low,index2,n_ghost)
!  stop
!  call wake%set_bc(grid,10,ig_low,index2,n_ghost)
!  call rear1%set_bc(grid,2,ig_high-1,ig_high,j_low,j_high)
!  call rear2%set_bc(grid,2,ig_low,ig_low+1,j_low,j_high)

!      do j = j_low-1,jg_low,-1
!        do i = index2+1, imax-index2-1
!          call reflect_vec( soln%V(2,i,j+1), soln%V(3,i,j+1), &
!                              grid%n_eta(i,j_low,1), grid%n_eta(i,j_low,2), &
!                              soln%V(2,i,j), soln%V(3,i,j) )
!          soln%V(4,i,j) = soln%V(4,i,j+1) + zero*epsM*( &
!                            soln%V(4,i,j+1) - soln%V(4,i,j+2) )
!          soln%V(1,i,j)   = soln%V(4,i,j)/( R_gas*soln%temp(i,j+1) )
!        end do
!      end do
!  call airfoil%set_bc(grid,5,index2+1,imax-index2-1,j_low,j_low)
!  call airfoil%set_bc(grid,5,(/1,0/),index2+1,imax-index2-1,jg_low,jg_low+1)
!  do j = 1,wake%N
!    do i = 1,wake%M
!      write(*,*) wake%i1(i,j), wake%j1(i,j), wake%i2(i,j), wake%j2(i,j)
!    end do
!  end do
!  stop
!  call wall%set_bc(grid,5,i_low,i_high,j_low,j_low)
  
  open(42,file='temp.txt',status='unknown')
!  call farfield%set_val(soln)
!  call rear1%set_val(soln)
!  call rear2%set_val(soln)
  counter = 0
  do k = 1,max_iter
!!==============================================================================
!  call Lbnd%set_val(soln)
!  call Rbnd%set_val(soln)
!  call Bbnd%set_val(soln)
!  call Tbnd%set_val(soln)
!
!  call Lbnd%enforce(soln)
!  call Rbnd%enforce(soln)
!  call Bbnd%enforce(soln)
!  call Tbnd%enforce(soln)
!  call rear1%set_val(soln)
!  call rear2%set_val(soln)
!  call airfoil%set_val(soln)
!  call wake%set_val(soln)
!  call farfield%set_val(soln)
  call bnds(1)%set_val(soln)
  call bnds(2)%set_val(soln)
  call bnds(3)%set_val(soln)
  call bnds(4)%set_val(soln)
  call bnds(5)%set_val(soln)
  !call farfield%enforce(soln)
  !call rear1%enforce(soln)
  !call rear2%enforce(soln)
!  call wake%enforce(soln)
!  call airfoil%enforce(soln)
  
  
!      do j = j_low-1,jg_low,-1
!        do i = index2+1, imax-index2-1
!          call reflect_vec( soln%V(2,i,j+1), soln%V(3,i,j+1), &
!                              grid%n_eta(i,j_low,1), grid%n_eta(i,j_low,2), &
!                              soln%V(2,i,j), soln%V(3,i,j) )
!          soln%V(4,i,j) = soln%V(4,i,j+1) + epsM*( &
!                            soln%V(4,i,j+1) - soln%V(4,i,j+2) )
!          soln%V(1,i,j)   = soln%V(4,i,j)/( R_gas*soln%temp(i,j+1) )
!        end do
!      end do
  
!  call wall%set_val(soln)
!  call inlet%set_val(soln)
!  call outlet%set_val(soln)
!  call inlet%enforce(soln)
!  call outlet%enforce(soln)
     
!  call limit_primitives(soln%V)
  call prim2cons(soln%U,soln%V)
  call update_states(soln)
  if ((limiter_freeze .eqv. .false.).and.(epsM>zero)) then
    call calculate_limiters(soln)
  end if
  call calc_flux_2D(grid,soln)
  
!  do j = j_low, j_high
!  do i = i_low-1, i_high
!    left = 0.5_prec*(soln%V(:,i,j) + soln%V(:,i+1,j))
!    call exact_flux(left,grid%n_xi(i,j,1),grid%n_xi(i,j,2),flux)
!    write(*,*) i, j
!    write(*,*) 'V   = ',left
!    write(*,*) 'F   = ',soln%Fxi(:,i,j)
!    write(*,*) 'Fex = ',flux
!    write(*,*) '!============================================================='
!  end do
!  end do
!  call output_flux(grid,soln,k)
!  stop  
  
  call bnds(4)%set_val(soln)
  call bnds(5)%set_val(soln)
!  call airfoil%set_val(soln)
!  call wake%set_val(soln)
!   call wall%enforce(soln)
   
!  call top_wall%set_val(soln)
!  call bottom_wall%set_val(soln)
!  call top_wall%enforce(soln)
!  call bottom_wall%enforce(soln)

!  j1 = j_high
!  do i1 = 21,i_high
!    press = soln%V(4,i1,j1) + epsM*(soln%V(4,i1,j1) - soln%V(4,i1,j1-1))
!    nx = grid%n_eta(i1,j1+1,1)
!    ny = grid%n_eta(i1,j1+1,2)
!    !soln%Feta(:,i1,j1) = (/ zero, nx*soln%V(4,i1,j1), ny*soln%V(4,i1,j1),zero/)
!    soln%Feta(:,i1,j1) = (/ zero, nx*press, ny*press,zero/)
!  end do
!  j1 = j_low
!  do i1 = 17,80
!    press = soln%V(4,i1,j1) + epsM*(soln%V(4,i1,j1) - soln%V(4,i1,j1+1))
!    nx = grid%n_eta(i1,j1,1)
!    ny = grid%n_eta(i1,j1,2)
!    !soln%Feta(:,i1,j1-1) = (/ zero, nx*soln%V(4,i1,j1), ny*soln%V(4,i1,j1),zero/)
!    soln%Feta(:,i1,j1-1) = (/ zero, nx*press, ny*press,zero/)
!  end do
!  i1 = i_low
!  do j1 = j_low,j_high
!    nx = grid%n_xi(i1,j1,1)
!    ny = grid%n_xi(i1,j1,2)
!    left = soln%Vmms(:,i1,j1)
!    !right(1) = rho_mms(Lmms,grid%x(i1+1,j1),grid%y(i1+1,j1))
!    !right(2) = uvel_mms(Lmms,grid%x(i1+1,j1),grid%y(i1+1,j1))
!    !right(3) = vvel_mms(Lmms,grid%x(i1+1,j1),grid%y(i1+1,j1))
!    !right(4) = press_mms(Lmms,grid%x(i1+1,j1),grid%y(i1+1,j1))
!    right = soln%V(:,i1-1,j1)
!    !call exact_flux(left,nx,ny,soln%Fxi(:,i1,j1))
!    call flux_fun(left,right,nx,ny,soln%Fxi(:,i1-1,j1))
!  end do
!  j1 = j_low
!  do i1 = i_low,i_high
!    nx = grid%n_eta(i1,j1,1)
!    ny = grid%n_eta(i1,j1,2)
!    left = soln%Vmms(:,i1,j1)
!    !right(1) = rho_mms(Lmms,grid%x(i1,j1+1),grid%y(i1,j1+1))
!    !right(2) = uvel_mms(Lmms,grid%x(i1,j1+1),grid%y(i1,j1+1))
!    !right(3) = vvel_mms(Lmms,grid%x(i1,j1+1),grid%y(i1,j1+1))
!    !right(4) = press_mms(Lmms,grid%x(i1,j1+1),grid%y(i1,j1+1))
!    right = soln%V(:,i1,j1-1)
!    !call exact_flux(left,nx,ny,soln%Feta(:,i1,j1))
!    call flux_fun(left,right,nx,ny,soln%Feta(:,i1,j1-1))
!  end do
    
    call calc_time_step(grid,soln)
    call explicit_RK(grid,soln)
    call cons2prim(soln%U,soln%V)
    if (k==2) then
      call residual_norms(soln%R,soln%Rnorm,2,soln%rinit)
      soln%rinit = soln%Rnorm
    else
      call residual_norms(soln%R,soln%Rnorm,2,soln%rinit)
    end if
    
!    if (all(soln%Rnorm<1.0e-2_prec)) then
!      !call SER(CFL,4.0_prec,soln%Rnorm,Rnorm2)
!      call RDM(CFL,4.0_prec,2.0_prec,soln%Rnorm,Rnorm2)
!    end if
    
    if (mod(k,res_save)==0) then
      if (isMMS) then
        call calc_DE(soln,soln%DE,soln%DEnorm,cons)
      end if
      call output_res(soln,k)
    end if
    
    if (mod(k,res_out)==0) then
      write(*,*) k, soln%Rnorm
    end if
    
    if (mod(k,soln_save)==0) then
      if (C_grid) then
        call output_airfoil_forces(grid,soln,k)
      end if
      if (isMMS) then
        call calc_DE(soln,soln%DE,soln%DEnorm,cons)
      end if
      call output_soln(grid,soln,k)
    !  call output_flux(grid,soln,k)
    end if
    
    if (all(soln%Rnorm<1.0e-3_prec).or.(k>20000)) then
      limiter_freeze = .true.
    end if
    if (all(abs(soln%Rnorm-Rnorm2)<1.0e-16_prec)) then
      counter = counter + 1
    else
      counter = 0
    end if
    if (counter>1000) then
      write(*,*) k, soln%Rnorm
      write(*,*) 'solution stalled'
      exit
    end if
    if (all(soln%Rnorm<tol)) then
      write(*,*) k, soln%Rnorm
      write(*,*) 'solution converged'
      exit
    end if
    Rnorm2 = soln%Rnorm
    
  end do
  
  write(42,*) i_high, 5, (soln%Rnorm(i),i=1,neq)
  if (C_grid) then
    call output_airfoil_forces(grid,soln,k)
  end if
  if (isMMS) then 
    call calc_DE(soln,soln%DE,soln%DEnorm,cons)
    write(42,*) i_high, 1,(soln%DEnorm(i,1),i=1,neq)
    write(42,*) i_high, 2,(soln%DEnorm(i,2),i=1,neq)
    write(42,*) i_high, 0,(soln%DEnorm(i,3),i=1,neq)
  end if
  close(42)
  call output_res(soln,k)
  call output_soln(grid,soln,k)
  
  call teardown_geometry(grid,soln)
  write(*,*) 'Program End'
end program main_program
  
