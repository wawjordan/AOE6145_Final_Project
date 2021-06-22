program main_program
  
  use set_precision, only : prec  
  use set_constants, only : zero, set_derived_constants
  use fluid_constants, only : set_fluid_constants
  use set_inputs, only : set_derived_inputs, out_dir, out_file, isMMS, cons, &
                         flux_out, limiter_freeze, C_grid, inlet, res_save,  &
                         res_out, soln_save, max_iter, tol, epsM, neq, imax, &
                         rho_inf, p_inf, u0, v0!,CFL
  use file_handling, only : grid_out, output_file_headers, output_flux, &
                            output_exact_soln, output_soln, output_res, &
                            output_airfoil_forces, output_pressure_loss
  use geometry, only : setup_geometry, teardown_geometry
  use variable_conversion, only : prim2cons, cons2prim, update_states, &
                                  limit_primitives
  use other_subroutines, only : calc_de, calc_sources, &
                                airfoil_forces, set_boundary_conditions
  use time_integration, only : calc_time_step, explicit_RK, residual_norms
  use limiter_calc, only : select_limiter, limiter_fun, calculate_limiters
!  use mms_functions, only : rho_mms, uvel_mms, vvel_mms, press_mms
  use init_problem, only : initialize_MMS, initialize_const
  use namelist, only : read_namelist
  use grid_type, only : grid_t
  use soln_type, only : soln_t
  use bc_type, only : bc_t
  use flux_calc, only : select_flux, calc_flux_2D
  implicit none
  
  type( grid_t )      :: grid
  type( soln_t )      :: soln
  type(bc_t), dimension(:), allocatable :: bnds
  real(prec), dimension(4) :: Rnorm2
  character(len=72) :: header_str1
  integer :: i, k, counter
  
  100 format (I6, T10, ES12.6, T26, ES12.6, T42, ES12.6, T58, ES12.6)
  write(header_str1,*) "Iter "//&
                   & "|     mass     "// &
                   & "|   x-momentum  "//&
                   & "|   y-momentum  "//&
                   & "|      energy    |"
   
  call set_derived_constants
  call set_fluid_constants
  call read_namelist('input.nml')
  call set_derived_inputs
  call setup_geometry(grid,soln)
  call select_flux()
  call select_limiter()
  if (isMMS) then
    call initialize_MMS(grid,soln)
    call update_states(soln)
  else
    call initialize_const(grid,soln,&
         (/rho_inf,u0,v0,p_inf/) )
    call update_states(soln)
  end if
  call output_file_headers(out_dir,out_file)
  if (isMMS) then
    call output_exact_soln(grid,soln,out_dir,out_file)
  end if
  call output_soln(grid,soln,0)
  call set_boundary_conditions(grid,bnds)
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
  
  open(42,file='temp.txt',status='unknown')
  
  counter = 0
  
  write(*,*) header_str1
  do k = 1,max_iter
!!==============================================================================
     
!  call limit_primitives(soln%V)
  call prim2cons(soln%U,soln%V)
  call update_states(soln)
  if ((limiter_freeze .eqv. .false.).and.(epsM>zero)) then
    call calculate_limiters(soln)
  end if
    
    call calc_time_step(grid,soln)
    call explicit_RK(grid,soln,bnds)
    call cons2prim(soln%U,soln%V)
    if (k==2) then
      call residual_norms(soln%R,soln%Rnorm,soln%rinit)
      soln%rinit = soln%Rnorm
      write(*,100) k, soln%Rnorm
    else
      call residual_norms(soln%R,soln%Rnorm,soln%rinit)
    end if
    
    !if (all(soln%Rnorm<1.0e-5_prec)) then
      !call SER(CFL,2.0_prec,soln%Rnorm,Rnorm2)
      !call RDM(CFL,2.0_prec,1.1_prec,soln%Rnorm,Rnorm2)
    !end if
    
    if (mod(k,res_save)==0) then
      if (isMMS) then
        call calc_DE(soln,soln%DE,soln%DEnorm,cons)
      end if
      call output_res(grid,soln,k)
    end if
    
    if (mod(k,res_out)==0) then
      write(*,100) k, soln%Rnorm
    end if
    
    if (mod(k,soln_save)==0) then
      if (C_grid) then
        call output_airfoil_forces(grid,soln,k)
      end if
      if (inlet) then
        call output_pressure_loss(grid,soln,k)
      end if
      if (isMMS) then
        call calc_DE(soln,soln%DE,soln%DEnorm,cons)
      end if
      call output_soln(grid,soln,k)
      if (flux_out) then
        call output_flux(grid,soln,k)
      end if
    end if
    
    if (all(soln%Rnorm<1.0e-3_prec).or.(k>5000)) then
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
      write(*,100) k, soln%Rnorm
      write(*,*) 'solution converged'
      exit
    end if
    Rnorm2 = soln%Rnorm
    
  end do
  
  write(42,*) imax-1, 5, (soln%Rnorm(i),i=1,neq)
  if (C_grid) then
    call output_airfoil_forces(grid,soln,k)
  end if
  if (inlet) then
    call output_pressure_loss(grid,soln,k)
  end if
  if (isMMS) then 
    call calc_DE(soln,soln%DE,soln%DEnorm,cons)
    write(42,*) imax-1, 1,(soln%DEnorm(i,1),i=1,neq)
    write(42,*) imax-1, 2,(soln%DEnorm(i,2),i=1,neq)
    write(42,*) imax-1, 0,(soln%DEnorm(i,3),i=1,neq)
  end if
  close(42)
  call output_res(grid,soln,k)
  call output_soln(grid,soln,k)
  
  call teardown_geometry(grid,soln)
  write(*,*) 'Program End'
end program main_program
  
