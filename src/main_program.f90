program main_program
  
  use set_constants, only : set_derived_constants
  use fluid_constants, only : set_fluid_constants
  use set_inputs, only : set_derived_inputs
  use set_inputs, only : max_iter, tol, soln_save, res_save
  use set_inputs, only : leftV, rightV, leftU, rightU, flux_scheme
  use set_inputs, only : limiter_freeze, res_out, cons
  use variable_conversion
  use time_integration
  use boundary_conds, only : enforce_bndry
  use limiter_calc, only : select_limiter
  use flux_calc, only : select_flux, flux_fun
  use other_subroutines
  use geometry, only : setup_geometry, teardown_geometry
  use init_problem, only : initialize
  use namelist, only : read_namelist
  use grid_type
  use soln_type
  implicit none
  
  character(len=100) :: header_str1
  character(len=100) :: niter
  integer :: j, pnorm
  type( grid_t )      :: grid
  type( soln_t )      :: soln
  
  !open(50,file='temp.txt',status='unknown')
  
  100 format(1(I0.8),3(G20.7))
  
  write(header_str1,*) " Iter  |   ||density||    |"// &
  & "   ||velocity||    |   ||pressure||    |"
  call set_derived_constants
  call set_fluid_constants
  call read_namelist('input.nml')
  call output_file_headers
  call set_derived_inputs
  call setup_geometry(grid,soln)
  
  call select_flux()
  call select_limiter()
  
  call initialize(grid,soln)
  
  call enforce_bndry( soln )
  call update_states( soln )
  call calculate_sources(soln%V(:,3),grid%dAc,soln%src)
  call calc_time_step(grid%dx,soln%V,soln%asnd,soln%lambda,soln%dt) 
  
  call MUSCL_extrap( soln%V, leftV, rightV )
  call prim2cons(leftU,leftV)
  call prim2cons(rightU,rightV)
  call flux_fun(leftU,rightU,soln%F)
  
  call explicit_euler(grid,soln%src,soln%dt,soln%F,soln%U,soln%R)
  call update_states( soln )
  
  call residual_norms(soln%R,soln%rinit,pnorm,(/one,one,one/))
  soln%rold = soln%rinit
  
  write(*,*) 'Residual Norms: Iteration 0'
  write(*,*) header_str1
  write(*,100) j, soln%rinit(1), soln%rinit(2), soln%rinit(3)
  write(*,*)
  write(*,*) 'Relative Residual Norms:'
  
  call teardown_geometry(grid,soln)
  close(50)
end program main_program
