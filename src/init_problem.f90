module init_problem
  
  use set_precision,   only : prec
  use set_constants,   only : zero, one, two, half
  use fluid_constants, only : R_gas, gamma
  use set_inputs,      only : p0, T0, eps, ig_low, ig_high
  use soln_type,       only : soln_t
  use grid_type,       only : grid_t
  use variable_conversion, only : prim2cons
  
  implicit none
  
  private
  
  public :: initialize, initialize_MMS
  
  contains
  
  !================================ initialize  ==============================80
  !>
  !! Description: 
  !!
  !! Inputs:      grid : 
  !!              soln : 
  !!
  !! Outputs:     grid :
  !!              soln :
  !<
  !===========================================================================80
  subroutine initialize( grid, soln )
    
    implicit none
    
    integer :: i
    type( grid_t ), intent(inout) :: grid
    type( soln_t ), intent(inout) :: soln
    
    

  end subroutine initialize
  
  subroutine initialize_MMS( grid, soln )
    
    use soln_type, only : calc_MMS
    type( grid_t ), intent(inout) :: grid
    type( soln_t ), intent(inout) :: soln
    call calc_mms(grid,soln)
    call prim2cons(soln%Umms,soln%Vmms)
    soln%V = soln%Vmms
    soln%U = soln%Umms
    soln%S = soln%Smms
    
  end subroutine initialize_MMS
end module init_problem
