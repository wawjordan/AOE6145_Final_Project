module init_problem
  
  use set_precision,   only : prec
  use set_constants,   only : zero, one, two, half
  use fluid_constants, only : R_gas, gamma
  use set_inputs,      only : p_inf, T_inf, eps, ig_low, ig_high
  use soln_type,       only : soln_t
  use grid_type,       only : grid_t
  use variable_conversion, only : prim2cons
  
  implicit none
  
  private
  
  public :: initialize, initialize_MMS, initialize_const
  
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
    integer :: i,j
    call calc_MMS(grid,soln)
    call prim2cons(soln%Umms,soln%Vmms)
    soln%V = soln%Vmms
    soln%U = soln%Umms
    soln%S = soln%Smms
    !do j = grid%j_low,grid%jg_high
    !do i = grid%i_low,grid%ig_high
    !soln%V(:,i,j) = soln%V(:,i-1,j-1)
    !soln%U(:,i,j) = soln%U(:,i-1,j-1)
    !end do
    !end do
    
  end subroutine initialize_MMS
  
  subroutine initialize_const( grid, soln, V )
    use other_subroutines, only : calc_sources
    type( grid_t ), intent(inout) :: grid
    type( soln_t ), intent(inout) :: soln
    real(prec), dimension(4), intent(in) :: V
    soln%V(1,:,:) = V(1)
    soln%V(2,:,:) = V(2)
    soln%V(3,:,:) = V(3)
    soln%V(4,:,:) = V(4)
    call calc_sources( soln, grid )
    call prim2cons(soln%U,soln%V)
    
  end subroutine initialize_const

  
end module init_problem
