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
  
  public :: initialize, initialize_MMS, initialize_const, initialize_lin_xi
  
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
  
  subroutine initialize_lin_xi( grid, soln, V1, V2 )
    
    use other_subroutines, only : calc_sources
    use variable_conversion, only : speed_of_sound
    use set_inputs, only : ig_low, ig_high, jg_low, jg_high, &
                           imax
    type( grid_t ), intent(inout) :: grid
    type( soln_t ), intent(inout) :: soln
    real(prec), dimension(4), intent(in) :: V1, V2
    integer :: i, j
    real(prec) :: Ni, Nj, asnd, mach, p0, temp
    Ni = real(imax-2,prec)
    call speed_of_sound( V1(4), V1(1), asnd )
    mach = sqrt(V1(2)**2 + V1(3)**2)/asnd
    p0 = V1(4)*(one+half*(gamma-one)*mach**2)**(gamma/(gamma-one))
    temp = asnd**2/(gamma*R_gas)
    do j = jg_low, jg_high
      do i = ig_low, ig_high
        !soln%V(2,i,j) = V1(2) + real(i-1,prec)/Ni*(V2(2)-V1(2))
        !soln%V(3,i,j) = V1(3) + real(i-1,prec)/Ni*(V2(3)-V1(3))
        !mach = sqrt(soln%V(2,i,j)**2 + soln%V(3,i,j)**2)/asnd
        !soln%V(4,i,j) = p0/&
        !     ((one+half*(gamma-one)*mach**2)**(gamma/(gamma-one)))
        soln%V(1,i,j) = soln%V(4,i,j)/(R_gas*temp)
        soln%V(1,i,j) = V1(1) + real(i-1,prec)/Ni*(V2(1)-V1(1))
        soln%V(2,i,j) = V1(2) + real(i-1,prec)/Ni*(V2(2)-V1(2))
        soln%V(3,i,j) = V1(3) + real(i-1,prec)/Ni*(V2(3)-V1(3))
        soln%V(4,i,j) = V1(4) + real(i-1,prec)/Ni*(V2(4)-V1(4))
      end do
    end do
    call calc_sources( soln, grid )
    call prim2cons(soln%U,soln%V)
    
  end subroutine initialize_lin_xi
  
end module init_problem
